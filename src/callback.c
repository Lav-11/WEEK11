#include "callback.h"
#include "cpx_utils.h"
#include "tsp_utils.h"
#include "heuristics.h"

// Funzione per salvare la soluzione per Gnuplot
void save_solution_for_gnuplot(const char *filename, double *heuristic_sol, instance *inst) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Unable to open file %s for writing.\n", filename);
        return;
    }
    // Scrivi gli archi della soluzione
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            if (heuristic_sol[xpos(i, j, inst)] > 0.5) {
                fprintf(file, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
                fprintf(file, "%f %f\n\n", inst->xcoord[j], inst->ycoord[j]); // Spezzone per Gnuplot
            }
        }
    }
    fclose(file);
    printf("Solution saved for Gnuplot in %s\n", filename);
}

// Callback function
int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    instance *inst = (instance *)userhandle;
    int ncols = inst->ncols;
    double *xstar = (double *)malloc(ncols * sizeof(double));
    if (!xstar) print_error("Memory allocation error");

    double objval = CPX_INFBOUND;
    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        if (CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval)) {
            free(xstar);
            print_error("CPXcallbackgetcandidatepoint error");
        }
    } else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
        if (CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols - 1, &objval)) {
            free(xstar);
            print_error("CPXcallbackgetrelaxationpoint error");
        }
    } else {
        free(xstar);
        return 0; // Non interessa questo contesto
    }

    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    if (!succ || !comp) print_error("Memory allocation error for succ/comp");

    int ncomp = 0;
    build_sol(xstar, inst, succ, comp, &ncomp);

    if (ncomp > 1 && contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        printf("Subtours detected (%d components), trying to post a heuristic solution...\n", ncomp);

        // Heuristic repair: nearest neighbor e VNS
        nearest_neighbor(inst, 0, true);
        double time_limit = 5;
        double learning_rate = 0.01;
        int max_jumps = 5;
        variable_neighborhood_search(inst->best_sol, time_limit, inst, learning_rate, max_jumps);

        // Stampa il tour generato
        printf("Generated tour:\n");
        for (int i = 0; i < inst->nnodes; i++) {
            printf("%d -> ", (int)(inst->best_sol->tour[i]));
        }
        printf("%d\n", (int)(inst->best_sol->tour[0])); // Chiude il tour

        // Prepara la soluzione euristica
        double *heuristic_sol = (double *)calloc(ncols, sizeof(double));
        if (!heuristic_sol) print_error("Memory allocation error for heuristic_sol");

        // Collega i nodi consecutivi del tour
        for (int i = 0; i < inst->nnodes - 1; i++) {
            int node_i = (int)(inst->best_sol->tour[i]) - 1;
            int node_j = (int)(inst->best_sol->tour[i + 1]) - 1;
            int index = xpos(node_i, node_j, inst);
            heuristic_sol[index] = 1.0;
        }
        heuristic_sol[xpos((int)(inst->best_sol->tour[inst->nnodes - 1]) - 1, (int)(inst->best_sol->tour[0]) - 1, inst)] = 1.0;

        // Calcola il valore obiettivo della soluzione euristica
        double heuristic_objval = 0.0;
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (heuristic_sol[xpos(i, j, inst)] > 0.5)
                    heuristic_objval += dist(i, j, inst);
            }
        }

        // Verifica i gradi dei nodi
        int *degree = (int *)calloc(inst->nnodes, sizeof(int));
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (heuristic_sol[xpos(i, j, inst)] > 0.5) {
                    degree[i]++;
                    degree[j]++;
                }
            }
        }
        // printf("Node degrees (debug):\n");
        // for (int i = 0; i < inst->nnodes; i++) {
        //     printf("Node %d: degree = %d\n", i + 1, degree[i]);
        // }
        // int valid = 1;
        // for (int i = 0; i < inst->nnodes; i++) {
        //     if (degree[i] != 2) {
        //         valid = 0;
        //         printf("Node %d has invalid degree %d\n", i + 1, degree[i]);
        //     }
        // }
        free(degree);

        // if (!valid) {
        //     printf("Heuristic solution is invalid (node degree issue).\n");
        //     printf("Edges in the heuristic solution:\n");
        //     for (int i = 0; i < inst->nnodes; i++) {
        //         for (int j = i + 1; j < inst->nnodes; j++) {
        //             if (heuristic_sol[xpos(i, j, inst)] > 0.5) {
        //                 printf("Edge (%d, %d)\n", i + 1, j + 1);
        //             }
        //         }
        //     }
        //     save_solution_for_gnuplot("../data/solution_gnuplot.dat", heuristic_sol, inst);
        //     free(heuristic_sol);
        //     free(succ);
        //     free(comp);
        //     free(xstar);
        //     return 0; // Non postare la soluzione
        // }

        // Prepara gli indici e i valori per postare la soluzione euristica
        int nzcnt = 0;
        for (int i = 0; i < ncols; i++) {
            if (heuristic_sol[i] > 0.5) nzcnt++;
        }
        printf("Number of non-zero variables in heuristic solution: %d\n", nzcnt);

        int *ind = (int *)malloc(inst->ncols * sizeof(int));
        if (!ind) {
            free(heuristic_sol);
            free(succ);
            free(comp);
            free(xstar);
            print_error("Memory allocation error for ind/val");
        }
        int pos = 0;
        for (int i = 0; i < ncols; i++) {
                ind[i] = i;
        }
        // printf("Posting heuristic solution:\n");
        // for (int i = 0; i < nzcnt; i++) {
        //     printf("  ind[%d] = %d\n", i, ind[i]);
        // }
        printf("Objective value of heuristic solution: %f\n", heuristic_objval);
        printf("nuber of cols: %d\n", inst->ncols);
        //print ind and val
        // for (int i = 0; i < nzcnt; i++) {
        //     printf("  ind[%d] = %d\n", i, ind[i]);
        // }
        // Print heuristic solution
        // for (int i = 0; i < inst->ncols; i++) {
        //     printf("heuristic_sol[%d] = %f\n", i, heuristic_sol[i]);
        // }
        int status = CPXcallbackpostheursoln(context, inst->ncols, ind, heuristic_sol, heuristic_objval, CPXCALLBACKSOLUTION_NOCHECK);
        if (status != 0) {
            printf("Warning: CPXcallbackpostheursoln failed with status %d. Check for feasibility issues.\n", status);
        } else {
            printf("Heuristic solution posted successfully!\n");
        }
        free(ind);
        free(heuristic_sol);
        printf("Heuristic solution posted.\n");
    }

    // Aggiungi sempre le SEC constraints
    add_SEC_constraints(inst, NULL, NULL, xstar, context, contextid);
    printf("SEC constraints added.\n");
    free(succ);
    free(comp);
    free(xstar);
    printf("Callback function completed.\n");
    return 0;
}

/*
How to view the solution in Gnuplot:
type

gnuplot

then

set title "TSP Solution"
set xlabel "X"
set ylabel "Y"
plot "../data/solution_gnuplot.dat" with linespoints lt 1 lw 2 pt 7

*/
