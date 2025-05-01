#include "callback.h"
#include "cpx_utils.h"
#include "tsp_utils.h"
#include "heuristics.h"

void save_solution_for_gnuplot(const char *filename, double *heuristic_sol, instance *inst) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Unable to open file %s for writing.\n", filename);
        return;
    }
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = i + 1; j < inst->nnodes; j++) {
            if (heuristic_sol[xpos(i, j, inst)] > 0.5) {
                fprintf(file, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
                fprintf(file, "%f %f\n\n", inst->xcoord[j], inst->ycoord[j]); 
            }
        }
    }
    fclose(file);
    printf("Solution saved for Gnuplot in %s\n", filename);
}


int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    instance *inst = (instance *)userhandle;
    int ncols = inst->ncols;
    double *xstar = (double *)malloc(ncols * sizeof(double));
    if (!xstar) print_error("Memory allocation error");

    double objval = CPX_INFBOUND;
    if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        my_candidate_callback(context, contextid, inst);
    } else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
        if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->ncols - 1, &objval)) {
            free(xstar);
            print_error("CPXcallbackgetrelaxationpoint error");
        }
    } else {
        free(xstar);
        return 0; 
    }

    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    if (!succ || !comp) print_error("Memory allocation error for succ/comp");

    int ncomp = 0;
    build_sol(xstar, inst, succ, comp, &ncomp);

    if (ncomp > 1 && contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
        printf("Subtours detected (%d components), trying to post a heuristic solution...\n", ncomp);

        nearest_neighbor(inst, 0, true);
        double time_limit = 5;
        double learning_rate = 0.01;
        int max_jumps = 5;
        variable_neighborhood_search(inst->best_sol, time_limit, inst, learning_rate, max_jumps);

        double *heuristic_sol = (double *)calloc(ncols, sizeof(double));
        if (!heuristic_sol) print_error("Memory allocation error for heuristic_sol");

        for (int i = 0; i < inst->nnodes - 1; i++) {
            int node_i = (int)(inst->best_sol->tour[i]) - 1;
            int node_j = (int)(inst->best_sol->tour[i + 1]) - 1;
            int index = xpos(node_i, node_j, inst);
            heuristic_sol[index] = 1.0;
        }
        heuristic_sol[xpos((int)(inst->best_sol->tour[inst->nnodes - 1]) - 1, (int)(inst->best_sol->tour[0]) - 1, inst)] = 1.0;

        double heuristic_objval = 0.0;
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (heuristic_sol[xpos(i, j, inst)] > 0.5)
                    heuristic_objval += dist(i, j, inst);
            }
        }

        int *degree = (int *)calloc(inst->nnodes, sizeof(int));
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (heuristic_sol[xpos(i, j, inst)] > 0.5) {
                    degree[i]++;
                    degree[j]++;
                }
            }
        }
        free(degree);

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
        printf("Objective value of heuristic solution: %f\n", heuristic_objval);
        printf("nuber of cols: %d\n", inst->ncols);
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

    add_SEC_constraints(inst, NULL, NULL, xstar, context, contextid);
    free(succ);
    free(comp);
    free(xstar);
    printf("Callback function completed.\n");
    return 0;
}


my_candidate_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, instance *inst) {
    if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) {
        free(xstar);
        print_error("CPXcallbackgetcandidatepoint error");
    }
}