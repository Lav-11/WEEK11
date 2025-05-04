#include "callback.h"
#include "cpx_utils.h"
#include "tsp_utils.h"
#include "heuristics.h"


int CPXPUBLIC candidate_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle) {
    cpx_data *data = (cpx_data *)userhandle;
    instance *inst = data->inst;
    ConfigParams *params = data->params;
   
    int ncols = inst->ncols;
    double *xstar = (double *)malloc(ncols * sizeof(double));
    if (!xstar) print_error("Memory allocation error");
    double objval = CPX_INFBOUND;    
    if(CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols-1, &objval)) print_error("CPXcallbackgetcandidatepoint error in candidate callback");

    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    if (!succ || !comp) print_error("Memory allocation error for succ/comp");

    int ncomp = 0;
    build_sol(xstar, inst, succ, comp, &ncomp);

    if (ncomp > 1) {
        //printf("Subtours detected (%d components), trying to post a heuristic solution...\n", ncomp);
        add_SEC_constraints(inst, NULL, NULL, xstar, context, contextid);

        // nearest_neighbor(inst, 0, true);
        // double time_limit = inst->time_limit;
        // double learning_rate = 0.01;
        // int max_jumps = 5;
        // variable_neighborhood_search(inst->best_sol, time_limit/10, inst, learning_rate, max_jumps);

        // double *posting_sol = (double *)calloc(ncols, sizeof(double));
        // if (!posting_sol) print_error("Memory allocation error for posting_sol");

        // for (int i = 0; i < inst->nnodes - 1; i++) {
        //     int node_i = (int)(inst->best_sol->tour[i]) - 1;
        //     int node_j = (int)(inst->best_sol->tour[i + 1]) - 1;
        //     int index = xpos(node_i, node_j, inst);
        //     posting_sol[index] = 1.0;
        // }
        // posting_sol[xpos((int)(inst->best_sol->tour[inst->nnodes - 1]) - 1, (int)(inst->best_sol->tour[0]) - 1, inst)] = 1.0;

        // double heuristic_objval = 0.0;
        // for (int i = 0; i < inst->nnodes; i++) {
        //     for (int j = i + 1; j < inst->nnodes; j++) {
        //         if (posting_sol[xpos(i, j, inst)] > 0.5)
        //             heuristic_objval += dist(i, j, inst);
        //     }
        // }

        // int *degree = (int *)calloc(inst->nnodes, sizeof(int));
        // for (int i = 0; i < inst->nnodes; i++) {
        //     for (int j = i + 1; j < inst->nnodes; j++) {
        //         if (posting_sol[xpos(i, j, inst)] > 0.5) {
        //             degree[i]++;
        //             degree[j]++;
        //         }
        //     }
        // }
        // free(degree);

        // int nzcnt = 0;
        // for (int i = 0; i < ncols; i++) {
        //     if (posting_sol[i] > 0.5) nzcnt++;
        // }
        // printf("Number of non-zero variables in heuristic solution: %d\n", nzcnt);

        // int *ind = (int *)malloc(inst->ncols * sizeof(int));
        // if (!ind) {
        //     free(posting_sol);
        //     free(succ);
        //     free(comp);
        //     free(xstar);
        //     print_error("Memory allocation error for ind/val");
        // }
        // int pos = 0;
        // for (int i = 0; i < ncols; i++) {
        //         ind[i] = i;
        // }
        // printf("Objective value of heuristic solution: %f\n", heuristic_objval);
        // printf("nuber of cols: %d\n", inst->ncols);
        // int status = CPXcallbackpostheursoln(context, inst->ncols, ind, posting_sol, heuristic_objval, CPXCALLBACKSOLUTION_NOCHECK);
        // if (status != 0) {
        //     printf("Warning: CPXcallbackpostheursoln failed with status %d. Check for feasibility issues.\n", status);
        // } else {
        //     printf("Heuristic solution posted successfully!\n");
        // }
        // free(ind);
        // free(posting_sol);
        //printf("Heuristic solution posted.\n");
    } else if (ncomp == 1 && params->posting == true) {
        solution *sol = (solution *)malloc(sizeof(solution));
        if (!sol) print_error("Memory allocation error for solution");
        sol->tour = (double *)malloc((inst->nnodes + 1) * sizeof(double));
        if (!sol->tour) print_error("Memory allocation error for solution tour");
        convert_succ_to_tour(succ, inst->nnodes, sol->tour);

        two_opt(sol, 1e20, inst);
        check_tour_feasability(sol, inst);

        double *posting_sol = (double *)calloc(ncols, sizeof(double));
        if (!posting_sol) print_error("Memory allocation error for posting_sol");

        for (int i = 0; i < inst->nnodes - 1; i++) {
            int node_i = (int)(sol->tour[i]) - 1;
            int node_j = (int)(sol->tour[i + 1]) - 1;
            int index = xpos(node_i, node_j, inst);
            posting_sol[index] = 1.0;
        }
        posting_sol[xpos((int)(sol->tour[inst->nnodes - 1]) - 1, (int)(sol->tour[0]) - 1, inst)] = 1.0;

        double heuristic_objval = 0.0;
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (posting_sol[xpos(i, j, inst)] > 0.5)
                    heuristic_objval += inst->distances[i * inst->nnodes + j];
            }
        }

        int nzcnt = 0;
        for (int i = 0; i < ncols; i++) {
            if (posting_sol[i] > 0.5) nzcnt++;
        }
        int *ind = (int *)malloc(inst->ncols * sizeof(int));
        if (!ind) {
            free(posting_sol);
            free(succ);
            free(comp);
            free(xstar);
            print_error("Memory allocation error for ind/val");
        }
        int pos = 0;
        for (int i = 0; i < ncols; i++) {
                ind[i] = i;
        }
        int status = CPXcallbackpostheursoln(context, inst->ncols, ind, posting_sol, heuristic_objval, CPXCALLBACKSOLUTION_NOCHECK);
        if (status != 0) {
            printf("Warning: CPXcallbackpostheursoln failed with status %d. Check for feasibility issues.\n", status);
        } else {
            printf("Heuristic solution posted successfully with cost: %f!\n", heuristic_objval);
        }
        free(ind);
        free(posting_sol);
    }

    free(succ);
    free(comp);
    free(xstar);
    return 0;
}


// my_candidate_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, instance *inst) {
//     if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) {
//         free(xstar);
//         print_error("CPXcallbackgetcandidatepoint error");
//     }
// }