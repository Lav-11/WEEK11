#include "cpx_utils.h"
#include "callback.h"  
#include "tsp_utils.h"     
#include "heuristics.h"


void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) {   
	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ ){
		for ( int j = i+1; j < inst->nnodes; j++ ){
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPSILON && fabs(xstar[k]-1.0) > EPSILON ){ 
                printf("xstar[%d] = %f\n", k, xstar[k]);
                print_error(" wrong xstar in build_sol()");
            }
			if ( xstar[k] > 0.5 ) {
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ ){
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}	
	free(degree);

	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ ){
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < inst->nnodes; start++ ){
		if ( comp[start] >= 0 ) continue; 
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done ) {
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ ) {
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ){
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start; 
	}
}

int TSPopt(instance *inst, ConfigParams *params) {

    if (params->use_benders && params->use_bc) {
        print_error("Multiple algorithms have been selected. Please choose only one.");
    }
    int error;
    CPXENVptr env = CPXopenCPLEX(&error);
    if (error) print_error("CPXopenCPLEX() error");
    CPXLPptr lp = CPXcreateprob(env, &error, "TSP model");
    if (error) print_error("CPXcreateprob() error");

    build_model(inst, env, lp);
    inst->ncols = CPXgetnumcols(env, lp);

    if (VERBOSE >= 50) CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
    if (VERBOSE >= 30) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); 
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed);
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);

    double start_time = 0.0;
    CPXgettime(env, &start_time);
    
    double time_limit = 0.0;
    CPXgetdblparam(env, CPX_PARAM_TILIM, &time_limit);
    printf("Time limit: %f seconds\n", time_limit);

    if (params->warmstart) {
        warmstart(env, lp, inst, start_time);
    }

    if (params->use_bc) {
        branch_and_cut(inst, env, lp, start_time, params);
    }

    if (params->use_benders) {
        benders_loop(inst, env, lp, start_time);
    }

    if (params->use_hardfixing) {
        hard_fixing(inst, env, lp, start_time, params);
    }

    // // Print the solution details if verbose
    // if (VERBOSE >= 50) {
    //     printf("Selected edges:\n");
    //     for (int i = 0; i < inst->nnodes; i++) {
    //         for (int j = i + 1; j < inst->nnodes; j++) {
    //             if (xstar[xpos(i, j, inst)] > 0.5)
    //                 printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
    //         }
    //     }

    //     printf("\nArray of successors:\n");
    //     for (int i = 0; i < inst->nnodes; i++) {
    //         printf("    Node %d -> Node %d\n", i + 1, succ[i] + 1);
    //     }
    // }

    // // Rebuild the solution and print the number of components
    // build_sol(xstar, inst, succ, comp, &ncomp);
    // printf("Number of components:  %d \n", ncomp);

    // // Get and print the final objective value
    // double objval;
    // if (CPXgetobjval(env, lp, &objval)) {
    //     print_error("CPXgetobjval() error");
    // }
    // printf("Final objective value: %f\n", objval);

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

int xpos(int i, int j, instance *inst){ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}
	

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp){    

	double zero = 0.0;  
	char binary = 'B'; 

	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));


	for ( int i = 0; i < inst->nnodes; i++ ){       
		for ( int j = i+1; j < inst->nnodes; j++ ){   
			sprintf(cname[0], "x(%d,%d)", i+1,j+1); 
			double obj = dist(i,j,inst);   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s");
		}
	} 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ ){
		double rhs = 2.0;
		char sense = 'E';                            
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ ){
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	if ( VERBOSE >= 100 ) CPXwriteprob(env, lp, "model.lp", NULL);   

}


void benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp, double start_time) {

    double *xstar = (double *)calloc(inst->ncols, sizeof(double));
    if (xstar == NULL) print_error("calloc() error for xstar in benders loop");

    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error in benders loop");
    if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1)) print_error("CPXgetx() error benders loop");

    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    if (succ == NULL) print_error("calloc() error for succ in benders loop");
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    if (comp == NULL) print_error("calloc() error for comp in benders loop");

    int ncomp = -1;
    build_sol(xstar, inst, succ, comp, &ncomp);

    while (ncomp > 1) {
        double time_now = 0.0;
        CPXgettime(env, &time_now);
        double elapsed_time = time_now - start_time;
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - elapsed_time);

        printf("Subtours detected: %d components. Adding SEC constraints...\n", ncomp);

        add_SEC_constraints(inst, env, lp, xstar, NULL, -1);

        free(xstar);
        xstar = (double *)calloc(inst->ncols, sizeof(double));
        if (xstar == NULL) print_error("calloc() error for xstar (loop)");

        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error (after SEC)");
        if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1)) print_error("CPXgetx() error (after SEC)");

        build_sol(xstar, inst, succ, comp, &ncomp);

        if (elapsed_time >= inst->time_limit) {
            printf("Time limit exceeded. Proceding with patching...\n");
            patch_solution(xstar, inst);
            break;
        }
    }
    build_sol(xstar, inst, succ, comp, &ncomp);
    solution *cpx_sol = (solution *)malloc(sizeof(solution));
    cpx_sol->tour = (double *)calloc(inst->nnodes + 1, sizeof(double));
    if (cpx_sol->tour == NULL) print_error("calloc() error for cpx_sol->tour in benders loop");
    convert_succ_to_tour(succ, inst->nnodes, cpx_sol->tour);
    two_opt(cpx_sol, 1e20, inst);

    check_tour_feasability(cpx_sol, inst);
    png_solution_for_gnuplot(cpx_sol, true, "../data/benders_solution", inst);

    free(xstar);
    free(succ);
    free(comp);
    free_solution(cpx_sol);
}

void add_SEC_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, double *xstar,
                         CPXCALLBACKCONTEXTptr context, int contextid) {
    int *succ = (int *)calloc(inst->nnodes, sizeof(int));
    int *comp = (int *)calloc(inst->nnodes, sizeof(int));
    int ncomp;

    build_sol(xstar, inst, succ, comp, &ncomp);

    if (ncomp <= 1) {
        free(succ);
        free(comp);
        return;
    }
    if (VERBOSE >= 50) {
        printf("Adding %d SEC constraints%s...\n", ncomp, context ? " via callback" : "");
    }

    for (int k = 1; k <= ncomp; k++) {
        int size = 0;
        for (int i = 0; i < inst->nnodes; i++)
            if (comp[i] == k) size++;
        if (size <= 1) continue;

        int max_nz = inst->ncols;
        int *index = (int *)malloc(max_nz * sizeof(int));
        double *value = (double *)malloc(max_nz * sizeof(double));
        int nz = 0;

        for (int i = 0; i < inst->nnodes; i++) {
            if (comp[i] != k) continue;
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (comp[j] != k) continue;
                index[nz] = xpos(i, j, inst);
                value[nz] = 1.0;
                nz++;
            }
        }

        double rhs = size - 1.0;
        char sense = 'L';
        int izero = 0;

        if (context) {
            if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
                if (CPXcallbackrejectcandidate(context, 1, nz, &rhs, &sense, &izero, index, value)) {
                    print_error("CPXcallbackrejectcandidate() error");
                }
            } else if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) {
                int purgeable = CPX_USECUT_FILTER;
                int local = 0;
                if (CPXcallbackaddusercuts(context, 1, nz, &rhs, &sense, &izero, index, value, &purgeable, &local)) {
                    print_error("CPXcallbackaddusercuts() error");
                }
            }
        } else {
            char *rowname = (char *)malloc(100 * sizeof(char));
            sprintf(rowname, "SEC_comp_%d", k);
            if (CPXaddrows(env, lp, 0, 1, nz, &rhs, &sense, &izero, index, value, NULL, &rowname)) {
                print_error("CPXaddrows() failed for SEC");
            }
            free(rowname);
        }
        free(index);
        free(value);
    }
    free(succ);
    free(comp);
}


void patch_solution(double *xstar, instance *inst) {
    int *succ = (int *) malloc(inst->nnodes * sizeof(int));
    int *comp = (int *) malloc(inst->nnodes * sizeof(int));
    int ncomp;

    build_sol(xstar, inst, succ, comp, &ncomp);

    if (ncomp <= 1) {
        if (VERBOSE >= 50) {
            printf("No patching needed. Solution is a single tour.\n");
        }
        free(succ);
        free(comp);
        return;
    }

    if (VERBOSE >= 30) {
        printf("Patching needed. Components found: %d\n", ncomp);
    }

    while (ncomp > 1) {
        double best_cost = DBL_MAX;
        int best_i = -1, best_j = -1;
        int succ_i = -1, succ_j = -1;
        bool swap = false;

        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (comp[i] != comp[j]) {
                    succ_i = succ[i];
                    succ_j = succ[j];
                    double cij = inst->distances[i * inst->nnodes + j] + inst->distances[succ_i * inst->nnodes + succ_j] - inst->distances[i * inst->nnodes + succ_i] - inst->distances[j * inst->nnodes + succ_j]; 
                    if (cij < best_cost) {
                        swap = false;
                        best_cost = cij;
                        best_i = i;
                        best_j = j;
                    }
                    double cij_swap = inst->distances[i * inst->nnodes + succ_j] + inst->distances[j * inst->nnodes + succ_i] - inst->distances[i * inst->nnodes + succ_i] - inst->distances[j * inst->nnodes + succ_j];
                    if (cij_swap < best_cost) {
                        best_cost = cij_swap;
                        best_i = i;
                        best_j = j;
                        swap = true; 
                    }
                }
            }
        }

        if (best_i == -1 || best_j == -1) {
            printf("No valid edge pairs found. Aborting patching.\n");
            break;
        }
        
        succ_i = succ[best_i];
        succ_j = succ[best_j];
        xstar[xpos(best_i, succ_i, inst)] = 0.0;
        xstar[xpos(best_j, succ_j, inst)] = 0.0;
        if  (swap == false){
            xstar[xpos(best_i, best_j, inst)] = 1.0;
            xstar[xpos(succ_i, succ_j, inst)] = 1.0;
        }
        else if (swap == true) { 
            xstar[xpos(best_i, succ_j, inst)] = 1.0;
            xstar[xpos(best_j, succ_i, inst)] = 1.0;
        }
        build_sol(xstar, inst, succ, comp, &ncomp);
    }
    
    free(succ);
    free(comp);
}

void warmstart(CPXENVptr env, CPXLPptr lp, instance *inst, double start_time) {

        nearest_neighbor(inst, 0, true);
        if (VERBOSE >= 30) printf("Best tour cost for warm start after nearest neighbor: %lf\n", inst->best_sol->tour_cost);
    
        double time_limit;
        if (CPXgetdblparam(env, CPX_PARAM_TILIM, &time_limit)) {
            print_error("CPXgetdblparam() error in warm start");
        }
        double learning_rate = 0.001;
        int max_jumps = 5;
        if (VERBOSE >= 30) printf("Starting Variable Neighborhood Search for warm start...\n");
        variable_neighborhood_search(inst->best_sol, (inst->time_limit / 10), inst, learning_rate, max_jumps);
    
        double *mip_start = (double *)calloc(inst->ncols, sizeof(double));
        if (mip_start == NULL) print_error("Memory allocation failed in warm start");
    
        for (int i = 0; i < inst->nnodes - 1; i++) {
            int node_i = (int)(inst->best_sol->tour[i]) - 1; 
            int node_j = (int)(inst->best_sol->tour[i + 1]) - 1;  
            mip_start[xpos(node_i, node_j, inst)] = 1.0;  
        }
        mip_start[xpos((int)(inst->best_sol->tour[inst->nnodes - 1]) - 1, (int)(inst->best_sol->tour[0]) - 1, inst)] = 1.0;
    
        int effort_level = CPX_MIPSTART_AUTO;
        int beg = 0; 	
        int *indices = (int *)malloc(inst->ncols * sizeof(int));
        if (indices == NULL) print_error("Memory allocation failed for indices in warm start");
    
        double *values = (double *)malloc(inst->ncols * sizeof(double));
        if (values == NULL) print_error("Memory allocation failed for values in warm start");
    
        for (int i = 0; i < inst->ncols; i++) {
            indices[i] = i;
            values[i] = mip_start[i];
        }

        int start_status = CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, indices, values, &effort_level, NULL);
        free(indices);
        free(values);
        free(mip_start);
        if (start_status) print_error("CPXaddmipstarts() error");

        double time_now = 0.0;
        CPXgettime(env, &time_now);
        double elapsed_time = time_now - start_time;
        CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - elapsed_time);
}

void convert_succ_to_tour(int *succ, int nnodes, double *tour) {
    if (tour == NULL) {
        print_error("Null pointer passed for tour in convert_succ_to_tour");
    }
    double *temp_tour = (double *)malloc((nnodes + 1) * sizeof(double));
    if (temp_tour == NULL) {
        print_error("Memory allocation failed for temp_tour in convert_succ_to_tour");
    }

    int current = 0; 
    for (int i = 0; i < nnodes; i++) {
        temp_tour[i] = current+1;
        current = succ[current]; 
    }
    temp_tour[nnodes] = temp_tour[0];  
    memcpy(tour, temp_tour, (nnodes + 1) * sizeof(double));
    free(temp_tour);
}

void branch_and_cut(instance *inst, CPXENVptr env, CPXLPptr lp, double start_time, ConfigParams *params) {

    cpx_data *data = (cpx_data *)malloc(sizeof(cpx_data));
    if (data == NULL) print_error("Memory allocation failed for cpx_data in branch_and_cut");
    data->inst = inst;
    data->params = params;
    if (params->fractional) {
        print_error("Fractional not implemented yet");
        // CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
        // if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
        //     print_error("Error while registering the callback");
    }
    else {
        CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
        if (CPXcallbacksetfunc(env, lp, contextid, candidate_callback, data))
            print_error("Error while registering the candidate callback");
    }

    if (CPXmipopt(env, lp)) print_error("CPXmipopt() error in branch_and_cut");

    free(data);
}

void hard_fixing(instance *inst, CPXENVptr env, CPXLPptr lp, double start_time, ConfigParams *params) {
    // Initialization with Nearest Neighbor
    nearest_neighbor(inst, 0, true);
    if (VERBOSE >= 30) {
        printf("Best tour cost after Nearest Neighbor: %lf\n", inst->best_sol->tour_cost);
    }

    // Allocate time for Variable Neighborhood Search (VNS)
    double time_now = 0.0;
    CPXgettime(env, &time_now);
    double time_remaining = inst->time_limit - (time_now - start_time);

    double vns_time_limit = time_remaining * 0.1; // Use 10% of the total time for VNS
    double hard_fixing_time_limit = time_remaining - vns_time_limit; // Remaining time for Hard Fixing

    // Execute Variable Neighborhood Search
    if (VERBOSE >= 30) printf("Starting Variable Neighborhood Search for warm start...\n");
    double learning_rate = 0.001;
    int max_jumps = 5;
    variable_neighborhood_search(inst->best_sol, vns_time_limit, inst, learning_rate, max_jumps);

    // Print the cost of the solution after VNS
    if (VERBOSE >= 30) {
        printf("Best tour cost after Variable Neighborhood Search: %lf\n", inst->best_sol->tour_cost);
    }

    // Prepare for Hard Fixing
    double *xstar = (double *)calloc(inst->ncols, sizeof(double));
    if (xstar == NULL) print_error("Memory allocation failed for xstar in hard fixing");

    for (int i = 0; i < inst->nnodes - 1; i++) {
        int node_i = (int)(inst->best_sol->tour[i]) - 1;
        int node_j = (int)(inst->best_sol->tour[i + 1]) - 1;
        xstar[xpos(node_i, node_j, inst)] = 1.0;
    }
    xstar[xpos((int)(inst->best_sol->tour[inst->nnodes - 1]) - 1, (int)(inst->best_sol->tour[0]) - 1, inst)] = 1.0;

    CPXgettime(env, &time_now);
    time_remaining = hard_fixing_time_limit - (time_now - start_time);

    double max_time_for_run = time_remaining / params->hf_num_of_run;
    double current_prob = params->hf_prob_lower_bound;

    int *modified_indices = (int *)malloc(inst->ncols * sizeof(int));
    if (modified_indices == NULL) print_error("Memory allocation failed for modified_indices");
    int modified_count = 0;

    // Hard Fixing Loop
    while (time_remaining > 0) {
        double iteration_time_limit = fmin(max_time_for_run, time_remaining / params->hf_num_of_run);
        CPXsetdblparam(env, CPX_PARAM_TILIM, iteration_time_limit);

        for (int i = 0; i < inst->ncols; i++) {
            double random_value = rand() / (double)RAND_MAX;
            if (xstar[i] > 0.5 && random_value < current_prob) {
                double fixed_value = 1.0;
                CPXchgbds(env, lp, 1, &i, "L", &fixed_value);
                CPXchgbds(env, lp, 1, &i, "U", &fixed_value);
                modified_indices[modified_count++] = i;
            }
        }

        if (CPXmipopt(env, lp)) {
            print_error("CPXmipopt() error during hard fixing");
        }

        double objval;
        int solstat = CPXgetstat(env, lp);
        if (solstat == CPXMIP_OPTIMAL || solstat == CPXMIP_FEASIBLE || solstat == CPXMIP_TIME_LIM_FEAS) {
            if (CPXgetobjval(env, lp, &objval)) {
                print_error("CPXgetobjval() error during hard fixing");
            }
            printf("Incumbent value after iteration: %f\n", objval);

            if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1)) {
                print_error("CPXgetx() error: Unable to retrieve solution");
            }
        } else {
            printf("No feasible solution found in this iteration. Status: %d\n", solstat);
        }

        for (int j = 0; j < modified_count; j++) {
            int var_index = modified_indices[j];
            double lb = 0.0;
            double ub = 1.0;
            CPXchgbds(env, lp, 1, &var_index, "L", &lb);
            CPXchgbds(env, lp, 1, &var_index, "U", &ub);
        }
        modified_count = 0;

        CPXgettime(env, &time_now);
        time_remaining = hard_fixing_time_limit - (time_now - start_time);

        // Update the probability dynamically
        current_prob = fmin(params->hf_prob_higher_bound, current_prob + params->hf_prob_delta);
    }

    free(modified_indices);
    free(xstar);
}