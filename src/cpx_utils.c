#include "cpx_utils.h"
#include "callback.h"  
#include "tsp_utils.h"     
#include "heuristics.h"


void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) {   
	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ ){
		for ( int j = i+1; j < inst->nnodes; j++ ){
			int k = xpos(i,j,inst);
			if ( fabs(xstar[k]) > EPSILON && fabs(xstar[k]-1.0) > EPSILON ) print_error(" wrong xstar in build_sol()");
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

    if (params->warmstart) {
        warmstart(env, lp, inst);
    }

    if (params->use_bc) {
        if (params->fractional) {
            print_error("Fractional not implemented yet");
            // CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
            // if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
            //     print_error("Error while registering the callback");
        }
        else {
            CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
            if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst))
                print_error("Error while registering the callback");
            //branch
        }
    }

    if (params->use_benders) {
        benders_loop(inst, env, lp, start_time);
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


void create_directory(const char *path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        #ifdef _WIN32
        mkdir(path);
        #else
        mkdir(path, 0700);
        #endif
    }
}

void plot_graph_to_image(int nnodes, double *xcoord, double *ycoord, double *xstar, instance *inst, double max_coord, double padding) {
    const char *dir_path = "../data";
    create_directory(dir_path);

    char data_file_path[256];
    snprintf(data_file_path, sizeof(data_file_path), "%s/graph_data.dat", dir_path);
    FILE *file = fopen(data_file_path, "w");
    if (!file) {
        printf("Error: Unable to open file for writing.\n");
        return;
    }

    fprintf(file, "# Selected edges (lines)\n");
    for (int i = 0; i < nnodes; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            if (xstar[xpos(i, j, inst)] > 0.5) { 
                fprintf(file, "%.2f %.2f\n", xcoord[i], ycoord[i]); 
                fprintf(file, "%.2f %.2f\n\n", xcoord[j], ycoord[j]); 
            }
        }
    }

    fprintf(file, "\n\n# Node coordinates (points)\n");
    for (int i = 0; i < nnodes; i++) {
        fprintf(file, "%.2f %.2f\n", xcoord[i], ycoord[i]); 
    }

    fclose(file);

    double xmin = 0 - padding;
    double xmax = max_coord + padding;
    double ymin = 0 - padding;
    double ymax = max_coord + padding;

    FILE *gnuplot = popen("gnuplot", "w");
    if (!gnuplot) {
        printf("Error: Unable to open GNUplot.\n");
        return;
    }

    char image_file_path[256];
    snprintf(image_file_path, sizeof(image_file_path), "%s/graph_plot.png", dir_path);

    fprintf(gnuplot, "set terminal png size 800,600\n"); 
    fprintf(gnuplot, "set output '%s'\n", image_file_path);
    fprintf(gnuplot, "set title 'Graph Visualization'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "set key off\n");
    fprintf(gnuplot, "set size square\n");
    fprintf(gnuplot, "plot \\\n");
    fprintf(gnuplot, "  '%s' index 0 using 1:2 with lines lc rgb 'red' title 'Edges', \\\n", data_file_path);
    fprintf(gnuplot, "  '%s' index 1 using 1:2 with points pt 7 lc rgb 'blue' title 'Nodes'\n", data_file_path);

    pclose(gnuplot);

    printf("Graph saved as an image: %s\n", image_file_path);
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
        double end_time = 0.0;
        CPXgettime(env, &end_time);
        double elapsed_time = end_time - start_time;

        printf("Subtours detected: %d components. Adding SEC constraints...\n", ncomp);

        add_SEC_constraints(inst, env, lp, xstar, NULL, -1);

        free(xstar);
        xstar = (double *)calloc(inst->ncols, sizeof(double));
        if (xstar == NULL) print_error("calloc() error for xstar (loop)");

        if (CPXmipopt(env, lp)) print_error("CPXmipopt() error (after SEC)");
        if (CPXgetx(env, lp, xstar, 0, inst->ncols - 1)) print_error("CPXgetx() error (after SEC)");

        build_sol(xstar, inst, succ, comp, &ncomp);

        if (elapsed_time >= inst->time_limit) {
            printf("Time limit exceeded. Initiating patching process...\n");
            return;
        }
    }
    solution *cpx_sol = (solution *)malloc(sizeof(solution));
    cpx_sol->tour = (double *)calloc(inst->nnodes + 1, sizeof(double));
    if (cpx_sol->tour == NULL) print_error("calloc() error for cpx_sol.tour in benders loop");
    convert_succ_to_tour(succ, inst->nnodes, cpx_sol->tour);

    png_solution_for_gnuplot(cpx_sol, true, "../data/cplex_solution", inst);

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

    printf("Adding %d SEC constraints%s...\n", ncomp, context ? " via callback" : "");

    for (int k = 1; k <= ncomp; k++) {
        int size = 0;
        for (int i = 0; i < inst->nnodes; i++)
            if (comp[i] == k) size++;
        if (size <= 1) continue;

        int max_nz = size * (size - 1) / 2;
        int *index = (int *)malloc(max_nz * sizeof(int));
        double *value = (double *)malloc(max_nz * sizeof(double));
        int nz = 0;

        for (int i = 0; i < inst->nnodes; i++) {
            if (comp[i] != k) continue;
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (comp[j] != k) continue;
                index[nz] = xpos(i, j, inst);
                value[nz++] = 1.0;
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

void invert_path(int start, int end, int *succ, double *xstar, instance *inst) {
    int current = start;
    int prev = -1;

    while (current != end) {
        int next = succ[current]; 
        succ[current] = prev;    
        prev = current;          
        current = next;          
    }

    succ[current] = prev;

    current = start;
    while (current != end) {
        int next = succ[current];
        xstar[xpos(current, next, inst)] = 1.0; 
        xstar[xpos(next, current, inst)] = 0.0; 
        current = next;
    }
}

void patch_solution(double *xstar, instance *inst) {
    int *succ = (int *) malloc(inst->nnodes * sizeof(int));
    int *comp = (int *) malloc(inst->nnodes * sizeof(int));
    int ncomp;

    build_sol(xstar, inst, succ, comp, &ncomp);

    if (ncomp <= 1) {
        printf("No patching needed. Solution is a single tour.\n");
        free(succ);
        free(comp);
        return;
    }

    printf("Patching needed. Components found: %d\n", ncomp);

    while (ncomp > 1) {

        double best_cost = DBL_MAX;
        int best_i = -1, best_j = -1;
        int succ_i = -1, succ_j = -1;
        bool swap = false;

        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                if (comp[i] != comp[j]) {
                    swap = false; 
                    succ_i = succ[best_i];
                    succ_j = succ[best_j];
                    double cij = dist(i, j, inst) + dist(succ_j, succ_i, inst); 
                    if (cij < best_cost) {
                        best_cost = cij;
                        best_i = i;
                        best_j = j;
                    }
                    double cij_swap = dist(i, succ_j, inst) + dist(j, succ_i, inst);
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

        xstar[xpos(best_i, succ_i, inst)] = 0.0;
        xstar[xpos(best_j, succ_j, inst)] = 0.0;

        if  (swap == false){
            xstar[xpos(best_i, best_j, inst)] = 1.0;
            xstar[xpos(succ_j, succ_i, inst)] = 1.0;
        }
        else if (swap == true) { 
            printf("Inverting path from %d to %d\n", best_j, succ_j);
            invert_path(best_j, succ_j, succ, xstar, inst);
        
            xstar[xpos(best_i, succ_j, inst)] = 1.0;
            xstar[xpos(best_j, succ_i, inst)] = 1.0;
        }
        build_sol(xstar, inst, succ, comp, &ncomp);
        printf("Components after patching: %d\n", ncomp);
    }

    free(succ);
    free(comp);
}

void warmstart(CPXENVptr env, CPXLPptr lp, instance *inst) {

        nearest_neighbor(inst, 0, true);
        if (VERBOSE >= 30) printf("Best tour cost for warm start after nearest neighbor: %lf\n", inst->best_sol->tour_cost);
    
        double time_limit;
        if (CPXgetdblparam(env, CPX_PARAM_TILIM, &time_limit)) {
            print_error("CPXgetdblparam() error in warm start");
        }
        double learning_rate = 0.001;
        int max_jumps = 5;
        if (VERBOSE >= 30) printf("Starting Variable Neighborhood Search for warm start...\n");
        variable_neighborhood_search(inst->best_sol, time_limit / 10, inst, learning_rate, max_jumps);
    
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