#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include "tsp_utils.h"
#include "heuristics.h"
#include "chrono.h"
#include "pthread.h"
#include "cpx_utils.h"
#include "callback.h"
#include "tsp_utils.h"
#include "solver.h"
#include <ilcplex/cplex.h>

BestResult global_best_result = { .best_cost = INFINITY, .best_params = "", .best_cost_mutex = PTHREAD_MUTEX_INITIALIZER };

void read_input(instance *inst);
void parse_command_line(int argc, char **argv, instance *inst,
	ConfigParams *params
); 
void *thread_function(void *arg);
void heuristics_multiparams(instance *inst, ConfigParams *params, int argc, char **argv, double t1);
void cpx_multiparams(instance *inst, ConfigParams *params, int argc, char **argv);

int main(int argc, char **argv) 
{ 	
    if (argc < 2) { 
        printf("Usage: %s -help for help\n", argv[0]); 
        exit(1); 
    }       
    if (VERBOSE >= 2) { 
        for (int a = 0; a < argc; a++) printf("%s ", argv[a]); 
        printf("\n"); 
    }

	ConfigParams params;
	set_default_params(&params);

	instance inst;
	printf("Starting TSP solver...\n");
	parse_command_line(argc, argv, &inst, &params);
	read_input(&inst);
	calculate_distances(&inst);

	if (inst.nnodes < 0) {
		printf("ERROR: Number of random nodes or input file not specified.\n");
		exit(1);
	}
    
	double t1 = second(); 

	if (inst.use_cplex) {
		cpx_multiparams(&inst, &params, argc, argv);
	}
	else if (inst.use_heuristics){
		heuristics_multiparams(&inst, &params, argc, argv, t1);
	}


	double t2 = second(); 

	if (VERBOSE >= 1)   
	{
		printf("TSP problem calculations terminated in %lf sec.s\n", t2 - t1);  
	}

	return 0;
}         


void read_input(instance *inst)  
{

    if (inst->nnodes < 0) {
        FILE *fp = fopen(inst->input_file, "r");  
        if (!fp) {
            print_error("Error opening file");
        }
    
        inst->nnodes = 0;
        inst->xcoord = NULL;
        inst->ycoord = NULL;
    
        char line[1024];
        int in_node_section = 0;  
        int count = 0;
    
        // Read the file line by line
        while (fgets(line, sizeof(line), fp)) {
            line[strcspn(line, "\r\n")] = 0; 
    
            // Extract the number of nodes from the "DIMENSION" line
            if (strncmp(line, "DIMENSION", 9) == 0) {
                char *ptr = strchr(line, ':');
                if (ptr)
                    inst->nnodes = atoi(ptr + 1);  
                else {
                    char dummy[100];
                    sscanf(line, "%s %d", dummy, &inst->nnodes);  
                }
                inst->xcoord = malloc(inst->nnodes * sizeof(double));
                inst->ycoord = malloc(inst->nnodes * sizeof(double));
                if (!inst->xcoord || !inst->ycoord)
                    print_error("Memory allocation error for coordinate arrays");
            }
            // When we reach the "NODE_COORD_SECTION", start reading coordinates
            else if (strcmp(line, "NODE_COORD_SECTION") == 0) {
                in_node_section = 1;
                continue;
            } else if (in_node_section) {
                if (strcmp(line, "EOF") == 0)
                    break;  
                int id;
                double x, y;
                // Read the node ID and coordinates
                if (sscanf(line, "%d %lf %lf", &id, &x, &y) < 3)
                    continue;  // Skip invalid lines
                if (count < inst->nnodes) {
                    inst->xcoord[count] = x;
                    inst->ycoord[count] = y;
                    count++;
                }
            }
        }
    
        fclose(fp);  
    } else {
        inst->integer_costs = false;
		inst->xcoord = malloc(inst->nnodes * sizeof(double));  
        inst->ycoord = malloc(inst->nnodes * sizeof(double));  
        strncpy(inst->input_file, "rand", 1000);  
        srand(inst->seed);
        if (VERBOSE >= 50)  
            printf("Generating random instance with %d nodes\n", inst->nnodes);
        for (int i = 0; i < inst->nnodes; i++) {
            inst->xcoord[i] = (double)((double)rand() / RAND_MAX) * inst->max_coord;  
            inst->ycoord[i] = (double)((double)rand() / RAND_MAX) * inst->max_coord; 
        }
    }       
}


void parse_command_line(int argc, char **argv, instance *inst, 
		ConfigParams *params) 
	{
	if (VERBOSE >= 60)
	printf("Running %s with %d parameters\n", argv[0], argc - 1);

	// Default values
	strcpy(inst->input_file, "NULL");
	inst->use_cplex = false;
	inst->use_heuristics = false;
	inst->seed = 1;
	inst->nnodes = -1;
	inst->time_limit = 10;
	inst->xcoord = NULL;
	inst->ycoord = NULL;
	inst->distances = NULL;
	inst->best_sol = NULL;
	inst->integer_costs = true;
	inst->max_coord = 10000.0; 
	int got_input_file = 0;

	// Parse arguments 
	for (int i = 1; i < argc; i++) {
	if ((strcmp(argv[i], "-file") == 0 || strcmp(argv[i], "-f") == 0) && i + 1 < argc) {
	strcpy(inst->input_file, argv[++i]);
	got_input_file = 1;
	} else if (strcmp(argv[i], "-use_cplex") == 0) {
	inst->use_cplex = true;
	} else if (strcmp(argv[i], "-warmstart") == 0) {
    params->warmstart = true;
    } else if (strcmp(argv[i], "-use_benders") == 0) {
    params->use_benders = true;
    } else if (strcmp(argv[i], "-use_bc") == 0) {
    params->use_bc = true;
    } else if (strcmp(argv[i], "-posting") == 0) {
    params->posting = true;
    } else if (strcmp(argv[i], "-fractional") == 0) {
    params->fractional = true;
    } else if (strcmp(argv[i], "-use_heu") == 0) {
	inst->use_heuristics = true;
	} else if (strcmp(argv[i], "-seed") == 0 && i + 1 < argc) {
	inst->seed = abs(atoi(argv[++i]));
	} else if (strcmp(argv[i], "-sequential_seed") == 0 && i + 1 < argc) {
	params->sequential_seed = abs(atoi(argv[++i]));
	} else if (strcmp(argv[i], "-num_nodes") == 0 && i + 1 < argc) {
	inst->nnodes = atoi(argv[++i]);
	} else if (strcmp(argv[i], "-tl") == 0 && i + 1 < argc) {
	inst->time_limit = atof(argv[++i]);
	}  else if (strcmp(argv[i], "-alg") == 0) {
		const char *alg = argv[++i];
		if (strcmp(alg, "tabu") == 0) {
			params->use_tabu_search = true;
			params->min_tenure_dimension_lower_bound = atof(argv[++i]);
			params->min_tenure_dimension_higher_bound = atof(argv[++i]);
			params->min_tenure_dimension_delta = atof(argv[++i]);
			params->max_tenure_dimension_lower_bound = atof(argv[++i]);
			params->max_tenure_dimension_higher_bound = atof(argv[++i]);
			params->max_tenure_dimension_delta = atof(argv[++i]);
		}
		else if (strcmp(alg, "vns") == 0) {
			params->use_vns_search = true;
			params->learning_rate_lower_bound = atof(argv[++i]);
			params->learning_rate_higher_bound = atof(argv[++i]);
			params->learning_rate_delta = atof(argv[++i]);
			params->max_jumps_lower_bound = atoi(argv[++i]);
			params->max_jumps_higher_bound = atoi(argv[++i]);
			params->max_jumps_delta = atoi(argv[++i]);
		}


	} else if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
	printf("Available options:\n");
	printf("-file <input_file>: Specify input file.\n");
	printf("-seed <value>: Random seed.\n");
	printf("-num_nodes <value>: Number of nodes.\n");
	printf("-tl <value>: Time limit.\n");
	printf("-alg <tabu/vns>: Specify algorithm to run.\n");
	printf("-min_tenure_dimension <value>: Minimum tenure dimension for Tabu Search.\n");
	printf("-max_tenure_dimension <value>: Maximum tenure dimension for Tabu Search.\n");
	printf("-increase_ten_dim_rate <value>: Rate for Tabu tenure increase.\n");
	printf("-vns_param1 <value>: First parameter for VNS.\n");
	printf("-vns_param2 <value>: Second parameter for VNS.\n");
	exit(0);
	}
	}

	// Validation
	if (inst->nnodes < 0 && !got_input_file) {
	printf("ERROR: Number of random nodes or input file not specified.\n");
	exit(1);
	}
}


// Function to start calculations of heuristics
void heuristics_multiparams(instance *inst, ConfigParams *params, int argc, char **argv, double t1) {
    int max_threads = sysconf(_SC_NPROCESSORS_ONLN) -2;
    if (max_threads < 1) max_threads = 1;

    if (VERBOSE >= 10) {
        printf("Multithreading in use\n");
        printf("Max threads: %d\n", max_threads + 2);
        printf("Threads used: %d\n", max_threads);
    }

    pthread_t threads[max_threads];
    ThreadData thread_data[max_threads];
    int thread_count = 0;
    int combination_count = 0;

    for (int i = 0; i <= params->sequential_seed; i++) {
        instance inst_copy;
        parse_command_line(argc, argv, &inst_copy, params);
        inst_copy.seed = inst->seed + i;
        read_input(&inst_copy);
        calculate_distances(&inst_copy);
        nearest_neighbor(&inst_copy, 0, false);

        if (VERBOSE >= 30) {
            if (params->sequential_seed > 0) {
                printf("Using seed: %d\n", inst_copy.seed);
            }
            printf("Best tour cost after nearest neighbor: %lf\n", inst_copy.best_sol->tour_cost);
        }

        if (params->use_tabu_search || params->use_vns_search) {
            if (params->use_tabu_search) {
                for (double min_tenure = params->min_tenure_dimension_lower_bound; min_tenure <= params->min_tenure_dimension_higher_bound + EPSILON; min_tenure += params->min_tenure_dimension_delta) {
                    for (double max_tenure = params->max_tenure_dimension_lower_bound; max_tenure <= params->max_tenure_dimension_higher_bound + EPSILON; max_tenure += params->max_tenure_dimension_delta) {
                        instance *inst_thread = copy_instance(&inst_copy);
                        thread_data[thread_count] = (ThreadData){
                            .inst = inst_thread,
                            .params = *params,
                            .sol = inst_thread->best_sol,
                            .min_tenure = min_tenure,
                            .max_tenure = max_tenure,
                            .is_tabu_search = true,
                            .is_vns_search = false
                        };
                        pthread_create(&threads[thread_count], NULL, thread_function, &thread_data[thread_count]);
                        thread_count++;
						combination_count++;

                        if (thread_count >= max_threads) {
                            for (int j = 0; j < thread_count; j++) {
                                pthread_join(threads[j], NULL);
                            }
                            thread_count = 0;
                        }
                    }
                }
            }
            if (params->use_vns_search) {
                for (double learning_rate = params->learning_rate_lower_bound; learning_rate <= params->learning_rate_higher_bound + EPSILON; learning_rate += params->learning_rate_delta) {
                    for (int max_jumps = params->max_jumps_lower_bound; max_jumps <= params->max_jumps_higher_bound + EPSILON; max_jumps += params->max_jumps_delta) {
                        instance *inst_thread = copy_instance(&inst_copy);
                        thread_data[thread_count] = (ThreadData){
                            .inst = inst_thread,
                            .params = *params,
                            .sol = inst_thread->best_sol,
                            .learning_rate = learning_rate,
                            .max_jumps = max_jumps,
                            .is_tabu_search = false,
                            .is_vns_search = true
                        };
                        pthread_create(&threads[thread_count], NULL, thread_function, &thread_data[thread_count]);
                        thread_count++;
						combination_count++;

                        if (thread_count >= max_threads) {
                            for (int j = 0; j < thread_count; j++) {
                                pthread_join(threads[j], NULL);
                            }
                            thread_count = 0;
                        }
                    }
                }
            }
        }
        free_instance(&inst_copy, false);
    }

    for (int i = 0; i < thread_count; i++) {
        pthread_join(threads[i], NULL);
    }
    double t2 = second();
    plot_costs("../data/tabu_costs.txt", "../data/tabu_costs");

    if (VERBOSE >= 1) {
        printf("TSP problem calculations terminated in %lf sec.s\n", t2 - t1);
    }

	printf("Best cost found: %.2f\n", global_best_result.best_cost);

    if (combination_count > 1 + EPSILON) {
        printf("Combinations tested: %d\n", combination_count);
        printf("Parameters of best solution: %s\n", global_best_result.best_params);
    }
}

// Function to run Cplex on multiple instances
void cpx_multiparams(instance *inst, ConfigParams *params, int argc, char **argv) {
    for (int i = 0; i <= params->sequential_seed; i++) {
        instance inst_copy;
        parse_command_line(argc, argv, &inst_copy, params);
        inst_copy.seed = inst->seed + i;
        read_input(&inst_copy);
        calculate_distances(&inst_copy);
        double t_start = second();
        TSPopt(&inst_copy, params);
        double t_end = second();
        double time_needed = t_end - t_start;
        if (VERBOSE >= 1) {
            printf("Cplex calculations terminated in %lf sec.s\n", t_end - t_start);
        }
        FILE *csv_file = fopen("../data/cpx_times.csv", "a");
        fprintf(csv_file, "cpx_%i_%i_%i,num_nodes_%d_seed_%d_time_limit_%.2f,%.2f\n",
            params->use_benders ? 0 : params->use_bc ? 1 : 2, params->warmstart ? 1 : 0, params->posting ? 1 : 0,
            inst->nnodes, inst_copy.seed, inst->time_limit, time_needed);
		fclose(csv_file);
        free_instance(&inst_copy, false);
    }
}
