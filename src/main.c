#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "tsp_utils.h"
#include "heuristics.h"
#include "chrono.h"
#include "pthread.h"
        

typedef struct {
    bool use_tabu_search;
    bool use_vns_search;
    double min_tenure_dimension;
    double max_tenure_dimension;
    double increase_ten_dim_rate;
    double learning_rate;
    int max_jumps;
} ConfigParams;

void read_input(instance *inst);
void parse_command_line(int argc, char **argv, instance *inst,
	ConfigParams *params
); 

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

    double t1 = second(); 
    instance inst;
    inst.start_time = t1;

	ConfigParams params = {
        .use_tabu_search = false,
        .use_vns_search = false,
        .min_tenure_dimension = 0.1,
        .max_tenure_dimension = 0.5,
        .increase_ten_dim_rate = 1.0,
        .learning_rate = 1.0,
        .max_jumps = 1,
    };

    // General calls for initial calculations and input reading
    parse_command_line(argc, argv, &inst, &params);  
    read_input(&inst);  
    calculate_distances(&inst);
    
    // Nearest neighbor heuristic and gnuplot output
    nearest_neighbor(&inst, false);
    png_solution_for_gnuplot(inst.best_sol, true, "../data/nearest_neighbor", &inst);

    // Tabu Search call with tenure parameters
	if(params.use_tabu_search){
    	tabu_search(inst.best_sol, 20, &inst, params.min_tenure_dimension, params.max_tenure_dimension, params.increase_ten_dim_rate);
	}

	// Variable Neighborhood Search call with parameters
	if(params.use_vns_search){
		variable_neighborhood_search(inst.best_sol, 20, &inst, params.learning_rate, params.max_jumps);
	}

    plot_costs("../data/tabu_costs.txt", "../data/tabu_costs");
    png_solution_for_gnuplot(inst.best_sol, true, "../data/tabu", &inst);
    
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
        int in_node_section = 0;  // Flag to track when we are in the node section
        int count = 0;
    
        // Read the file line by line
        while (fgets(line, sizeof(line), fp)) {
            line[strcspn(line, "\r\n")] = 0; // Removes newline characters from the line
    
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
        // Set the number of nodes (random value between 50 and 100)
        inst->xcoord = malloc(inst->nnodes * sizeof(double));  
        inst->ycoord = malloc(inst->nnodes * sizeof(double));  
        strncpy(inst->input_file, "randomly generated", 1000);  
        int seed = inst->seed;  
        srand(seed);  

        // Generate random coordinates between 0 and 9999
        for (int i = 0; i < inst->nnodes; i++) {
            inst->xcoord[i] = (double)((double)rand() / RAND_MAX) * 10000;  
            inst->ycoord[i] = (double)((double)rand() / RAND_MAX) * 10000; 
        }
    }       
}


void parse_command_line(int argc, char **argv, instance *inst, 
		ConfigParams *params) 
	{
	if (VERBOSE >= 100)
	printf("Running %s with %d parameters\n", argv[0], argc - 1);

	// Default values
	strcpy(inst->input_file, "NULL");
	inst->seed = 0;
	inst->nnodes = -1;
	inst->time_limit = 30;
	inst->best_sol = (solution *)malloc(sizeof(solution));
	if (!inst->best_sol) {
	print_error("Memory allocation failed for best_sol");
	}
	inst->best_sol->tour_cost = 1e+20;
	int got_input_file = 0;

	// Parse arguments dynamically
	for (int i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-file") == 0 && i + 1 < argc) {
	strcpy(inst->input_file, argv[++i]);
	got_input_file = 1;
	} else if (strcmp(argv[i], "-seed") == 0 && i + 1 < argc) {
	inst->seed = abs(atoi(argv[++i]));
	} else if (strcmp(argv[i], "-num_nodes") == 0 && i + 1 < argc) {
	inst->nnodes = atoi(argv[++i]);
	} else if (strcmp(argv[i], "-tl") == 0 && i + 1 < argc) {
	inst->time_limit = atof(argv[++i]);
	}  else if (strcmp(argv[i], "-alg") == 0) {
		const char *alg = argv[++i];
		if (strcmp(alg, "tabu_search") == 0) {
			params->use_tabu_search = true;
			params->min_tenure_dimension = atof(argv[++i]);
			params->max_tenure_dimension = atof(argv[++i]);
			params->increase_ten_dim_rate = atof(argv[++i]);
		}
		else if (strcmp(alg, "vns") == 0) {
			params->use_vns_search = true;
			params->learning_rate = atof(argv[++i]);
			params->max_jumps = atoi(argv[++i]);
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
