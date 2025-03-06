#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "tsp_functions.h"  
#include "chrono.h"         

double second();
double random01();     
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst); 

int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	double t1 = second(); 
	instance inst;

	parse_command_line(argc,argv, &inst);    	  
	read_input(&inst);  
	simple_solution(&inst);
	export_solution_for_gnuplot("../data/solution.dat", &inst);
	png_solution_for_gnuplot("../data/solution.dat", "../data/solution.png");

    double t2 = second(); 

	if ( VERBOSE >= 1 )   
	{
		printf("TSP problem solved in %lf sec.s\n", t2-t1);  
	}
	
	return 0; 
}         


void read_input(instance *inst)  
{
	if (inst->nnodes < 0){
		FILE *fp = fopen(inst->input_file, "r");  // Open the file for reading
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
					inst->nnodes = atoi(ptr + 1);  // Extract the number of nodes
				else {
					char dummy[100];
					sscanf(line, "%s %d", dummy, &inst->nnodes);  // Alternative way to extract the number of nodes
				}
				// Allocate memory for the coordinates of the nodes
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
					break;  // End of the node section
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
	
		fclose(fp);  // Close the file
	}
	else {
		// Set the number of nodes (random value between 50 and 100)
		inst->xcoord = malloc(inst->nnodes * sizeof(double));  // Allocate memory for x-coordinates
		inst->ycoord = malloc(inst->nnodes * sizeof(double));  // Allocate memory for y-coordinates
		strncpy(inst->input_file, "randomly generated", 1000);  // Mark as a random instance
		int seed = inst->seed;  // Random seed
		srand(seed);  // Initialize the random number generator
		// Generate random coordinates between 0 and 9999
		for (int i = 0; i < inst->nnodes; i++) {
			inst->xcoord[i] = (double)((double)rand() / RAND_MAX)*10000;  // Random x-coordinate
			inst->ycoord[i] = (double)((double)rand() / RAND_MAX)*10000;  // Random y-coordinate
		}
	
	}		
}


void parse_command_line(int argc, char** argv, instance *inst) 
{ 
	
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	strcpy(inst->input_file, "NULL");
	inst->seed = 0; 
	inst->nnodes = -1;
	inst->timelimit = 1000000; 
	//inst->integer_costs = 0;
	int got_input_file = 0;

    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); got_input_file=1; continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); got_input_file=1; continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); got_input_file=1; continue; } 				// input file
		if ( strcmp(argv[i],"-tl") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		// if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		// if ( strcmp(argv[i],"-old_benders") == 0 ) { inst->old_benders = atoi(argv[++i]); continue; } 	// old benders
		// if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->seed = abs(atoi(argv[++i])); continue; } 				// random seed
		if ( strcmp(argv[i],"-num_nodes") == 0 ) { inst->nnodes = atoi(argv[++i]); continue; } 				// random seed
		// if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		// if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		// if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		// if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		// if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		// if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\navailable parameters (vers. 16-may-2015) --------------------------------------------------\n");
		printf("-file %s\n", inst->input_file); 
		printf("-time_limit %lf\n", inst->timelimit);
		printf("-seed %d\n", inst->seed); 
		printf("-number of nodes %d\n", inst->nnodes); 
		// printf("-model_type %d\n", inst->model_type); 
		// printf("-old_benders %d\n", inst->old_benders); 
		// printf("-seed %d\n", inst->randomseed); 
		// printf("-threads %d\n", inst->num_threads);  
		// printf("-max_nodes %d\n", inst->max_nodes); 
		// printf("-memory %d\n", inst->available_memory); 
		// printf("-int %d\n", inst->integer_costs); 
		// printf("-node_file %s\n", inst->node_file);
		// printf("-cutoff %lf\n", inst->cutoff); 
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if ( inst->nnodes < 0 && !got_input_file ) {
		printf("ERROR: number of nodes of input file not specified\n");
		help = 1; 
		exit(1); 
	} 		

	if ( help ) exit(1);

}    





 

