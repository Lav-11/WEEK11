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

	//General calls for initial calculations and input reading
	parse_command_line(argc,argv, &inst);  
	read_input(&inst);  
	calculate_distances(&inst);
	
	//Nearest neighbor heuristic and gnuplot output
	nearest_neighbor(&inst);
	export_solution_for_gnuplot("../data/solution.dat", &inst);
	png_solution_for_gnuplot("../data/solution.dat", "../data/solution.png");

	//2-opt heuristic and gnuplot output
	export_solution_for_gnuplot("../data/final_solution.dat", &inst);
	png_solution_for_gnuplot("../data/final_solution.dat", "../data/final_solution.png");

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
	}
	else {
		// Set the number of nodes (random value between 50 and 100)
		inst->xcoord = malloc(inst->nnodes * sizeof(double));  
		inst->ycoord = malloc(inst->nnodes * sizeof(double));  
		strncpy(inst->input_file, "randomly generated", 1000);  
		int seed = inst->seed;  
		srand(seed);  

		// Generate random coordinates between 0 and 9999
		for (int i = 0; i < inst->nnodes; i++) {
			inst->xcoord[i] = (double)((double)rand() / RAND_MAX)*10000;  
			inst->ycoord[i] = (double)((double)rand() / RAND_MAX)*10000; 
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
	inst->time_limit = 60; 
	inst->time_left = inst->time_limit;
	inst->best_sol_cost = 1e+20;
	int got_input_file = 0;

    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); got_input_file=1; continue; } 	// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); got_input_file=1; continue; } 	// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); got_input_file=1; continue; } 		// input file
		if ( strcmp(argv[i],"-tl") == 0 ) { inst->time_limit = atof(argv[++i]); continue; }							// time limit
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->seed = abs(atoi(argv[++i])); continue; } 						// random seed
		if ( strcmp(argv[i],"-num_nodes") == 0 ) { inst->nnodes = atoi(argv[++i]); continue; } 						// random seed
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 												// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 												// help
		help = 1;
    }      

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\navailable parameters --------------------------------------------------\n");
		printf("-file: %s\n", inst->input_file); 
		printf("-time_limit: %lf\n", inst->time_limit);
		printf("-seed %d\n", inst->seed); 
		printf("-number of random nodes (-1 if a file is provided): %d\n", inst->nnodes); 
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





 

