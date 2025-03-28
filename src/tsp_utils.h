#ifndef TSP_UTILS_H
#define TSP_UTILS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <stdbool.h>  

//#include <cplex.h>  
//#include <pthread.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ



typedef struct {
    double *tour;
    double tour_cost;
} solution;


// Structure for the TSP instance
typedef struct {
    int nnodes;                     // Number of nodes in the TSP problem
    double *xcoord;                 // Array of x coordinates for each node
    double *ycoord;                 // Array of y coordinates for each node
    int seed;                       // Seed used to generate the random instance (if applicable)
    double* distances;              // Distance matrix for the TSP instance
    double time_limit;				// overall time limit, in sec.s
    double start_time;               // time left, in sec.s
    char input_file[1000];          // The name of the input file (for debugging or reference)
    solution *best_sol;             // Best known solution  for the TSP instance
} instance;


// Function prototypes

// Function to print an error message and terminate the program
void print_error(const char *err);

// Function to generate a random number between 0 and 1
double random01(unsigned int *seed);

// Function to calculate the Euclidean distance between two nodes i and j in the TSP instance
double dist(int i, int j, instance *inst);

// Function to calculate the distance matrix for the TSP instance
void calculate_distances(instance *inst);

// Function to save the solution in a PNG file using gnuplot
void png_solution_for_gnuplot(solution *sol, const bool save_dat_file, char *output_filename, instance *inst);

// Function to check the feasibility of a tour
bool check_tour_feasability(solution *sol, instance *inst);

// Function to check if the current solution is better than the best solution found so far
void check_if_best_solution(solution *sol, instance *inst);

// Function to update the best solution found so far
void update_best_solution(solution *sol, instance *inst);

// Function to calculate the cost of a tour for the TSP instance
void calculate_tour_cost(solution *sol,  const instance *inst);

// Function to free the memory of a solution struct
void free_solution(solution *sol);

// Function to read a file containing costs and plot a graph of it using gnuplot
void plot_costs(char *input_filename, char *output_filename);

#endif // TSP_UTILS_H
