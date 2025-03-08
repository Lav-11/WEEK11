#ifndef TSP_FUNCTIONS_H
#define TSP_FUNCTIONS_H

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



// Structure for the TSP instance
typedef struct {
    int nnodes;                     // Number of nodes in the TSP problem
    double *xcoord;                 // Array of x coordinates for each node
    double *ycoord;                 // Array of y coordinates for each node
    int seed;                       // Seed used to generate the random instance (if applicable)
    double* distances;              // Distance matrix for the TSP instance
    double time_limit;				// overall time limit, in sec.s
    double time_left;               // time left, in sec.s
    char input_file[1000];          // The name of the input file (for debugging or reference)
    double *best_sol;               // Best known solution (tour) for the TSP instance
    double best_sol_cost;           // Cost of the best known solution
} instance;

// Function prototypes

// Function to print an error message and terminate the program
void print_error(const char *err);

// Function to calculate the Euclidean distance between two nodes i and j in the TSP instance
double dist(int i, int j, instance *inst);

// Function to calculate the distance matrix for the TSP instance
void calculate_distances(instance *inst);

// Function to save the solution in a file that will be used by gnuplot
void export_solution_for_gnuplot(const char *filename, const instance *inst);

// Function to save the solution in a PNG file using gnuplot
void png_solution_for_gnuplot(const char *input_filename, const char *output_filename);

// Function to check the feasibility of a tour
bool check_tour_feasability(const double *tour, instance *inst);

// Function to check if the current solution is better than the best solution found so far
void check_solution(double* tour, double cur_sol_cost, instance *inst);

// Function to update the best solution found so far
void update_best_solution(double* tour, double cur_sol_cost, instance *inst);

// Function to find the nearest neighbor tour for the TSP
void nearest_neighbor(instance *inst);

// Function to calculate the cost of a tour for the TSP instance
double calculate_tour_cost(const double *tour,  instance *inst);

// Function to implement the 2-opt heuristic for the TSP
void two_opt(double *solution, instance *inst);

#endif // TSP_FUNCTIONS_H
