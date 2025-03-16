#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <stdbool.h>  

//#include <cplex.h>  
//#include <pthread.h>  

// Function prototypes

// Function to find the nearest neighbor tour for the TSP
void nearest_neighbor(instance *inst, bool use_two_opt);

// Function to implement the variable neighborhood search for the TSP
void variable_neighborhood_search(instance *inst, double learning_rate, bool exponential_learning_rate);

// Function to implement the 2-opt heuristic for the TSP
void two_opt(double *solution, instance *inst);

// Function to implement the 3-opt heuristic for the TSP
void three_opt(double *solution, instance *inst);

#endif // HEURISTICS_H