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

// Function to implement the GRASP heuristic for the TSP
void grasp(instance *inst, bool use_two_opt, double deviating_probability, bool prob_proportional_to_cost);

// Function to implement the variable neighborhood search for the TSP
void variable_neighborhood_search(instance *inst, double learning_rate, int max_jumps);

// Function to implement the tabu search for the TSP
void tabu_search(instance *inst, double tenure_dimension);

// Function to implement the 2-opt heuristic for the TSP
void two_opt(solution *sol, instance *inst);

// Function to implement the 3-opt heuristic for the TSP
void three_opt(solution *sol, instance *inst);

#endif // HEURISTICS_H