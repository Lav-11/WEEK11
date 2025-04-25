#ifndef HEURISTICS_H
#define HEURISTICS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <stdbool.h>  
#include <pthread.h>  
#include "tsp_utils.h"
#include "multithread_utils.h"


// Function to find the nearest neighbor tour for the TSP
void nearest_neighbor(instance *inst, int starting_node, bool use_two_opt);

// Function to calculate the best nearest neighbor tour for the TSP
void best_nearest_neighbor(instance *inst, bool use_two_opt);

// Function to implement the GRASP heuristic for the TSP
void grasp(double timelimit, instance *inst, double deviating_probability);

// Function to implement the variable neighborhood search for the TSP
void variable_neighborhood_search(solution *sol, double timelimit, const  instance *inst, double learning_rate, int max_jumps);

// Function to implement the tabu search for the TSP
void tabu_search(solution *sol, double timelimit, const instance *inst, double min_tenure_dimension, double max_tenure_dimension);

// Function to implement the 2-opt heuristic for the TSP
void two_opt(solution *sol, double timelimit, const instance *inst);

// Function to implement the 3-opt heuristic for the TSP
void three_opt(solution *sol, const instance *inst);

#endif // HEURISTICS_H