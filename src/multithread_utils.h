#ifndef MULTITHREAD_UTILS_H
#define MULTITHREAD_UTILS_H

#include "tsp_utils.h"
#include "heuristics.h"
#include <pthread.h>

typedef struct {
	int sequential_seed;


	bool use_tabu_search;

    double min_tenure_dimension_lower_bound;
	double min_tenure_dimension_higher_bound;
	double min_tenure_dimension_delta;

    double max_tenure_dimension_lower_bound;
	double max_tenure_dimension_higher_bound;
	double max_tenure_dimension_delta;

    double increase_ten_dim_rate_lower_bound;
	double increase_ten_dim_rate_higher_bound;
	double increase_ten_dim_rate_delta;


    bool use_vns_search;

    double learning_rate_lower_bound;
	double learning_rate_higher_bound;
	double learning_rate_delta;

    int max_jumps_lower_bound;
	int max_jumps_higher_bound;
	int max_jumps_delta;
} ConfigParams;

typedef struct {
	instance *inst;
	ConfigParams params;
	solution *sol;
	double min_tenure;
	double max_tenure;
	double increase_rate;
	double learning_rate;
	int max_jumps;
	bool is_tabu_search;
	bool is_vns_search;
} ThreadData;

// Function to set default parameters for the configuration
void set_default_params(ConfigParams *params);

// Thread function to be executed by each thread
void *thread_function(void *arg);

#endif // MULTITHREAD_UTILS_H