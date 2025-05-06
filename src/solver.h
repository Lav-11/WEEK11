#ifndef SOLVER_H
#define SOLVER_H

#include "tsp_utils.h"
#include "heuristics.h"
#include <pthread.h>

typedef struct {
	int sequential_seed;


	//HEURISTICS
	bool use_tabu_search;

    double min_tenure_dimension_lower_bound;
	double min_tenure_dimension_higher_bound;
	double min_tenure_dimension_delta;

    double max_tenure_dimension_lower_bound;
	double max_tenure_dimension_higher_bound;
	double max_tenure_dimension_delta;

    bool use_vns_search;

    double learning_rate_lower_bound;
	double learning_rate_higher_bound;
	double learning_rate_delta;

    int max_jumps_lower_bound;
	int max_jumps_higher_bound;
	int max_jumps_delta;


	//CPLEX
	bool warmstart;
	bool posting;

	bool use_benders;

	bool use_bc;
	bool fractional;

	bool use_hardfixing;
	double hf_prob_lower_bound;
	double hf_prob_higher_bound;
	double hf_prob_delta;
	int hf_num_of_run;

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

typedef struct {
    double best_cost;
    char best_params[256]; 
    pthread_mutex_t best_cost_mutex; 
} BestResult;

extern BestResult global_best_result;

void set_default_params(ConfigParams *params);

void *thread_function(void *arg);

#endif // MULTITHREAD_UTILS_H