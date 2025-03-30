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
#include "multithread_utils.h"

#include "tsp_utils.h"



void set_default_params(ConfigParams *params) {
    params->sequential_seed = 0;
    params->use_tabu_search = false;
    params->min_tenure_dimension_lower_bound = 0.1;
    params->min_tenure_dimension_higher_bound = 0.1;
    params->min_tenure_dimension_delta = 0.1;
    params->max_tenure_dimension_lower_bound = 0.5;
    params->max_tenure_dimension_higher_bound = 0.5;
    params->max_tenure_dimension_delta = 0.1;
    params->increase_ten_dim_rate_lower_bound = 1.0;
    params->increase_ten_dim_rate_higher_bound = 1.0;
    params->increase_ten_dim_rate_delta = 0.5;
    params->use_vns_search = false;
    params->learning_rate_lower_bound = 0.1;
    params->learning_rate_higher_bound = 0.1;
    params->learning_rate_delta = 3;
    params->max_jumps_lower_bound = 2;
    params->max_jumps_higher_bound = 3;
    params->max_jumps_delta = 1;
}


void *thread_function(void *arg) {
	ThreadData *data = (ThreadData *)arg;
	FILE *csv_file = fopen("../data/results.csv", "a");
	if (!csv_file) {
		print_error("Failed to open results.csv for writing");
	}

	if (data->is_tabu_search) {
		tabu_search(data->sol, data->inst->time_limit, data->inst, data->min_tenure, data->max_tenure, data->increase_rate);

		fprintf(csv_file, "tabu_%.2f_%.2f_%.3f,num_nodes_%d_seed_%d_time_limit_%.2f,%.2f\n",
			data->min_tenure, data->max_tenure, data->increase_rate,
			data->inst->nnodes, data->inst->seed, data->inst->time_limit, data->sol->tour_cost);

		free_instance(data->inst, true);

	} else if (data->is_vns_search) {
		variable_neighborhood_search(data->sol, data->inst->time_limit, data->inst, data->learning_rate, data->max_jumps);

		fprintf(csv_file, "vns_%.3f_%d,num_nodes_%d_seed_%d_time_limit_%.2f,%.2f\n",
			data->learning_rate, data->max_jumps,
			data->inst->nnodes, data->inst->seed, data->inst->time_limit, data->sol->tour_cost);
		
		free_instance(data->inst, true);
	}

	fclose(csv_file);
	pthread_exit(NULL);
}