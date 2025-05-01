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
#include "solver.h"
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
    params->use_vns_search = false;
    params->learning_rate_lower_bound = 0.1;
    params->learning_rate_higher_bound = 0.1;
    params->learning_rate_delta = 3;
    params->max_jumps_lower_bound = 2;
    params->max_jumps_higher_bound = 3;
    params->max_jumps_delta = 1;
    params->use_benders = false;
    params->use_bc = false;
    params->posting = false;
    params->fractional = false;
    params->warmstart = false;
}


pthread_mutex_t file_mutex = PTHREAD_MUTEX_INITIALIZER;

void *thread_function(void *arg) {
    ThreadData *data = (ThreadData *)arg;

    double current_cost = 0.0;

    if (data->is_tabu_search) {
        tabu_search(data->sol, data->inst->time_limit, data->inst, data->min_tenure, data->max_tenure);

        current_cost = data->sol->tour_cost;

		pthread_mutex_lock(&file_mutex);
		FILE *csv_file = fopen("../data/results.csv", "a");
		if (!csv_file) {
			print_error("Failed to open results.csv for writing");
			pthread_mutex_unlock(&file_mutex); // Unlock before exiting
			pthread_exit(NULL);
		}

        fprintf(csv_file, "tabu_2_%.2f_%.2f,num_nodes_%d_seed_%d_time_limit_%.2f,%.2f\n",
            data->min_tenure, data->max_tenure,
            data->inst->nnodes, data->inst->seed, data->inst->time_limit, current_cost);
		
		fclose(csv_file);
	    pthread_mutex_unlock(&file_mutex);
        free_instance(data->inst, true);

    } else if (data->is_vns_search) {
        variable_neighborhood_search(data->sol, data->inst->time_limit, data->inst, data->learning_rate, data->max_jumps);

        current_cost = data->sol->tour_cost;

		pthread_mutex_lock(&file_mutex);
		FILE *csv_file = fopen("../data/results.csv", "a");
		if (!csv_file) {
			print_error("Failed to open results.csv for writing");
			pthread_mutex_unlock(&file_mutex); // Unlock before exiting
			pthread_exit(NULL);
		}

        fprintf(csv_file, "vns_%.3f_%d,num_nodes_%d_seed_%d_time_limit_%.2f,%.2f\n",
            data->learning_rate, data->max_jumps,
            data->inst->nnodes, data->inst->seed, data->inst->time_limit, current_cost);
        
		fclose(csv_file);
		pthread_mutex_unlock(&file_mutex);
        free_instance(data->inst, true);
    }


    // Update the global best result
    pthread_mutex_lock(&global_best_result.best_cost_mutex);
    if (current_cost < global_best_result.best_cost) {
        global_best_result.best_cost = current_cost;
        if (data->is_tabu_search) {
            snprintf(global_best_result.best_params, sizeof(global_best_result.best_params),
                     "tabu_%.2f_%.2f", data->min_tenure, data->max_tenure);
        } else if (data->is_vns_search) {
            snprintf(global_best_result.best_params, sizeof(global_best_result.best_params),
                     "vns_%.3f_%d", data->learning_rate, data->max_jumps);
        }
    }
    pthread_mutex_unlock(&global_best_result.best_cost_mutex);

    pthread_exit(NULL);
}