#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "tsp_utils.h"
#include "heuristics.h"
#include "chrono.h"


// Function to find the nearest neighbor tour for the TSP
void nearest_neighbor(instance *inst, int starting_node, bool use_two_opt) {
    
    if (inst->best_sol == NULL) {
        inst->best_sol = (solution *)malloc(sizeof(solution));
        if (!inst->best_sol) {
            print_error("Memory allocation failed for best_sol in nearest_neighbor");
        }
        inst->best_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
        if (inst->best_sol->tour == NULL){
            print_error("Memory allocation failed in nearest_neighbor");
        }
        inst->best_sol->tour_cost = 1e20;
    }

    int start = starting_node % inst->nnodes;
    solution *sol = (solution *) malloc(sizeof(solution));
    sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (sol->tour == NULL){
        print_error("Memory allocation failed for tour in nearest_neighbor");
    }
    bool *visited = calloc(inst->nnodes, sizeof(bool));
    if (!visited){
        print_error("Memory allocation error for visited in nearest_neighbor");
    }

    int current = start;
    sol->tour[0] = current + 1;
    visited[current] = true;
    for (int i = 1; i < inst->nnodes; i++) {
        double minDist = 1e20;
        int nextNode = -1;
        for (int j = 0; j < inst->nnodes; j++) {
            if (!visited[j]) {
                double d = inst->distances[current * inst->nnodes + j];
                if (d < minDist) {
                    minDist = d;
                    nextNode = j;
                }
            }
        }

        if (nextNode == -1) print_error("Error constructing the tour in nearest_neighbor");

        sol->tour[i] = (double)nextNode + 1;
        visited[nextNode] = true;
        current = nextNode;
    }
    sol->tour[inst->nnodes] = sol->tour[0];
    calculate_tour_cost(sol, inst);
    check_if_best_solution(sol, inst);
    double time_now = second();
    if (time_now-inst->start_time > inst->time_limit){
        if (VERBOSE >= 50) printf("Time limit reached in NN\n");
        if (VERBOSE >= 50){
            printf("NEAREST NEIGHBOR FINAL COST: %lf\n", inst->best_sol->tour_cost);
        }
        free_solution(sol);
        free(visited);
        return;
    }
    if (use_two_opt) {
        two_opt(sol, 1e20, inst);
        check_if_best_solution(sol, inst);
    }


    free_solution(sol);
    free(visited);
    
    if (VERBOSE >= 50){
        printf("NEAREST NEIGHBOR FINAL COST: %lf\n", inst->best_sol->tour_cost);
    }
}

// Function to calculate the best nearest_neigbbor tour for the TSP
void best_nearest_neighbor(instance *inst, bool use_two_opt) {
    inst->best_sol = (solution *)malloc(sizeof(solution));
	if (!inst->best_sol) {
	print_error("Memory allocation failed for best_sol");
	};
    inst->best_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (inst->best_sol->tour == NULL){
        print_error("Memory allocation failed in best_nearest_neighbor");
    }
    inst->best_sol->tour_cost = 1e+20;

    if (VERBOSE >= 50){
        printf("Nearest Neighbor calculations\n");
    }
    for (int i = 0; i < inst->nnodes; i++) {
        nearest_neighbor(inst, i, use_two_opt);
        double time_now = second();
        if (time_now-inst->start_time > inst->time_limit){
            if (VERBOSE >= 20) printf("Time limit reached in best_nearest_neighbor\n");
            break;
        }
    }
}

// Function to implement GRASP for TSP
// void grasp(instance *inst, bool use_two_opt, double deviating_probability, bool prob_proportional_to_cost) {
//     inst->best_sol = (double *) calloc(inst->nnodes+1, sizeof(double));
//     if (inst->best_sol == NULL){
//         print_error("Memory allocation failed");
//     }
//     double *tour = (double *) calloc(inst->nnodes+1, sizeof(double));
//     if (tour == NULL){
//         print_error("Memory allocation failed");
//     }
//     if (VERBOSE >= 50){
//         printf("Nearest Neighbor calculations\n");
//     }
//     double nn_tour_cost = 1e20;  // Initialize the best tour cost for nearest neighbor

//     if (!prob_proportional_to_cost) {
//         for (int start = 0; start < inst->nnodes; start++) {
//             bool *visited = calloc(inst->nnodes, sizeof(bool));
//             if (!visited){
//                 print_error("Memory allocation error for visited");
//             }

//             int current = start;
//             tour[0] = current + 1;
//             visited[current] = true;
//             for (int i = 1; i < inst->nnodes; i++) {
//                 double minDist[4] = {1e20, 1e20, 1e20, 1e20};
//                 int nextNode[4] = {-1, -1, -1, -1};
//                 for (int j = 0; j < inst->nnodes; j++) {
//                     if (!visited[j]) {
//                     double d = inst->distances[current * inst->nnodes + j];
//                     if (d < minDist[0]) {
//                         minDist[3] = minDist[2];
//                         nextNode[3] = nextNode[2];
//                         minDist[2] = minDist[1];
//                         nextNode[2] = nextNode[1];
//                         minDist[1] = minDist[0];
//                         nextNode[0] = j;
//                         minDist[0] = d;
//                     } else if (d < minDist[1]) {
//                         minDist[3] = minDist[2];
//                         nextNode[3] = nextNode[2];
//                         minDist[2] = minDist[1];
//                         nextNode[2] = nextNode[1];
//                         minDist[1] = d;
//                         nextNode[1] = j;
//                     } else if (d < minDist[2]) {
//                         minDist[3] = minDist[2];
//                         nextNode[3] = nextNode[2];
//                         minDist[2] = d;
//                         nextNode[2] = j;
//                     } else if (d < minDist[3]) {
//                         minDist[3] = d;
//                         nextNode[3] = j;
//                     }
//                     }
//                 }
//                 int chose_node_index = 0; // Default to the first node
//                 if (random01(&inst->seed) < deviating_probability && nextNode[3] != -1) {
//                     chose_node_index = (rand() % 3) + 1; // Randomly choose 1, 2, or 3
//                 }
//                 int chosen_node = nextNode[chose_node_index];
//                 printf("Chose node %d\n", chosen_node);
//                 if (chosen_node == -1) print_error("Error constructing the tour here");

//                 tour[i] = (double)chosen_node + 1;
//                 visited[chosen_node] = true;
//                 current = chosen_node;
//             }
//         }
//     }
//         tour[inst->nnodes] = tour[0];

//         double cur_sol_cost = calculate_tour_cost(tour, inst);
//         if (cur_sol_cost < nn_tour_cost) {
//             nn_tour_cost = cur_sol_cost;
//         }

//         check_if_best_solution(tour, cur_sol_cost, inst);
//         double time_now = second();
//         if (time_now-inst->start_time > inst->time_limit){
//             if (VERBOSE >= 50) printf("Time limit reached in NN\n");
//             if (VERBOSE >= 30){
//                 printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", nn_tour_cost);
//                 if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
//             }
//             free(tour);
//             return;
//         }
//         if (use_two_opt) two_opt(tour, inst);


//         free(tour);
//     if (VERBOSE >= 30){
//         printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", nn_tour_cost);
//         if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
//     }
// }


//Function to implement variable neighborhood search(VNS) for TSP
void variable_neighborhood_search(solution *sol, double timelimit, const instance *inst, double learning_rate, int max_jumps) {

    if (inst->nnodes < 3) {
        print_error("VNS is not applicable for instances with less than 3 nodes.");
    }

    double t_start = second();  
    double time_now = second();

    // Open a file to save the costs
    FILE *cost_file = fopen("../data/vns_costs.txt", "w");
    if (cost_file == NULL) {
        print_error("Failed to open file for writing costs in VNS");
    }

    // Initialize the tour with the best solution
    if (VERBOSE >= 50){
        printf("----------------------------------------------------------------------------------------------\n\n");
        printf("Variable Neighborhood Search calculations\n");
    }

    solution *best_sol = (solution *) malloc(sizeof(solution));
    best_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (best_sol->tour == NULL){
        print_error("Memory allocation failed in VNS for best solution tour");
    }
    memcpy(best_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
    best_sol->tour_cost = sol->tour_cost;

    double jumps_to_perform = 1.0;
    bool skip_two_opt = false;

    solution *temp_sol = (solution *) malloc(sizeof(solution));
    temp_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (temp_sol->tour == NULL){
        print_error("Memory allocation failed in VNS for temp solution tour");
    }
    memcpy(temp_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
    temp_sol->tour_cost = sol->tour_cost;

    while(time_now - t_start < timelimit) {
        if (!skip_two_opt) {
            two_opt(temp_sol, timelimit-(time_now-t_start), inst);
            // Save the cost after 2-opt
            fprintf(cost_file, "%lf\n", temp_sol->tour_cost);
        }
        double running_tour_cost = temp_sol->tour_cost;
        if(running_tour_cost >= best_sol->tour_cost-EPSILON) {
            memcpy(temp_sol->tour, best_sol->tour, (inst->nnodes + 1) * sizeof(double));
            for(int i = 0; i < ((int)jumps_to_perform % max_jumps); i++) {
                three_opt(temp_sol, inst);
            }
            // Save the cost after 3-opt
            fprintf(cost_file, "%lf\n", temp_sol->tour_cost);
            jumps_to_perform+=learning_rate;
            skip_two_opt = false;
        }
        else {
            memcpy(best_sol->tour, temp_sol->tour, (inst->nnodes + 1) * sizeof(double));
            best_sol->tour_cost = running_tour_cost;

            if (VERBOSE >=60) fprintf(stdout, "Jumps performed to find better solution:  %d\n", (int)jumps_to_perform % max_jumps);
            jumps_to_perform = 1.0;
            skip_two_opt = true;
        }
        time_now = second();
    }
    if (VERBOSE >= 50){
        calculate_tour_cost(best_sol, inst);
        printf("VNS FINAL COST: %lf\n", best_sol->tour_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");
    }
    memcpy(sol->tour, best_sol->tour, (inst->nnodes + 1) * sizeof(double));
    sol->tour_cost = best_sol->tour_cost;

    free_solution(best_sol);
    free_solution(temp_sol);

    // Close the file
    fclose(cost_file);
}

void tabu_search(solution *sol, double timelimit, const instance *inst, double min_tenure_dimension, double max_tenure_dimension, double increase_ten_dim_rate) {

    if (inst->nnodes < 10) {
        print_error("Tabu search functioning is at risk for instances with less than 10 nodes.");
    }
    if (max_tenure_dimension < min_tenure_dimension-EPSILON) {
        print_error("Max tenure dimension must be greater than min tenure dimension.");
    }


    double t_start = second();  // Start time
    double time_now = second();

    if (VERBOSE >= 50){
        printf("----------------------------------------------------------------------------------------------\n\n");
        printf("Tabu Search calculations\n");
    }
    // Tabu list initialization as a 2D matrix
    int max_tenure_length = (int)(max_tenure_dimension * inst->nnodes);
    int min_tenure_length = (int)(min_tenure_dimension * inst->nnodes);
    int current_tenure_length = min_tenure_length;
    bool is_tenure_increasing = true;
    int iteration_count = 0;
    int **tenure_matrix = (int **)calloc(inst->nnodes, sizeof(int *));
    for (int i = 0; i < inst->nnodes; i++) {
        tenure_matrix[i] = (int *)calloc(inst->nnodes, sizeof(int));
        if (tenure_matrix[i] == NULL) {
            print_error("Memory allocation failed for tabu matrix");
        }
    }

    // Solution to store the best outcome found so far
    solution *best_sol = (solution *) malloc(sizeof(solution));
    best_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (best_sol->tour == NULL){
        print_error("Memory allocation failed for best solution tour in tabu_search");
    }
    memcpy(best_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
    best_sol->tour_cost = sol->tour_cost;

    // Open a file to save the costs
    FILE *cost_file = fopen("../data/tabu_costs.txt", "w");
    if (cost_file == NULL) {
        print_error("Failed to open file for writing costs in tabu_search");
    }

    while (time_now - t_start < timelimit) {
        double best_neighbor_cost = 1e20;
        int best_i = -1, best_j = -1;

        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 1; j < inst->nnodes; j++) {
                int i1 = (int)(sol->tour[i] - 1);
                int i2 = (int)(sol->tour[i + 1] - 1);
                int j1 = (int)(sol->tour[j] - 1);
                int j2 = (int)(sol->tour[j + 1] - 1);

                bool is_tabu = (tenure_matrix[i1][j1] > 0 || tenure_matrix[i2][j2] > 0 || tenure_matrix[j1][i1] > 0 || tenure_matrix[j2][i2] > 0);

                double cost_delta = inst->distances[i1 * inst->nnodes + j1] +
                                    inst->distances[i2 * inst->nnodes + j2] -
                                    inst->distances[i1 * inst->nnodes + i2] -
                                    inst->distances[j1 * inst->nnodes + j2];

                // Accept the move if it's not tabu or if it improves the best solution. It must also be different from the current solution.
                if (sol->tour_cost + cost_delta < best_neighbor_cost - EPSILON && !is_tabu && abs(cost_delta) > EPSILON) {
                    best_neighbor_cost = sol->tour_cost + cost_delta;
                    best_i = i;
                    best_j = j;
                }
            }
        }

        if (best_i != -1 && best_j != -1) {
            // Update the tabu matrix
            int i1 = (int)(sol->tour[best_i] - 1);
            int i2 = (int)(sol->tour[best_i + 1] - 1);
            int j1 = (int)(sol->tour[best_j] - 1);
            int j2 = (int)(sol->tour[best_j + 1] - 1);

            tenure_matrix[i1][i2] = current_tenure_length;
            tenure_matrix[j1][j2] = current_tenure_length;

            // Reverse the segment between best_i+1 and best_j
            int start = best_i + 1;
            int end = best_j;
            while (start < end) {
                double temp = sol->tour[start];
                sol->tour[start] = sol->tour[end];
                sol->tour[end] = temp;
                start++;
                end--;
            }

            // Update the solution cost
            calculate_tour_cost(sol, inst);
            fprintf(cost_file, "%lf\n", sol->tour_cost);
            if (sol->tour_cost < best_sol->tour_cost - EPSILON) {
                memcpy(best_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
                best_sol->tour_cost = sol->tour_cost;
            }

        } else {
            printf("Error occured in tabu_search in finding best i, j \n");
        }

        // Decrease tenure values by 1
        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = 0; j < inst->nnodes; j++) {
                if (tenure_matrix[i][j] > 0) {
                    tenure_matrix[i][j]--;
                }
            }
        }

        // Adjust the tabu list length based on increase_ten_dim_rate
        iteration_count++;
        if (iteration_count >= (int)(1.0 / increase_ten_dim_rate)) {
            iteration_count = 0;
            if (is_tenure_increasing) {
                current_tenure_length++;
                if (current_tenure_length > max_tenure_length) {
                    is_tenure_increasing = false;
                }
            } else {
                current_tenure_length--;
                if (current_tenure_length < min_tenure_length) {
                    is_tenure_increasing = true;
                }
            }
        }

        time_now = second();
    }

    if (VERBOSE >= 50) {
        printf("TABU SEARCH FINAL COST: %lf\n", best_sol->tour_cost);
    }

    for (int i = 0; i < inst->nnodes; i++) {
        free(tenure_matrix[i]);
    }
    free(tenure_matrix);
    free_solution(best_sol);
    fclose(cost_file);
}


// Function to implement 2-opt heuristic for TSP that uses instance's best solution
void two_opt(solution *sol, double timelimit, const instance *inst) {
    double t_start = second();  
    double time_now = second();

    solution *temp_sol = (solution *) malloc(sizeof(solution));
    temp_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (temp_sol->tour == NULL){
        print_error("Memory allocation failed in 2-opt for temp solution tour");
    }

    memcpy(temp_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt INITIAL COST: %lf\n", sol->tour_cost);
    }

    bool improved = true;
    while(improved) {
        improved = false;
        for (int i = 0; i < inst->nnodes ; i++) {
            for (int j = i+1; j < inst->nnodes; j++) {

                // Calculate the new tour cost
                double cost_delta = inst->distances[(int)(temp_sol->tour[i]-1) * inst->nnodes + (int)(temp_sol->tour[j]-1)] +
                                    inst->distances[(int)(temp_sol->tour[i+1]-1) * inst->nnodes + (int)(temp_sol->tour[j+1]-1)] -
                                    inst->distances[(int)(temp_sol->tour[i]-1) * inst->nnodes + (int)(temp_sol->tour[i+1]-1)] -
                                    inst->distances[(int)(temp_sol->tour[j]-1) * inst->nnodes + (int)(temp_sol->tour[j+1]-1)];

                if (cost_delta < -EPSILON) {
                    // Reverse the segment between i+1 and j in place
                    int start = i + 1;
                    int end = j;
                    while (start < end) {
                        double temp = temp_sol->tour[start];
                        temp_sol->tour[start] = temp_sol->tour[end];
                        temp_sol->tour[end] = temp;
                        start++;
                        end--;
                    }
                    improved = true;
                    if (VERBOSE >= 80){
                        fprintf(stdout, "MODIFIED COST: %lf\n", temp_sol->tour_cost);
                    }
                    calculate_tour_cost(temp_sol, inst);

                    // Update time left only if a new best solution is found
                    double time_now = second();

                    if (time_now-t_start > timelimit) {
                        if (VERBOSE >= 50) printf("Time limit reached in two_opt\n");
                        memcpy(sol->tour, temp_sol->tour, (inst->nnodes + 1) * sizeof(double));
                        calculate_tour_cost(sol, inst);
                        free(temp_sol);
                        return;
                    }
                }
            }
        }
    }
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt FINAL COST: %lf\n", temp_sol->tour_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");

    }
    memcpy(sol->tour, temp_sol->tour, (inst->nnodes + 1) * sizeof(double));
    calculate_tour_cost(sol, inst);
    free_solution(temp_sol);
}

// Function to implement 3-opt jump
void three_opt(solution *sol, const instance *inst) {
    if (VERBOSE >= 60){
        printf("3-opt jump called");
    }

    // Generate randomly 3 integers i, j, k such that i < j < k < nnodes
    int i, j, k;
    srand(time(NULL) + rand()); 
    do {
        i = rand() % inst->nnodes;
        j = rand() % inst->nnodes;
        k = rand() % inst->nnodes;
    } while (!(i < j && j < k && i != j && j != k && i != k));

    // Perform the swap between the block between i+1 and j and the block between j+1 and k
    solution *new_sol = (solution *) malloc(sizeof(solution));
    new_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (new_sol->tour == NULL){
        print_error("Memory allocation failed in 3-opt for new solution tour");
    }

    int count = 0;
    for (int l = 0; l <= i; l++) {
        new_sol->tour[count] = sol->tour[l];
        count++;
    }
    for (int l = j+1; l <= k; l++) {
        new_sol->tour[count] = sol->tour[l];
        count++;
    }
    for (int l = i+1; l <= j; l++) {
        new_sol->tour[count] = sol->tour[l];
        count++;
    }
    for (int l = k+1; l < inst->nnodes; l++) {
        new_sol->tour[count] = sol->tour[l];
        count++;
    }
    new_sol->tour[inst->nnodes] = new_sol->tour[0];

    memcpy(sol->tour, new_sol->tour, (inst->nnodes+1) * sizeof(double));
    calculate_tour_cost(sol, inst);
    free_solution(new_sol);
}
