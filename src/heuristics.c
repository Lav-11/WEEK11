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
void nearest_neighbor(instance *inst, bool use_two_opt) {
    inst->best_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (inst->best_sol->tour == NULL){
        print_error("Memory allocation failed");
    }
    if (VERBOSE >= 50){
        printf("Nearest Neighbor calculations\n");
    }
    double best_nn_tour_cost = 1e20;  // Store best tour cost for nearest neighbor

    for (int start = 0; start < inst->nnodes; start++) {

        solution *sol = (solution *) malloc(sizeof(solution));
        sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
        if (sol->tour == NULL){
            print_error("Memory allocation failed");
        }
        bool *visited = calloc(inst->nnodes, sizeof(bool));
        if (!visited){
            print_error("Memory allocation error for visited");
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

            if (nextNode == -1) print_error("Error constructing the tour");

            sol->tour[i] = (double)nextNode + 1;
            visited[nextNode] = true;
            current = nextNode;
        }
        sol->tour[inst->nnodes] = sol->tour[0];

        calculate_tour_cost(sol, inst);
        if (sol->tour_cost < best_nn_tour_cost) {
            best_nn_tour_cost = sol->tour_cost;
        }

        check_if_best_solution(sol, inst);
        double time_now = second();
        if (time_now-inst->start_time > inst->time_limit){
            if (VERBOSE >= 50) printf("Time limit reached in NN\n");
            if (VERBOSE >= 30){
                printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
                if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol->tour_cost);
            }
            free(sol);
            free(visited);
            return;
        }
        if (use_two_opt) two_opt(sol, inst);


        free(sol);
        free(visited);
    }
    if (VERBOSE >= 30){
        printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
        if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol->tour_cost);
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
//     double best_nn_tour_cost = 1e20;  // Initialize the best tour cost for nearest neighbor

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
//         if (cur_sol_cost < best_nn_tour_cost) {
//             best_nn_tour_cost = cur_sol_cost;
//         }

//         check_if_best_solution(tour, cur_sol_cost, inst);
//         double time_now = second();
//         if (time_now-inst->start_time > inst->time_limit){
//             if (VERBOSE >= 50) printf("Time limit reached in NN\n");
//             if (VERBOSE >= 30){
//                 printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
//                 if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
//             }
//             free(tour);
//             return;
//         }
//         if (use_two_opt) two_opt(tour, inst);


//         free(tour);
//     if (VERBOSE >= 30){
//         printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
//         if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
//     }
// }

//Function to immpement variable neighborhood search(VNS) for TSP
void variable_neighborhood_search(instance *inst, double learning_rate, int max_jumps) {
    double t1 = second();  // Start time
    solution *sol = (solution *) malloc(sizeof(solution));
    sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (sol->tour == NULL){
        print_error("Memory allocation failed");
    }

    // Initialize the tour with the best solution
    memcpy(sol->tour, inst->best_sol->tour, (inst->nnodes + 1) * sizeof(double));
    if (VERBOSE >= 50){
        printf("----------------------------------------------------------------------------------------------\n\n");
        printf("Variable Neighborhood Search calculations\n");
    }

    solution *best_sol = (solution *) malloc(sizeof(solution));
    best_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (best_sol->tour == NULL){
        print_error("Memory allocation failed");
    }
    memcpy(best_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
    double best_tour_cost = inst->best_sol->tour_cost;
    double jumps_to_perform = 1.0;
    bool skip_two_opt = false;
    while(t1-inst->start_time < inst->time_limit) {
        if (!skip_two_opt) two_opt(sol, inst);
        double running_tour_cost = inst->best_sol->tour_cost;
        if(running_tour_cost >= best_tour_cost-EPSILON) {
            memcpy(sol->tour, best_sol->tour, (inst->nnodes + 1) * sizeof(double));
            for(int i = 0; i < ((int)jumps_to_perform % max_jumps); i++) {
                three_opt(sol, inst);
            }
            jumps_to_perform+=learning_rate;
            skip_two_opt = false;
        }
        else {
            memcpy(best_sol->tour, sol->tour, (inst->nnodes + 1) * sizeof(double));
            best_tour_cost = running_tour_cost;
            if (VERBOSE >=50) fprintf(stdout, "Jumps performed to find better solution:  %d\n", (int)jumps_to_perform % max_jumps);
            jumps_to_perform = 1.0;
            skip_two_opt = true;
        }
        t1 = second();
    }
    if (VERBOSE >= 30){
        printf("VNS FINAL COST: %lf\n", inst->best_sol->tour_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");
    }
    free(sol);
    free(best_sol);
}

// Function to implement 2-opt heuristic for TSP that uses instance's best solution
void two_opt(solution *sol, instance *inst) {
    double t1 = second();  // Start time

    solution *temp_sol = (solution *) malloc(sizeof(solution));
    temp_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (temp_sol->tour == NULL){
        print_error("Memory allocation failed");
    }

    // Initialize the tour with the best solution
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

                // Check if the new tour is better
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
                    check_if_best_solution(temp_sol, inst);

                    if (inst->best_sol->tour_cost == temp_sol->tour_cost) {
                        // Update time left only if a new best solution is found
                        double time_now = second();

                        // Check if the time limit has been reached
                        if (time_now-inst->start_time > inst->time_limit) {
                            if (VERBOSE >= 50) printf("Time limit reached in two_opt\n");
                            memcpy(sol->tour, temp_sol->tour, (inst->nnodes + 1) * sizeof(double));
                            free(temp_sol);
                            return;
                        }
                    }
                }
            }
        }
    }
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt FINAL COST: %lf\n", inst->best_sol->tour_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");

    }
    memcpy(sol->tour, temp_sol->tour, (inst->nnodes + 1) * sizeof(double));
    free(temp_sol);
}

// Function to implement 3-opt jump
void three_opt(solution *sol, instance *inst) {
    if (VERBOSE >= 60){
        printf("3-opt jump called");
    }

    // Generate randomly 3 integers i, j, k such that i < j < k < nnodes
    int i, j, k;
    srand(time(NULL)); // Seed the random number generator with the current time
    do {
        i = rand() % inst->nnodes;
        j = rand() % inst->nnodes;
        k = rand() % inst->nnodes;
    } while (!(i < j && j < k && i != j && j != k && i != k));

    // Perform the swap between the block between i+1 and j and the block between j+1 and k
    solution *new_sol = (solution *) malloc(sizeof(solution));
    new_sol->tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (new_sol->tour == NULL){
        print_error("Memory allocation failed");
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
    free(new_sol->tour);
    free(new_sol);
}
