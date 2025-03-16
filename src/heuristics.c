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
    inst->best_sol = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (inst->best_sol == NULL){
        print_error("Memory allocation failed");
    }
    if (VERBOSE >= 50){
        printf("Nearest Neighbor calculations\n");
    }
    double best_nn_tour_cost = 1e20;  // Initialize the best tour cost for nearest neighbor

    for (int start = 0; start < inst->nnodes; start++) {
        double *tour = (double *) calloc(inst->nnodes+1, sizeof(double));  
        bool *visited = calloc(inst->nnodes, sizeof(bool));  
        if (!visited){
            print_error("Memory allocation error for visited");
        } 

        int current = start;  
        tour[0] = current + 1;  
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

            tour[i] = (double)nextNode + 1;  
            visited[nextNode] = true;  
            current = nextNode;  
        }
        tour[inst->nnodes] = tour[0];

        double cur_sol_cost = calculate_tour_cost(tour, inst);  
        if (cur_sol_cost < best_nn_tour_cost) {
            best_nn_tour_cost = cur_sol_cost;  
        }

        check_if_best_solution(tour, cur_sol_cost, inst);
        double time_now = second();
        if (time_now-inst->start_time > inst->time_limit){
            if (VERBOSE >= 50) printf("Time limit reached in NN\n");
            if (VERBOSE >= 30){
                printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
                if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
            }            
            free(tour);  
            free(visited);  
            return;
        }
        if (use_two_opt) two_opt(tour, inst);
        

        free(tour);
        free(visited);  
    }
    if (VERBOSE >= 30){
        printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
        if (use_two_opt) printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
    }
}

//Function to immpement variable neighborhood search(VNS) for TSP
void variable_neighborhood_search(instance *inst, double learning_rate, bool exponential_learning_rate) {
    double t1 = second();  // Start time
    double *tour = (double *) calloc(inst->nnodes+1, sizeof(double));

    if (tour == NULL){
        print_error("Memory allocation failed");
    }

    // Initialize the tour with the best solution
    memmove(tour, inst->best_sol, (inst->nnodes + 1) * sizeof(double));
    if (VERBOSE >= 50){
        printf("----------------------------------------------------------------------------------------------\n\n");
        printf("Variable Neighborhood Search calculations\n");
    }
    
    double best_tour_cost = calculate_tour_cost(tour, inst);
    double jumps_to_perform = 1.0;
    while(t1-inst->start_time < inst->time_limit) {
        two_opt(tour, inst);
        double running_tour_cost = calculate_tour_cost(tour, inst);
        if(running_tour_cost >= best_tour_cost) {
            for(int i = 0; i < (int)jumps_to_perform; i++) {
                three_opt(tour, inst);
            }
            if (exponential_learning_rate) jumps_to_perform*=learning_rate;
            else jumps_to_perform+=learning_rate;
        }
        else {
            best_tour_cost = running_tour_cost;
            if (VERBOSE >=50) fprintf(stdout, "Jumps performed to find better solution: %d\n", (int) jumps_to_perform);
            jumps_to_perform = 1;
        }
        t1 = second();
    }
    if (VERBOSE >= 30){
        printf("VNS FINAL COST: %lf\n", inst->best_sol_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");
    }
    free(tour);  
}

// Function to implement 2-opt heuristic for TSP that uses instance's best solution
void two_opt(double *solution, instance *inst) {
    double t1 = second();  // Start time
    double *tour = (double *) calloc(inst->nnodes+1, sizeof(double));  
    if (tour == NULL){
        print_error("Memory allocation failed");
    }

    // Initialize the tour with the best solution
    memmove(tour, solution, (inst->nnodes + 1) * sizeof(double));
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt INITIAL COST: %lf\n", calculate_tour_cost(tour, inst));
    }

    bool improved = true;  
    while(improved) {
        improved = false;
        for (int i = 0; i < inst->nnodes ; i++) {
            for (int j = i+1; j < inst->nnodes; j++) {
    
                // Calculate the new tour cost
                double cost_delta = inst->distances[(int)(tour[i]-1) * inst->nnodes + (int)(tour[j]-1)] + inst->distances[(int)(tour[i+1]-1) * inst->nnodes + (int)(tour[j+1]-1)] - 
                                    inst->distances[(int)(tour[i]-1) * inst->nnodes + (int)(tour[i+1]-1)] - inst->distances[(int)(tour[j]-1) * inst->nnodes + (int)(tour[j+1]-1)];

                // Check if the new tour is better
                if (cost_delta < -EPSILON) {
                    // Reverse the segment between i+1 and j in place
                    int start = i + 1;
                    int end = j;
                    while (start < end) {
                        double temp = tour[start];
                        tour[start] = tour[end];
                        tour[end] = temp;
                        start++;
                        end--;
                    }
                    improved = true;  
                    if (VERBOSE >= 80){
                        fprintf(stdout, "MODIFIED COST: %lf\n", calculate_tour_cost(tour, inst));
                    }
                    check_if_best_solution(tour, calculate_tour_cost(tour, inst), inst);

                    if (inst->best_sol_cost == calculate_tour_cost(tour, inst)) {
                        // Update time left only if a new best solution is found
                        double time_now = second();

                        // Check if the time limit has been reached
                        if (time_now-inst->start_time > inst->time_limit) {
                            if (VERBOSE >= 50) printf("Time limit reached in two_opt\n");
                            memmove(solution, tour, (inst->nnodes+1) * sizeof(double));
                            free(tour); 
                            return; 
                        }
                    }
                }
            }
        }
    }
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt FINAL COST: %lf\n", inst->best_sol_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");

    }
    memmove(solution, tour, (inst->nnodes+1) * sizeof(double));
    free(tour);  
}

// Function to implement 3-opt jump
void three_opt(double *solution, instance *inst) {
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
    double *new_tour = (double *) calloc(inst->nnodes+1, sizeof(double));
    if (new_tour == NULL){
        print_error("Memory allocation failed");
    }
    int count = 0;
    for (int l = 0; l <= i; l++) {
        new_tour[count] = solution[l];
        count++;
    }
    for (int l = j+1; l <= k; l++) {
        new_tour[count] = solution[l];
        count++;
    }
    for (int l = i+1; l <= j; l++) {
        new_tour[count] = solution[l];
        count++;
    }
    for (int l = k+1; l < inst->nnodes; l++) {
        new_tour[count] = solution[l];
        count++;
    }
    new_tour[inst->nnodes] = new_tour[0];

    memmove(solution, new_tour, (inst->nnodes+1) * sizeof(double));
    free(new_tour); 
}
