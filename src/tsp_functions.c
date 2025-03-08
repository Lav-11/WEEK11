#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "tsp_functions.h"
#include "chrono.h"

// Function to print an error message and terminate the program
void print_error(const char *err) {
    printf("\n\n ERROR: %s \n\n", err);  
    fflush(NULL);  
    exit(1);  
}

// Function to generate a random number between 0 and 1
double random01(unsigned int *seed) {
    return ((double) rand_r(seed) / RAND_MAX);
}

// Function to calculate the Euclidean distance between two nodes
double dist(int i, int j, instance *inst) {
    double dx = inst->xcoord[i] - inst->xcoord[j];  
    double dy = inst->ycoord[i] - inst->ycoord[j];  
    return (double)(sqrt(dx * dx + dy * dy));
}

// Function that calculates a matrix that stores the distance between each pair of nodes
void calculate_distances(instance *inst) {
    inst->distances = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));
    if (inst->distances == NULL){
        print_error("Memory allocation failed");
        return;
    }
    for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            inst->distances[i * inst->nnodes + j] = dist(i, j, inst);
        }
    }
}

// Function to save the solution in a file that will be used by gnuplot
void export_solution_for_gnuplot(const char *filename, const instance *inst) {
    FILE *fp = fopen(filename, "w");  
    if (!fp) {
        print_error("Error opening the output file");
        return;
    }
    for (int i = 0; i < inst->nnodes; i++) {
        int idx = inst->best_sol[i] - 1;  // Adjust for 1-based indexing
        if (idx < 0 || idx >= inst->nnodes) continue;
        fprintf(fp, "%lf %lf %lf\n", inst->xcoord[idx], inst->ycoord[idx], inst->best_sol[i]);
    }
    // If the tour is a cycle, return to the starting node
    if (inst->nnodes > 0) {
        int idx = inst->best_sol[0] - 1;
        fprintf(fp, "%lf %lf %lf\n", inst->xcoord[idx], inst->ycoord[idx], inst->best_sol[0]);
    }

    fclose(fp);  // Close the file
}

// Function to create a PNG image of the solution using gnuplot
void png_solution_for_gnuplot(const char *dat_filename, const char *png_filename) {
    char command[256];
    snprintf(command, sizeof(command),
             "gnuplot -e \"set terminal png size 1280,900; set output '%s'; plot '%s' using 1:2 with lines linecolor 'blue' title 'TSP Path', '' using 1:2:3 with points pointtype 7 pointsize 1.5 linecolor 'red' notitle\"",
             png_filename, dat_filename);
    int ret = system(command);
    if (ret == -1) {
        print_error("Error executing gnuplot command");
    }
}

// Function to calculate instance's best solution cost
double calculate_tour_cost(const double *tour, instance *inst) {
    double cost = 0.0;
    for (int i = 0; i < inst->nnodes; i++) {
        int j = (i + 1) % inst->nnodes;  
        cost += dist(tour[i] - 1, tour[j] - 1, inst); 
    }
    return cost;
}

// Function to check the feasibility of a tour. It checks if nodes compare only once in the tour
bool check_tour_feasability(const double *tour, instance *inst) {
    bool *visited = calloc(inst->nnodes, sizeof(bool));  
    if (!visited){
        print_error("Memory allocation error for visited");
        } 

    for (int i = 0; i < inst->nnodes; i++) {
        if (tour[i] < 1 || tour[i] > inst->nnodes) {
            free(visited);  
            printf("Node %d out of bounds\n", (int)tour[i]);
            return false;
        }
        if (visited[(int)tour[i] - 1]) {
            free(visited);  
            printf("Node %d visited more than once\n", (int)tour[i]);
            return false;
        }
        visited[(int)tour[i] - 1] = true;
    }

    free(visited);  
    return true;
}

// Function to compare input tour to best solution
void check_solution(double* tour, double cur_sol_cost, instance *inst) {
    double cost = inst->best_sol_cost;  
    if (cur_sol_cost < cost) {  
        if (VERBOSE >= 60){
        printf("New best solution found with cost: %.10f\n", cur_sol_cost);  
        }
    update_best_solution(tour, cur_sol_cost, inst);  
    }
}

// Function to update the best solution
void update_best_solution(double* tour, double cur_sol_cost, instance *inst) {
    inst->best_sol_cost = cur_sol_cost;  
    memmove(inst->best_sol, tour, inst->nnodes * sizeof(double)); 
    }

// Function to find the nearest neighbor tour for the TSP
void nearest_neighbor(instance *inst) {
    inst->best_sol = (double *) calloc(inst->nnodes, sizeof(double));
    if (inst->best_sol == NULL){
        print_error("Memory allocation failed");
    }
    if (VERBOSE >= 50){
        printf("----------------------------------------------------------------------------------------------\n\n");
        printf("Nearest Neighbor calculations\n");
    }
    double t1 = second();  // Start time
    double best_nn_tour_cost = 1e20;  // Initialize the best tour cost for nearest neighbor

    for (int start = 0; start < inst->nnodes; start++) {
        double *tour = (double *) calloc(inst->nnodes, sizeof(double));  
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
                    double d = dist(current, j, inst);  
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

        double cur_sol_cost = calculate_tour_cost(tour, inst);  
        if (cur_sol_cost < best_nn_tour_cost) {
            best_nn_tour_cost = cur_sol_cost;  
        }

        bool solution_is_fiseable = check_tour_feasability(tour, inst);  
        if (solution_is_fiseable) {
            // Update time left
            double t2 = second();
            t1 = t2;

            // Check if the time limit has been reached
            if (inst->time_left <= 0){
                if (VERBOSE >= 50) printf("Time limit reached\n");
                if (VERBOSE >= 30){
                    printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
                    printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
                }            
                free(tour);  
                free(visited);  
                return;
            }
            check_solution(tour, cur_sol_cost, inst);
            two_opt(tour, inst);
        }

        free(tour);
        free(visited);  
    }
    if (VERBOSE >= 30){
        printf("NEAREST NEIGHBOR BEST FINAL COST: %lf\n", best_nn_tour_cost);
        printf("UPDATED COST AFTER 2-OPT: %lf\n", inst->best_sol_cost);
    }
}

// Function to implement 2-opt heuristic for TSP that uses instance's best solution
void two_opt(double *solution, instance *inst) {

    double *tour = (double *) calloc(inst->nnodes, sizeof(double));  
    if (tour == NULL){
        print_error("Memory allocation failed");
    }

    // Initialize the tour with the best solution
    memmove(tour, solution, inst->nnodes * sizeof(double));
    double t1 = second();  // Start time
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt INITIAL COST: %lf\n", calculate_tour_cost(tour, inst));
    }

    bool improved = true;  
    while(improved) {
        improved = false;
        for (int i = 0; i < inst->nnodes ; i++) {
            for (int j = i+1; j < inst->nnodes; j++) {
    
                // Calculate the new tour cost
                int next_j = (j + 1 == inst->nnodes) ? 0 : j + 1; // Node j+1 is node 0 if node j is the last node of the solution
                double cost_delta = inst->distances[(int)(tour[i]-1) * inst->nnodes + (int)(tour[j]-1)] + inst->distances[(int)(tour[i+1]-1) * inst->nnodes + (int)(tour[next_j]-1)] - 
                                    inst->distances[(int)(tour[i]-1) * inst->nnodes + (int)(tour[i+1]-1)] - inst->distances[(int)(tour[j]-1) * inst->nnodes + (int)(tour[next_j]-1)];

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
                    if (VERBOSE >= 60){
                        fprintf(stdout, "MODIFIED COST: %lf\n", calculate_tour_cost(tour, inst));
                    }
                    check_solution(tour, calculate_tour_cost(tour, inst), inst);

                    // Update time left
                    double t2 = second();
                    inst->time_left -= (t2 - t1);
                    t1 = t2;

                    // Check if the time limit has been reached
                    if (inst->time_left <= 0) {
                        if (VERBOSE >= 50) printf("Time limit reached in two_opt\n");
                        free(tour);  
                        return;
                    }
                }
            }
        }
    }
    if (VERBOSE >= 60){
        fprintf(stdout, "2-opt FINAL COST: %lf\n", inst->best_sol_cost);
        printf("----------------------------------------------------------------------------------------------\n\n");

    }

    free(tour);  
}