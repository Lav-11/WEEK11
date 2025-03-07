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
    printf("\n\n ERROR: %s \n\n", err);  // Print the error message
    fflush(NULL);  // Ensure all output is flushed
    exit(1);  // Exit the program with a non-zero status to indicate error
}

// Function to generate a random number between 0 and 1
double random01(unsigned int *seed) {
    return ((double) rand_r(seed) / RAND_MAX);
}

// Function to calculate the Euclidean distance between two nodes
double dist(int i, int j, instance *inst) {
    double dx = inst->xcoord[i] - inst->xcoord[j];  // Calculate difference in x-coordinates
    double dy = inst->ycoord[i] - inst->ycoord[j];  // Calculate difference in y-coordinates
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
    FILE *fp = fopen(filename, "w");  // Open the output file for writing
    if (!fp) {
        perror("Error opening the output file");
        return;
    }
    // Write the node coordinates and their labels (node number)
    for (int i = 0; i < inst->nnodes; i++) {
        int idx = inst->best_sol[i] - 1;  // Adjust for 1-based indexing
        if (idx < 0 || idx >= inst->nnodes) continue;

        // Write the coordinates and node number to the file
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
             "gnuplot -e \"set terminal png size 800,600; set output '%s'; plot '%s' using 1:2 with linespoints title 'TSP Solution'\"",
             png_filename, dat_filename);
    int ret = system(command);
    if (ret == -1) {
        perror("Error executing gnuplot command");
    }
}

// Function to calculate instance's best solution cost
double calculate_tour_cost(const double *tour, instance *inst) {
    double cost = 0.0;
    for (int i = 0; i < inst->nnodes; i++) {
        int j = (i + 1) % inst->nnodes;  // Next node in the tour
        cost += dist(tour[i] - 1, tour[j] - 1, inst);  // Calculate the distance
    }
    return cost;
}

// Function to compare input tour to best solution
void check_solution(double* tour, double cur_sol_cost, instance *inst) {
    double cost = inst->best_sol_cost;  // Get the best solution cost
    if (cur_sol_cost < cost) {  
        if (VERBOSE >= 50){
        printf("New best solution found\n");  // Print a message
        printf("New cost: %.10f\n", cur_sol_cost);  // Print the new cost
        }
    update_best_solution(tour, cur_sol_cost, inst);  // Update the best solution
    }
}

// Function to update the best solution
void update_best_solution(double* tour, double cur_sol_cost, instance *inst) {
    inst->best_sol_cost = cur_sol_cost;  // Update the best solution cost
    memmove(inst->best_sol, tour, inst->nnodes * sizeof(double));  // Copy the tour to the best solution
    }

// Function to find the nearest neighbor tour for the TSP
void nearest_neighbor(instance *inst) {
    inst->best_sol = (double *) calloc(inst->nnodes, sizeof(double));
    if (inst->best_sol == NULL){
        print_error("Memory allocation failed");
        return;
    }
    double t1 = second();  // Start time

    for (int start = 0; start < inst->nnodes; start++) {
        double *tour = (double *) calloc(inst->nnodes, sizeof(double));  // Allocate memory for the tour
        bool *visited = calloc(inst->nnodes, sizeof(bool));  // Array to track visited nodes
        if (!visited){
            print_error("Memory allocation error for visited");
            return;
        } 

        int current = start;  // Start at the current node
        tour[0] = current + 1;  // First node in the tour
        visited[current] = true;
        // For each node, find the nearest unvisited neighbor
        for (int i = 1; i < inst->nnodes; i++) {
            double minDist = 1e20;
            int nextNode = -1;

            // Search for the nearest neighbor
            for (int j = 0; j < inst->nnodes; j++) {
                if (!visited[j]) {
                    double d = dist(current, j, inst);  // Calculate distance
                    if (d < minDist) {
                        minDist = d;
                        nextNode = j;
                    }
                }
            }

            if (nextNode == -1) print_error("Error constructing the tour");

            tour[i] = (double)nextNode + 1;  // Assign the next node to the tour
            visited[nextNode] = true;  // Mark it as visited
            current = nextNode;  // Move to the next node
        }

        double cur_sol_cost = calculate_tour_cost(tour, inst);  // Calculate the tour cost
        check_solution(tour, cur_sol_cost, inst);  // Check if this is the best solution

        // Update time left
        double t2 = second();
        inst->time_left -= (t2 - t1);
        t1 = t2;

        // Check if the time limit has been reached
        if (inst->time_left <= 0){
            if (VERBOSE >= 50) printf("Time limit reached\n");            
            free(tour);  // Free the tour array
            free(visited);  // Free the visited array
            return;
        }

        free(tour);  // Free the tour array
        free(visited);  // Free the visited array
    }
}

// Function to implement 2-opt heuristic for TSP that uses instance's best solution
void two_opt(instance *inst) {
    // Check if the time limit has already been reached
    if (inst->time_left <= 0){
        if (VERBOSE >= 70) printf("Time limit reached before starting 2-opt\n");
        return;
    }

    double *tour = (double *) calloc(inst->nnodes, sizeof(double));  // Allocate memory for the tour
    if (tour == NULL){
        print_error("Memory allocation failed");
        return;
    }
    double t1 = second(); 

    // Initialize the tour with the best solution
    memmove(tour, inst->best_sol, inst->nnodes * sizeof(double));

    if (VERBOSE >= 50){
        printf("----------------------------------------------------------------------------------------------\n\n");
        printf("Starting 2-opt heuristic\n");
        fprintf(stdout, "INITIAL COST: %lf\n", calculate_tour_cost(tour, inst));
    }

    bool improved = true;  // Flag to track if the tour has been improved
    while(improved) {
        improved = false;
        for (int i = 0; i < inst->nnodes ; i++) {
            for (int j = i; j <= fmin(inst->nnodes-2+i, inst->nnodes-1); j++) {
                if (abs(i-j) < 1 || abs(i-j-inst->nnodes)<1) continue;  // Skip adjacent nodes
    
                // Calculate the new tour cost
                int next_j = (j + 1 == inst->nnodes) ? 0 : j + 1;
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
                    improved = true;  // Mark the tour as improved
                    if (VERBOSE >= 50){
                        fprintf(stdout, "MODIFIED COST: %lf\n", calculate_tour_cost(tour, inst));
                    }
                    update_best_solution(tour, calculate_tour_cost(tour, inst), inst);  // Update the best solution

                    // Update time left
                    double t2 = second();
                    inst->time_left -= (t2 - t1);
                    t1 = t2;

                    // Check if the time limit has been reached
                    if (inst->time_left <= 0) {
                        if (VERBOSE >= 50) printf("Time limit reached in two_opt\n");
                        free(tour);  // Free the tour array
                        return;
                    }
                }
            }
        }
    }

    free(tour);  // Free the tour array
}