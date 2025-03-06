#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "tsp_functions.h"

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
    if (!inst->integer_costs)
        return sqrt(dx * dx + dy * dy);  // Return the Euclidean distance (non-integer)
    return (int)(sqrt(dx * dx + dy * dy) + 0.499999999); // Rounding to integer distance
}

// Function that create a simple solution for the TSP by visiting nodes in order
void simple_solution(instance *inst) {
    // Allocate memory for the best solution
    inst->best_sol = (double *) calloc(inst->nnodes, sizeof(double));
    if (inst->best_sol == NULL) print_error("Memory allocation failed");

    // Initialize the best solution with the order of nodes
    for (int i = 0; i < inst->nnodes; i++) {
        inst->best_sol[i] = i+1;
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
