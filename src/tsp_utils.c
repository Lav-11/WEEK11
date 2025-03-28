#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "tsp_utils.h"
#include "heuristics.h"
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
            inst->distances[i * inst->nnodes + j] = round(dist(i, j, inst));
        }
    }
}

// Function to save the solution in a png using gnuplot
void png_solution_for_gnuplot(solution *sol, const bool save_dat_file, char *output_filename, instance *inst) {
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (!gnuplotPipe) {
        print_error("Error opening pipe to gnuplot");
        return;
    }

    fprintf(gnuplotPipe, "set terminal png size 1280,900\n");
    fprintf(gnuplotPipe, "set output '%s.png'\n", output_filename);
    fprintf(gnuplotPipe, "set title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "plot '-' using 1:2 with lines linewidth 2 linecolor 'blue' title 'TSP Path', '-' using 1:2 with points pointtype 7 pointsize 2 linecolor 'red' notitle\n");

    for (int i = 0; i < inst->nnodes+1; i++) {
        int idx = (int)sol->tour[i] - 1;
        fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
    }

    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < inst->nnodes+1; i++) {
        int idx = (int)sol->tour[i] - 1;
        fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
    }

    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);
    pclose(gnuplotPipe);

    if (save_dat_file) {
        char dat_filename[256];
        snprintf(dat_filename, sizeof(dat_filename), "%s.dat", output_filename);
        FILE *datFile = fopen(dat_filename, "w");
        if (!datFile) {
            print_error("Error opening solution .dat file");
            return;
        }
        for (int i = 0; i < inst->nnodes+1; i++) {
            int idx = (int)sol->tour[i] - 1;
            fprintf(datFile, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
        }
        fclose(datFile);
    }
}


// Function to calculate tour cost
void calculate_tour_cost(solution *sol, const instance *inst) {
    double total_cost = 0.0;
    for (int i = 0; i < inst->nnodes; i++) {
        int from = (int)sol->tour[i] - 1;
        int to = (int)sol->tour[i + 1] - 1;
        total_cost += inst->distances[from * inst->nnodes + to];
    }
    sol->tour_cost = total_cost;
}


// Function to check the feasibility of a tour. It checks if nodes compare only once in the tour
bool check_tour_feasability(solution *sol, instance *inst) {
    bool *visited = calloc(inst->nnodes, sizeof(bool));  
    if (!visited){
        print_error("Memory allocation error for visited");
        } 

    for (int i = 0; i < inst->nnodes; i++) {
        if (sol->tour[i] < 1 || sol->tour[i] > inst->nnodes) {
            free(visited);  
            printf("Node %d out of bounds\n", (int)sol->tour[i]);
            return false;
        }
        if (visited[(int)sol->tour[i] - 1]) {
            free(visited);  
            printf("Node %d visited more than once\n", (int)sol->tour[i]);
            return false;
        }
        visited[(int)sol->tour[i] - 1] = true;
    }

    for (int i = 0; i < inst->nnodes; i++) {
        if (!visited[i]) {
            free(visited);  
            printf("Node %d not visited\n", (int)sol->tour[i]);
            return false;
        }
    }

    free(visited);  
    return true;
}

// Function to compare input tour to best solution
void check_if_best_solution(solution *sol, instance *inst) {
    if (sol->tour_cost < inst->best_sol->tour_cost) {  
        if (check_tour_feasability(sol, inst)){      //Divide the 2 if for better efficiency  
            if (VERBOSE >= 60){
            printf("New best solution found with cost: %.10f\n", sol->tour_cost);  
            }
            update_best_solution(sol, inst);  
        }
    }
}

// Function to update the best solution
void update_best_solution(solution *sol, instance *inst) {
    memcpy(inst->best_sol->tour, sol->tour, (inst->nnodes+1) * sizeof(double)); 
    inst->best_sol->tour_cost = sol->tour_cost; 
    }

// Function to free memory of a solution struct
void free_solution(solution *sol) {
    free(sol->tour);
    free(sol);
}

// Function that read a file containing costs and create a png file containing a graph showing them of it using gnuplot
void plot_costs(char *input_filename, char *output_filename) {
    FILE *cost_file = fopen(input_filename, "r");
    if (!cost_file) {
        print_error("Error opening cost file");
        return;
    }

    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (!gnuplotPipe) {
        print_error("Error opening pipe to gnuplot");
        fclose(cost_file);
        return;
    }

    fprintf(gnuplotPipe, "set terminal png size 1280,900\n");
    fprintf(gnuplotPipe, "set output '%s.png'\n", output_filename);
    fprintf(gnuplotPipe, "set title 'Costs'\n");
    fprintf(gnuplotPipe, "set xlabel 'Iteration'\n");
    fprintf(gnuplotPipe, "set ylabel 'Cost'\n");
    fprintf(gnuplotPipe, "plot '-' with lines linewidth 2 linecolor 'blue' title 'Cost', '-' with lines linewidth 2 linecolor 'red' title 'Minimum Cost'\n");

    char line[1024];
    int count = 0;
    double min_cost = INFINITY;

    // Plot the costs and track the minimum cost dynamically
    while (fgets(line, sizeof(line), cost_file)) {
        line[strcspn(line, "\r\n")] = 0;
        double current_cost = atof(line);
        if (current_cost < min_cost) {
            min_cost = current_cost;
        }
        fprintf(gnuplotPipe, "%d %lf\n", count, current_cost);
        count++;
    }
    fprintf(gnuplotPipe, "e\n");

    // Reset file pointer to re-read the costs for plotting the dynamic minimum cost
    rewind(cost_file);
    count = 0;
    min_cost = INFINITY;

    while (fgets(line, sizeof(line), cost_file)) {
        line[strcspn(line, "\r\n")] = 0;
        double current_cost = atof(line);
        if (current_cost < min_cost) {
            min_cost = current_cost;
        }
        fprintf(gnuplotPipe, "%d %lf\n", count, min_cost);
        count++;
    }
    fprintf(gnuplotPipe, "e\n");

    fflush(gnuplotPipe);
    pclose(gnuplotPipe);

    if (fclose(cost_file) != 0) {
        print_error("Error closing cost file");
    }
}

