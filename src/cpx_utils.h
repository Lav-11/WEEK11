#ifndef CPX_UTILS_H
#define CPX_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <ilcplex/cplex.h>
#include "tsp_utils.h"
#include "solver.h"


// Function prototypes

// Compute position index for x(i, j)
int xpos(int i, int j, instance *inst);

// Build the CPLEX model
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);

void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);

// Solve the TSP using CPLEX
int TSPopt(instance *inst, ConfigParams *params);

// Function to ensure the directory exists
void create_directory(const char *path);

// Function to plot the graph and save it as an image
void plot_graph_to_image(int nnodes, double *xcoord, double *ycoord, double *xstar, instance *inst, double max_coord, double padding);

// Function to use benders
void benders_loop(instance *inst, CPXENVptr env, CPXLPptr lp, double start_time);

// Function to build the solution from xstar
void add_SEC_constraints(instance *inst, CPXENVptr env, CPXLPptr lp, double *xstar,
CPXCALLBACKCONTEXTptr context, int contextid);

void patch_solution(double *xstar, instance *inst);

void invert_path(int start, int end, int *succ, double *xstar, instance *inst);

// warm start function
void warmstart(CPXENVptr env, CPXLPptr lp, instance *inst);

// convert successor array to tour
void convert_succ_to_tour(int *succ, int nnodes, double *tour);


#endif // CPX_UTILS_H