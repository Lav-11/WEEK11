#ifndef CALLBACK_H
#define CALLBACK_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <ilcplex/cplex.h>
#include "tsp_utils.h"
#include "heuristics.h"



int CPXPUBLIC candidate_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);

#endif