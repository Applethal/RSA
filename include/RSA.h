#ifndef RSA_H
#define RSA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> // Do I really need this to setup the Big Ms? lol

typedef struct
{
  char *objective;     // MINIMIZE or MAXIMIZE
  int num_constraints; // Number of constraints
  int num_vars;        // Number of variables
  double **columns;    // Constraints x variables coeffs matrix
  double *coeffs;      // Variable objective coefficients
  double *rhs_vector;  // Right hand side vector of the constraints
  int *basics_vector;  // Vector where I keep track of the basic variables
  double objective_function;
  char *constraints_symbols; // Tracks the constraints' symbols, for debugging purposes only
  int inequalities_count;    // Counts the number of inequalities
  int equalities_count;      // Counts the number of equality constraints
  int *equalities_vector;    // Contains the indices of the equality constraints
  int solver_iterations;
  int *non_basics;      // Contains indices of non-basic variables
  int non_basics_count; // Counts the number of non-basic variables
} Model;

Model *ReadCsv(FILE *textfile);
void FreeModel(Model *model);
void InvertMatrix(double **matrix, int n);
void RevisedSimplex(Model *model);
double **Get_BasicsMatrix(Model *model);
void PrintColumns(Model *model);
double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector);
double *Get_SimplexMultiplier(Model *model, double **B_inv);
double *Get_pivot_column(double **B_inv, Model *model, int best_cost_idx);
void UpdateRhs(Model *model, double *rhs_vector_copy, double **B);
void Get_ObjectiveFunction(Model *model, double *rhs_vector);
void RevisedSimplex_Debug(Model *model);
void ValidateModelPointers(Model *model);
void PrintHelp();
#endif
