#ifndef RSA_H
#define RSA_H

#include <stddef.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define ARTIFICIAL 2
#define STANDARD 1
#define SLACK 3

typedef struct {
  double value;
  int type; // 1, 2 and 3 for standard, artificial and slack variables respectively
} Variable;

typedef struct
{
  char *objective;           // MINIMIZE or MAXIMIZE
  size_t num_constraints;    // Number of constraints
  size_t num_vars;           // Number of variables
  double **lhs_matrix;       // Constraints Left hand side
  Variable *coeffs;          // Variable objective coefficients
  double *rhs_vector;        // Right hand side vector of the constraints
  int *basics_vector;        // Vector where I keep track of the basic variables
  double objective_function; // Self explanatory, stores the model's objective function
  char *constraints_symbols; // Tracks the constraints' symbols, for debugging purposes only
  int slacks_surplus_count;  // Counts the number of slack and surplus vars
  int artificials_count;     // Counts the number of artificial vars
  int *artificials_vector;   // Contains the indices of artificial vars
  int solver_iterations;     // Self explanatory, tracks the solver iterations count 
  int *non_basics;           // Contains indices of non-basic variables
  int non_basics_count;      // Counts the number of non-basic variables
} Model;

// Function declarations
void PrintHelp();
Model *ReadCsv(FILE *csvfile);
void TransformModel(Model *model);
void Printlhs_matrix(Model *model);
double **Get_BasisInverse(Model *model, int iteration);
void InvertMatrix(double **matrix, size_t n);
void RevisedSimplex(Model *model);
void RevisedSimplex_Debug(Model *model);
double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector);
double *Get_SimplexMultiplier(Model *model, double **B_inv);
double *Get_pivot_column(double **B_inv, Model *model, int best_cost_idx);
void UpdateRhs(Model *model, double *rhs_vector_copy, double **B);
void Get_ObjectiveFunction(Model *model, double *rhs_vector);
void FreeModel(Model *model);
void ValidateModelPointers(Model *model);

#endif // RSA_H
