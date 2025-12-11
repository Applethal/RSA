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
  size_t num_constraints; // Number of constraints
  size_t num_vars;        // Number of variables
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
  double BIG_M; // Big constant for the artificial variables, will equal max(coeffs) * 2
} Model;

Model *ReadCsv(FILE *textfile); // Reads CSV file and returns a Model
void TransformModel(Model *model); // Transforms the model to the canonical form 
void FreeModel(Model *model); // Free the model pointers
void InvertMatrix(double **matrix, size_t n); // Inverts a matrix, used only to invert the basis starting from the second RSA iteration
void RevisedSimplex(Model *model); // RSA Loop 
double **Get_BasisInverse(Model *model, int iteration); // returns the B matrix
void PrintColumns(Model *model); // Prints the constraints matrix columns 
double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector); // returns the reduces prices 
double *Get_SimplexMultiplier(Model *model, double **B_inv); // returns simplex multipliers for the exiting variables
double *Get_pivot_column(double **B_inv, Model *model, int best_cost_idx); // Returns the pivot column
void UpdateRhs(Model *model, double *rhs_vector_copy, double **B); // Updates the RHS vector at each RSA iteration 
void Get_ObjectiveFunction(Model *model, double *rhs_vector); // returns the Objective function at the end, if the solving succeeds
void RevisedSimplex_Debug(Model *model); // If the "-Debug" is used, the RSA loop will print the algorithm's output at each iteration
void ValidateModelPointers(Model *model); // Makes sure that each model pointer is valid before executing the RSA loop
void PrintHelp(); // If no arguments are provided, a mock-up man page will be printed and the program will exit

#endif
