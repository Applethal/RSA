#include "RSA.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> // Do I really need this to setup the Big Ms? lol

#define MAX_LINE_LENGTH 1024 
#define MAX_VARS 2000
#define MAX_CONSTRAINTS 2000 
#define MAX_MEM_BYTES 1024 // this will limit the constraints matrix memory size, this should be enough for 2000 vars and 2000 constraints



void PrintHelp(){

  printf("NAME\n");
  printf("     RSA - Linear programming solver using the Revised Simplex Algorithm with the Big M method\n\n");

  printf("SYNOPSIS\n");
  printf("     ./RSA csv file path [-Debug]\n\n");

  printf("DESCRIPTION\n");
  printf("     RSA is a linear programming solver that implements the Revised\n");
  printf("     Simplex Algorithm. It reads a linear programming model from a CSV file\n");
  printf("     and computes the optimal solution for the given objective function and\n");
  printf("     constraints.\n\n");
  printf("     The program supports both maximization and minimization objectives and\n");
  printf("     automatically handles equality and inequality constraints by adding slack\n");
  printf("     and artificial variables as needed.\n\n");
  printf("     Visit the Github repo https://github.com/Applethal/RSA to see an example of how the CSV data input should look like.\n");

  printf("OPTIONS\n");
  printf("     csv file path\n");
  printf("             Path to the input CSV file containing the linear programming\n");
  printf("             model. This argument is required. The CSV file should be properly\n");
  printf("             formatted with the objective function, variables, and constraints.\n\n");
  printf("     -Debug  Enables debug mode. When this flag is provided, the program\n");
  printf("             displays detailed information about the model and shows iterative\n");
  printf("             steps during the solving process. \n");

  printf("EXIT STATUS\n");
  printf("     0       Optimal solution obtained.\n");
  printf("     1       File opening error (file does not exist or cannot be accessed).\n");
  printf("     Other   Model is infeasible or unbounded.\n\n");


  printf("OUTPUT\n");
  printf("     The program displays:\n");
  printf("     - Start message with usage hint\n");
  printf("     - Debug information (if -Debug flag is used)\n");
  printf("     - RSA iterations count\n");
  printf("     - Final objective function value and the values for each variables in the final basis.\n");
  printf("     - Solving time in seconds\n");
  printf("     - Model's size in bytes\n");

  printf("AUTHOR\n");
  printf("     Written by Applethal / Saad Nordine\n\n");

  printf("REPORTING BUGS\n");
  printf("     Please report any bugs in the issues page in https://github.com/Applethal/RSA\n\n");





}

Model *ReadCsv(FILE *csvfile)
{
  Model *model = (Model *)malloc(sizeof(Model));
  char line[MAX_LINE_LENGTH];

  // Line 1: Read objective type, num_vars, num_constraints
  if (!fgets(line, sizeof(line), csvfile))
  {
    fprintf(stderr, "Error: Could not read header line\n");
    free(model);
    return NULL;
  }

  size_t len = strlen(line);
  if (len == sizeof(line) - 1 && line[len-1] != '\n') {
    fprintf(stderr, "Error: Line too long (max %d chars)\n", MAX_LINE_LENGTH);
    free(model);
    return NULL;
  }

  // Remove newline
  line[strcspn(line, "\r\n")] = 0;

  // Parse: OBJECTIVE,num_vars,num_constraints
  char *token = strtok(line, ",");
  if (!token) {
    fprintf(stderr, "Error: Invalid header format\n");
    free(model);
    return NULL;
  }

  model->objective = (char *)malloc(strlen(token) + 1);
  strcpy(model->objective, token);

  token = strtok(NULL, ",");
  if (!token) {
    fprintf(stderr, "Error: Missing variables count in header\n");
    free(model->objective);
    free(model);
    return NULL;
  }
  model->num_vars = atoi(token);

  token = strtok(NULL, ",");
  if (!token) {
    fprintf(stderr, "Error: Missing constraints count in header\n");
    free(model->objective);
    free(model);
    return NULL;
  }
  model->num_constraints = atoi(token);

  // Validate limits
  if (model->num_vars > MAX_VARS) {
    fprintf(stderr, "Error: Too many variables (%zu). Maximum allowed: %d\n", 
            model->num_vars, MAX_VARS);
    free(model->objective);
    free(model);
    return NULL;
  }

  if (model->num_constraints > MAX_CONSTRAINTS) {
    fprintf(stderr, "Error: Too many constraints (%zu). Maximum allowed: %d\n",
            model->num_constraints, MAX_CONSTRAINTS);
    free(model->objective);
    free(model);
    return NULL;
  }

  // Line 2: Read objective coefficients
  if (!fgets(line, sizeof(line), csvfile))
  {
    fprintf(stderr, "Error: Could not read objective coefficients\n");
    free(model->objective);
    free(model);
    return NULL;
  }

  // Allocate and parse coefficients
  model->coeffs = (double *)malloc(model->num_vars * sizeof(double));
  token = strtok(line, ",");
  int idx = 0;
  while (token != NULL && idx < model->num_vars)
  {
    model->coeffs[idx++] = atof(token);
    token = strtok(NULL, ",");
  }

  if (idx != model->num_vars) {
    fprintf(stderr, "Error: Expected %zu variable coefficients but got %d\n", 
            model->num_vars, idx);
    free(model->coeffs);
    free(model->objective);
    free(model);
    return NULL;
  }

  // Allocate arrays for raw constraint data
  model->lhs_matrix = (double **)malloc(model->num_constraints * sizeof(double *));
  for (int i = 0; i < model->num_constraints; i++)
  {
    model->lhs_matrix[i] = (double *)malloc(model->num_vars * sizeof(double));
  }

  model->rhs_vector = (double *)malloc(model->num_constraints * sizeof(double));
  model->constraints_symbols = (char *)malloc(model->num_constraints * sizeof(char));
  int num_le = 0, num_ge = 0, num_eq = 0;

  // Read all constraints in one pass
  for (int i = 0; i < model->num_constraints; i++)
  {
    if (!fgets(line, sizeof(line), csvfile))
    {
      fprintf(stderr, "Error: Unexpected end of file at constraint %d. This might indicate a missing constraint or wrong constraints count in the header\n", i);
      // Cleanup on error
      for (int j = 0; j <= i; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model->objective);
      free(model);
      return NULL;
    }

    // Parse constraint coefficients
    token = strtok(line, ",");
    int col_idx = 0;

    // Read variable coefficients
    while (token != NULL && col_idx < model->num_vars)
    {
      model->lhs_matrix[i][col_idx] = atof(token);
      col_idx++;
      token = strtok(NULL, ",");
    }

    if (col_idx != model->num_vars) {
      fprintf(stderr, "Error: Constraint %d has %d coefficients, expected %zu variables\n", 
              i, col_idx, model->num_vars);
      // Cleanup
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model->objective);
      free(model);
      return NULL;
    }

    // Read operator (<=, >=, =)
    if (token == NULL)
    {
      fprintf(stderr, "Error: Missing operator in constraint %d\n", i);
      // Cleanup
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model->objective);
      free(model);
      return NULL;
    }

    // Remove spaces from operator token
    char op[4] = "";
    int op_idx = 0;
    for (char *p = token; *p && op_idx < 3; p++)
    {
      if (*p != ' ' && *p != '\t')
      {
        op[op_idx++] = *p;
      }
    }
    op[op_idx] = '\0';

        // Store constraint type
    if (strcmp(op, "<=") == 0)
    {
      model->constraints_symbols[i] = 'L';
      num_le++;

    }
    else if (strcmp(op, ">=") == 0)
    {
      model->constraints_symbols[i] = 'G';
      num_ge++;

    }
    else if (strcmp(op, "=") == 0)
    {
      model->constraints_symbols[i] = 'E';
      num_eq++;

    }
    

    else
    {
      fprintf(stderr, "Error: Invalid operator '%s' in constraint %d\n", op, i);
      // Cleanup
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model->objective);
      free(model);
      return NULL;
    }
     
    // Read RHS
    token = strtok(NULL, ",");
    if (token == NULL)
    {
      fprintf(stderr, "Error: Missing RHS in constraint %d\n", i);
      // Cleanup
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model->objective);
      free(model);
      return NULL;
    }
    model->rhs_vector[i] = atof(token);
  }

  model->artificials_vector = NULL;
  model->basics_vector = NULL;
  model->non_basics = NULL;
  model->non_basics_count = 0;
  model->solver_iterations = 1;
  model->objective_function = 0.0;
  model->BIG_M = 0.0;
  model->slacks_surplus_count = num_le + num_ge;  // slack + surplus
  model->artificials_count = num_eq + num_ge;    // artificial variables

  if (model->num_constraints > model->num_vars)
  {
    printf("Warning: More constraints than variables! The revised simplex algorithm works best when there are more variables than constraints!\n");
  }

  return model;
}

void TransformModel(Model *model)
{
  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
  }
  int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;

  // Check memory requirements
  size_t memory_needed = (size_t)model->num_constraints * total_cols * sizeof(double);
  size_t max_memory = (size_t)MAX_MEM_BYTES * 1024 * 1024;

  if (memory_needed > max_memory) {
    fprintf(stderr, "Error: Model too large. Requires %zu MB, maximum is %d MB\n",
            memory_needed / (1024 * 1024), MAX_MEM_BYTES);
    exit(1);
  }

  // Store old lhs_matrix and coeffs
  double **old_lhs_matrix = model->lhs_matrix;
  double *old_coeffs = model->coeffs;
  size_t old_num_vars = model->num_vars;

  // Allocate new expanded lhs_matrix matrix
  model->lhs_matrix = (double **)malloc(model->num_constraints * sizeof(double *));
  for (int i = 0; i < model->num_constraints; i++)
  {
    model->lhs_matrix[i] = (double *)calloc(total_cols, sizeof(double));

    // Copy original constraint coefficients
    for (int j = 0; j < old_num_vars; j++)
    {
      model->lhs_matrix[i][j] = old_lhs_matrix[i][j];
    }
  }

  // Allocate new expanded coefficients array
  model->coeffs = (double *)calloc(total_cols, sizeof(double));

  // Copy original objective coefficients
  for (int i = 0; i < old_num_vars; i++)
  {
    model->coeffs[i] = old_coeffs[i];
  }

  // Free old arrays
  for (int i = 0; i < model->num_constraints; i++)
  {
    free(old_lhs_matrix[i]);
  }
  free(old_lhs_matrix);
  free(old_coeffs);

  // Allocate tracking arrays
  model->artificials_vector = (int *)malloc(model->artificials_count * sizeof(int));
  model->basics_vector = (int *)calloc(model->num_constraints, sizeof(int));

  // Add slack, surplus, and artificial variables
  int slack_col = model->num_vars;
  int artificial_col = model->num_vars + model->slacks_surplus_count;
  int slack_idx = 0;
  int artificial_idx = 0;

  for (int i = 0; i < model->num_constraints; i++)
  {
    if (model->constraints_symbols[i] == 'L')
    {
      // Less-than: add slack variable
      int col = slack_col + slack_idx;
      model->lhs_matrix[i][col] = 1.0;
      model->basics_vector[i] = col;
      model->coeffs[col] = 0.0;
      slack_idx++;
    }
    else if (model->constraints_symbols[i] == 'G')
    {
      // Greater-than: add surplus and artificial
      int surplus_col = slack_col + slack_idx;
      int artif_col = artificial_col + artificial_idx;

      model->lhs_matrix[i][surplus_col] = -1.0;
      model->lhs_matrix[i][artif_col] = 1.0;
      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      slack_idx++;
      artificial_idx++;
    }
    else if (model->constraints_symbols[i] == 'E')
    {
      // Equality: add artificial variable
      int artif_col = artificial_col + artificial_idx;
      model->lhs_matrix[i][artif_col] = 1.0;
      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      artificial_idx++;
    }
  }

  // Calculate Big M based on largest objective coefficient
  double biggest_coeff = 0;
  for (int i = 0; i < model->num_vars; i++)
  {
    if (fabs(model->coeffs[i]) > biggest_coeff) {
      biggest_coeff = fabs(model->coeffs[i]);
    }
  }
  model->BIG_M = biggest_coeff * 2;

  // Apply Big M penalties to artificial variables
  int artificial_start = model->num_vars + model->slacks_surplus_count;
  for (int i = 0; i < model->artificials_count; i++)
  {
    if (strcmp(model->objective, "MAXIMIZE") == 0)
    {
      model->coeffs[artificial_start + i] = -model->BIG_M;
    }
    else
  {
      model->coeffs[artificial_start + i] = model->BIG_M;
    }
  }

  // Convert MINIMIZE to MAXIMIZE by negating all coefficients
  if (strcmp(model->objective, "MINIMIZE") == 0)
  {
    for (int i = 0; i < total_cols; i++)
    {
      model->coeffs[i] *= -1;
    }
  }

  // Build non-basics list
  model->non_basics = (int *)malloc(total_cols * sizeof(int));
  int non_basic_idx = 0;
  for (int col = 0; col < total_cols; col++)
  {
    int is_basic = 0;
    for (int j = 0; j < model->num_constraints; j++)
    {
      if (model->basics_vector[j] == col)
      {
        is_basic = 1;
        break;
      }
    }
    if (!is_basic)
    {
      model->non_basics[non_basic_idx++] = col;
    }
  }
  model->non_basics_count = non_basic_idx;
}



void Printlhs_matrix(Model *model)
{
  int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;

  printf("Constraints in canonical form:\n");
  for (int i = 0; i < model->num_constraints; i++)
  {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < total_cols; j++)
    {
      printf("%5.1f ", model->lhs_matrix[i][j]);
    }
    printf("| RHS: %5.1f", model->rhs_vector[i]);
    printf("\n");
  }
  printf("\n");

  printf("Constraint symbols in the original problem:\n");

  for (size_t i = 0; i < model->num_constraints; i++)
  {
    printf(" %c ", model->constraints_symbols[i]);
  }
  printf("\n");
}

double **Get_BasisInverse(Model *model, int iteration) {
  size_t n = model->num_constraints;
  double **B_inv = (double **)malloc(n * sizeof(double *));

  for (size_t i = 0; i < n; i++) {
    B_inv[i] = (double *)malloc(n * sizeof(double));
  }

  if (iteration == 1) {

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        B_inv[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }
  } else {
    // Build B then invert it
    double **B = (double **)malloc(n * sizeof(double *));
    for (size_t i = 0; i < n; i++) {
      B[i] = (double *)malloc(n * sizeof(double));
      for (size_t j = 0; j < n; j++) {
        int basis_col = model->basics_vector[j];
        B[i][j] = model->lhs_matrix[i][basis_col];
      }
    }

    // Copy B to B_inv for in-place inversion
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        B_inv[i][j] = B[i][j];
      }
    }

    // Free B immediately
    for (size_t i = 0; i < n; i++) {
      free(B[i]);
    }
    free(B);

    InvertMatrix(B_inv, n);
  }

  return B_inv;
}

void InvertMatrix(double **matrix, size_t n)
{
  const double EPS = 1e-12;

  // Single block allocation for the augmented matrix
  double *block = malloc(n * (2 * n) * sizeof(double));
  double **aug = malloc(n * sizeof(double *));
  if (!block || !aug) {
    printf("ERROR: Out of memory\n");
    exit(1);
  }

  for (size_t i = 0; i < n; i++)
    aug[i] = block + i * (2 * n);

  // Build augmented matrix
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++)
      aug[i][j] = matrix[i][j];
    for (size_t j = n; j < 2 * n; j++)
      aug[i][j] = (j - n == i);
  }

  // Gauss-Jordan elimination with partial pivoting
  for (size_t i = 0; i < n; i++) {

    // Find pivot
    size_t max_row = i;
    double max_val = fabs(aug[i][i]);

    for (size_t k = i + 1; k < n; k++) {
      double v = fabs(aug[k][i]);
      if (v > max_val) {
        max_val = v;
        max_row = k;
      }
    }

    // Swap
    if (max_row != i) {
      double *tmp = aug[i];
      aug[i] = aug[max_row];
      aug[max_row] = tmp;
    }

    // Check singularity
    double pivot = aug[i][i];
    if (fabs(pivot) < EPS) {
      printf("ERROR: Singular matrix - basis is not invertible!\n");
      free(aug);
      free(block);
      exit(1);
    }

    // Normalize pivot row
    double inv_pivot = 1.0 / pivot;
    for (size_t j = 0; j < 2 * n; j++)
      aug[i][j] *= inv_pivot;

    // Eliminate
    for (size_t k = 0; k < n; k++) {
      if (k == i) continue;
      double factor = aug[k][i];
      if (factor == 0.0) continue;

      for (size_t j = 0; j < 2 * n; j++)
        aug[k][j] -= factor * aug[i][j];
    }
  }

  // Copy back result
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      matrix[i][j] = aug[i][n + j];

  free(aug);
  free(block);
}

void RevisedSimplex(Model *model) {

  TransformModel(model);

  ValidateModelPointers(model); 


  int termination = 0;
  size_t n = model->num_constraints;
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1;

  while (termination != 1) {
    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      FreeModel(model);
      exit(0);
    }

    int feasibility_check = 0;
    printf("Beginning solver iteration %i ... \n", model->solver_iterations);


    double **B_inv = Get_BasisInverse(model, model->solver_iterations);

    double original_RHS[n];
    for (size_t i = 0; i < n; i++) {
      original_RHS[i] = model->rhs_vector[i];
    }

    double *Simplex_multiplier = Get_SimplexMultiplier(model, B_inv);

    int entering_var_idx = 0;
    int entering_var = 0;
    UpdateRhs(model, original_RHS, B_inv);

    double best_reduced_cost = -DBL_MAX;

    // Getting the entering variable
    for (size_t i = 0; i < model->non_basics_count; i++) {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);

      if (reduced_cost > best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0) {
        feasibility_check++;
      }
    }

    // If all reduced costs are <= 0, optimal solution found
    if (feasibility_check == model->non_basics_count) {
      printf("Solver loop terminated!\n");
      termination++;
      Get_ObjectiveFunction(model, original_RHS);

      for (size_t i = 0; i < n; i++) {
        free(B_inv[i]);
      }
      free(B_inv);
      free(Simplex_multiplier);
      break;
    }

    double *Pivot = Get_pivot_column(B_inv, model, entering_var);

    // Getting the exiting variable
    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = DBL_MAX;

    for (size_t i = 0; i < n; i++) {
      if (Pivot[i] <= 1e-6) {
        continue;
      }
      double ratio = original_RHS[i] / Pivot[i];
      if (ratio < best_ratio && ratio >= 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }

    if (exiting_var_idx == -1) {
      printf("LP is unbounded! Terminating!\n");
      exit(0);
    }

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;
    model->solver_iterations++;

    free(Pivot);
    for (size_t i = 0; i < n; i++) {
      free(B_inv[i]);
    }
    free(B_inv);
    free(Simplex_multiplier);
  }
}


void RevisedSimplex_Debug(Model *model) {


  TransformModel(model);

  ValidateModelPointers(model); 



  printf("Debug Mode: On\n");
  printf("Model after transformation:\n");
  printf("Objective function mode: %s\n", model->objective);
  printf("Number of variables: %d\n", model->num_vars);
  printf("Number of constraints: %d\n", model->num_constraints);
  printf("Objective coefficients:");
  printf("Big M: %.1f", model->BIG_M);
  for (int i = 0; i < model->num_vars; i++) {
    printf(" %.1f", model->coeffs[i]);  
  }
  printf("\n\n");
  Printlhs_matrix(model);



  int termination = 0;
  size_t n = model->num_constraints;
  printf("\n");
  printf("===================================================================\n");

  printf("Initial Non basic variable indices:\n");
  for (size_t i = 0; i < model->non_basics_count; i++) {
    printf(" %i ", model->non_basics[i]);
  }
  printf("\n");
  printf("Initial Basic variable indices:\n");

  for (size_t i = 0; i < model->num_constraints; i++) {
    printf(" %i ", model->basics_vector[i]);
  }
  printf("\n");
  printf("===================================================================\n");
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1;

  while (termination != 1) {
    int feasibility_check = 0;
    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      FreeModel(model);
      exit(0);
    }

    printf("Beginning solver iteration %i ... \n", model->solver_iterations);

    if (model->solver_iterations == 1) {
      printf("Solver iteration 1, B^-1 is identity matrix\n");
    } else {
      printf("Getting B^-1 matrix:\n");
    }

    double **B_inv = Get_BasisInverse(model, model->solver_iterations);

    printf("===================================================================\n");
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        printf(" %f ", B_inv[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    printf("===================================================================\n");

    double original_RHS[n];
    for (size_t i = 0; i < n; i++) {
      original_RHS[i] = model->rhs_vector[i];
    }

    printf("Getting the simplex multiplier P vector:\n");
    double *Simplex_multiplier = Get_SimplexMultiplier(model, B_inv);
    printf("P vector elements:\n");

    for (size_t i = 0; i < n; i++) {
      printf(" %f ", Simplex_multiplier[i]);
    }
    printf("\n");
    printf("===================================================================\n");

    int entering_var_idx = 0;
    int entering_var = 0;
    UpdateRhs(model, original_RHS, B_inv);

    printf("Updating the Right hand side vector\n");
    for (size_t i = 0; i < n; i++) {
      printf(" %f ", model->rhs_vector[i]);
    }
    printf("\n");

    double best_reduced_cost = -DBL_MAX;
    printf("MAXIMIZATION problem (after conversion), choosing the highest reduced cost value\n");
    printf("===================================================================\n");
    printf("Getting the entering variable:\n");
    printf("===================================================================\n");

    for (size_t i = 0; i < model->non_basics_count; i++) {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
      printf("The reduced cost of non basic variable %i is %f \n", non_basic_idx, reduced_cost);

      if (reduced_cost > best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0) {
        feasibility_check++;
      }
    }

    printf("Entering variable index: %i, cost: %f \n", entering_var, best_reduced_cost);
    printf("===================================================================\n");

    if (feasibility_check == model->non_basics_count) {
      printf("Solver loop terminated!\n");
      termination++;
      Get_ObjectiveFunction(model, original_RHS);

      for (size_t i = 0; i < n; i++) {
        free(B_inv[i]);
      }
      free(B_inv);
      free(Simplex_multiplier);
      break;
    }

    double *Pivot = Get_pivot_column(B_inv, model, entering_var);
    printf("Getting the exiting variable:\n");
    printf("===================================================================\n");

    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = DBL_MAX;

    for (size_t i = 0; i < n; i++) {
      if (Pivot[i] <= 1e-6) {
        continue;
      }
      double ratio = original_RHS[i] / Pivot[i];
      printf("Ratio of variable %i is %f which has a pivot value of %f \n", 
             model->basics_vector[i], ratio, Pivot[i]);
      if (ratio < best_ratio && ratio >= 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }

    if (exiting_var_idx == -1) {
      printf("No exiting variable found, LP is unbounded! Terminating!\n");
      exit(0);
    }

    printf("The exiting variable is %i with ratio of %f \n", exiting_var, best_ratio);
    printf("===================================================================\n");

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;

    printf("Entering variable:%i at position %i \n", entering_var, entering_var_idx);
    printf("Exiting variable:%i at position: %i \n", exiting_var, exiting_var_idx);
    printf("Non basics updated:\n");
    for (size_t i = 0; i < model->non_basics_count; i++) {
      printf(" %i ", model->non_basics[i]);
    }
    printf("\n");
    printf("Basics vector:\n");
    for (size_t i = 0; i < model->num_constraints; i++) {
      printf(" %i ", model->basics_vector[i]);
    }

    model->solver_iterations++;
    free(Pivot);
    for (size_t i = 0; i < n; i++) {
      free(B_inv[i]);
    }
    free(B_inv);
    free(Simplex_multiplier);

    printf("\n");
    printf("Iteration %i ends here, press any key to continue!\n", 
           model->solver_iterations - 1);
    printf("===================================================================\n");
    getchar();
  }
}

double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector)
{
  size_t n = model->num_constraints;
  double dot_product = 0.0;

  for (int i = 0; i < n; i++)
  {
    dot_product += multiplier_vector[i] * model->lhs_matrix[i][var_col];
  }

  double reduced_cost = model->coeffs[var_col] - dot_product;
  return reduced_cost;
}

double *Get_SimplexMultiplier(Model *model, double **B_inv)
{
  size_t n = model->num_constraints;
  double *multiplier_vector = (double *)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++)
  {
    double sum = 0.0;
    for (int j = 0; j < n; j++)
    {
      int basic_col_idx = model->basics_vector[j];
      sum += model->coeffs[basic_col_idx] * B_inv[j][i];
    }
    multiplier_vector[i] = sum;
  }

  return multiplier_vector;
}

double *Get_pivot_column(double **B_inv, Model *model, int best_cost_idx)
{
  size_t n = model->num_constraints;
  double *Pivot = (double *)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++)
  {
    double sum = 0.0;
    for (int j = 0; j < n; j++)
    {
      sum += B_inv[i][j] * model->lhs_matrix[j][best_cost_idx];
    }
    Pivot[i] = sum;
  }

  return Pivot;
}

void UpdateRhs(Model *model, double *rhs_vector_copy, double **B)
{
  size_t n = model->num_constraints;
  double *temp = (double *)malloc(n * sizeof(double));

  for (int i = 0; i < n; i++)
  {
    temp[i] = 0.0;
    for (int j = 0; j < n; j++)
    {
      temp[i] += B[i][j] * rhs_vector_copy[j];
    }
  }

  for (int i = 0; i < n; i++)
  {
    rhs_vector_copy[i] = temp[i];
  }

  free(temp);
}

void Get_ObjectiveFunction(Model *model, double *rhs_vector)
{
  size_t n = model->num_constraints;

  for (int i = 0; i < n; i++)
  {
    int basic_idx = model->basics_vector[i];

    // Infeasibility check!
    if (fabs(model->coeffs[basic_idx]) == model->BIG_M  && rhs_vector[i] > 1e-6)
    {
      printf("Model Infeasible. Artificial basic variable has a positive RHS value. Terminating!\n");
      exit(0);
    }

    model->objective_function += model->coeffs[basic_idx] * rhs_vector[i];

    if (fabs(model->coeffs[basic_idx]) > 1e-6 &&
      fabs(model->coeffs[basic_idx]) < model->BIG_M &&
      rhs_vector[i] > 1e-6)
    {

      // Convert back to original coefficient if it's a minimization problem originally
      double original_coeff = model->coeffs[basic_idx];
      if (strcmp(model->objective, "MINIMIZE") == 0)
      {
        original_coeff *= -1;
      }

      printf("Value of variable x%i is %f with coefficient %f\n",
             basic_idx, rhs_vector[i], original_coeff);
    }
  }

  // Same logic applies
  if (strcmp(model->objective, "MINIMIZE") == 0)
  {
    model->objective_function *= -1;
  }

  printf("Optimal solution found! Objective value: %f\n", model->objective_function);
}



void FreeModel(Model *model)
{
  size_t size = 0;

  size += sizeof(char) * 8; // objective string
  size += sizeof(int) * model->num_constraints * 2;
  size += sizeof(int) * model->num_vars;

  int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;

  for (int i = 0; i < model->num_constraints; i++)
  {
    free(model->lhs_matrix[i]);
    size += sizeof(double) * total_cols;
  }

  size += sizeof(double) * total_cols;             // Coefficients
  size += sizeof(double) * model->num_constraints; // RHS vector
  size += sizeof(int) * model->num_constraints;    // Basics vector
  size += sizeof(double);                          // Objective function
  size += sizeof(char) * model->num_constraints;   // Constraints symbols
  size += sizeof(int) * 2;                         // Inequalities and equalities count
  size += sizeof(int) * model->artificials_count;   // Equalities vector
  size += sizeof(int);                             // Solver iterations
  size += sizeof(int) * model->non_basics_count;
  size += sizeof(double); // Big-M definition

  size_t n = model->num_constraints;
  size_t inversion_memory = n * (2 * n) * sizeof(double) + n * sizeof(double*); // Thanks Gauss!
  
  printf("Model's estimated memory usage: %zu bytes (%.2f KB)\n", size, size / 1024.0);
  printf("Matrix inversion memory usage per iteration: %zu bytes (%.2f KB)\n", 
         inversion_memory, inversion_memory / 1024.0);
  printf("Total iterations: %d\n", model->solver_iterations);

  free(model->lhs_matrix);
  free(model->objective);
  free(model->coeffs);
  free(model->rhs_vector);
  free(model->basics_vector);
  free(model->constraints_symbols);
  free(model->artificials_vector);
  free(model->non_basics);
  free(model);


  printf("-------------------------------\n");
  printf("Model free'd from the heap\n");
}

void ValidateModelPointers(Model *model)
{
  if (!model)
  {
    fprintf(stderr, "Fatal Error: Model pointer is NULL.\n");
    exit(1);
  }

  if (!model->lhs_matrix)
  {
    fprintf(stderr, "Fatal Error: model->lhs_matrix is NULL.\n");
    exit(1);
  }

  if (!model->objective)
  {
    fprintf(stderr, "Fatal Error: model->objective is NULL.\n");
    exit(1);
  }

  if (!model->coeffs)
  {
    fprintf(stderr, "Fatal Error: model->coeffs is NULL.\n");
    exit(1);
  }

  if (!model->rhs_vector)
  {
    fprintf(stderr, "Fatal Error: model->rhs_vector is NULL.\n");
    exit(1);
  }

  if (!model->basics_vector)
  {
    fprintf(stderr, "Fatal Error: model->basics_vector is NULL.\n");
    exit(1);
  }

  if (!model->constraints_symbols)
  {
    fprintf(stderr, "Fatal Error: model->constraints_symbols is NULL.\n");
    exit(1);
  }

  if (!model->artificials_vector)
  {
    fprintf(stderr, "Fatal Error: model->artificials_vector is NULL.\n");
    exit(1);
  }

  if (!model->non_basics)
  {
    fprintf(stderr, "Fatal Error: model->non_basics is NULL.\n");
    exit(1);
  }

  for (int i = 0; i < model->num_constraints; i++)
  {
    if (!model->lhs_matrix[i])
    {
      fprintf(stderr, "Fatal Error: model->lhs_matrix[%d] is NULL.\n", i);
      exit(1);
    }
  }
}

