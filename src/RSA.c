#include "RSA.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> // Do I really need this to setup the Big Ms? lol

Model *ReadCsv(FILE *csvfile)
{
  Model *model = (Model *)malloc(sizeof(Model));
  char line[1024];

  // Line 1: Read objective type (MINIMIZE or MAXIMIZE)
  if (!fgets(line, sizeof(line), csvfile))
  {
    fprintf(stderr, "Error: Could not read objective type\n");
    free(model);
    return NULL;
  }

  // Remove newline
  line[strcspn(line, "\r\n")] = 0;

  // Storing objective
  model->objective = (char *)malloc(strlen(line) + 1);
  strcpy(model->objective, line);

  // Line 2: Read objective coefficients and count variables
  if (!fgets(line, sizeof(line), csvfile))
  {
    fprintf(stderr, "Error: Could not read objective coefficients\n");
    free(model->objective);
    free(model);
    return NULL;
  }

  // Count variables by counting commas + 1
  model->num_vars = 1;
  for (char *p = line; *p; p++)
  {
    if (*p == ',')
      model->num_vars++;
  }

  // Parse objective coefficients
  double *objective_coeffs = (double *)malloc(model->num_vars * sizeof(double));
  char *token = strtok(line, ",");
  int idx = 0;
  while (token != NULL && idx < model->num_vars)
  {
    objective_coeffs[idx++] = atof(token);
    token = strtok(NULL, ",");
  }

  // First pass: count constraints and their types
  long constraint_start = ftell(csvfile);
  model->num_constraints = 0;
  int num_le = 0, num_ge = 0, num_eq = 0;

  while (fgets(line, sizeof(line), csvfile))
  {
    if (strstr(line, "<="))
      num_le++;
    else if (strstr(line, ">="))
      num_ge++;
    else if (strstr(line, "="))
      num_eq++;
    model->num_constraints++;
  }

  int num_slack_surplus = num_le + num_ge;
  int num_artificial = num_eq + num_ge;
  int total_cols = model->num_vars + num_slack_surplus + num_artificial;

  // Allocate model arrays
  model->columns = (double **)malloc(model->num_constraints * sizeof(double *));
  for (int i = 0; i < model->num_constraints; i++)
  {
    model->columns[i] = (double *)malloc(total_cols * sizeof(double));
    for (int j = 0; j < total_cols; j++)
    {
      model->columns[i][j] = 0.0;
    }
  }

  model->rhs_vector = (double *)malloc(model->num_constraints * sizeof(double));
  model->constraints_symbols = (char *)malloc(model->num_constraints * sizeof(char));
  model->equalities_vector = (int *)malloc(num_artificial * sizeof(int));
  model->basics_vector = (int *)calloc(model->num_constraints, sizeof(int));
  model->coeffs = (double *)calloc(total_cols, sizeof(double));

  model->inequalities_count = num_slack_surplus;
  model->equalities_count = num_artificial;

  // Second pass: read constraints
  fseek(csvfile, constraint_start, SEEK_SET);

  int slack_col = model->num_vars;
  int artificial_col = model->num_vars + num_slack_surplus;
  int slack_idx = 0;
  int artificial_idx = 0;

  for (int i = 0; i < model->num_constraints; i++)
  {
    if (!fgets(line, sizeof(line), csvfile))
    {
      fprintf(stderr, "Error: Unexpected end of file at constraint %d\n", i);
      break;
    }

    // Parse constraint coefficients
    token = strtok(line, ",");
    int col_idx = 0;

    // Read variable coefficients
    while (token != NULL && col_idx < model->num_vars)
    {
      model->columns[i][col_idx] = atof(token);
      col_idx++;
      token = strtok(NULL, ",");
    }

    // Read operator (<=, >=, =)
    if (token == NULL)
    {
      fprintf(stderr, "Error: Missing operator or variable coefficient in constraint %d\n, make sure all the variables have exactly %i variables", i, model->num_vars);
      exit(0);
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

    // Read RHS
    token = strtok(NULL, ",");
    if (token == NULL)
    {
      fprintf(stderr, "Error: Missing RHS in constraint %d\n", i);
      exit(0);
    }
    model->rhs_vector[i] = atof(token);

    // Process constraint type
    if (strcmp(op, "<=") == 0)
    {
      model->constraints_symbols[i] = 'L';
      int col = slack_col + slack_idx;
      model->columns[i][col] = 1.0;
      model->basics_vector[i] = col;
      model->coeffs[col] = 0.0; // Slack variables have 0 coefficient
      slack_idx++;
    }
    else if (strcmp(op, ">=") == 0)
    {
      model->constraints_symbols[i] = 'G';
      int surplus_col = slack_col + slack_idx;
      int artif_col = artificial_col + artificial_idx;

      model->columns[i][surplus_col] = -1.0;
      model->columns[i][artif_col] = 1.0;
      model->equalities_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;
      // Big M will be set later
      slack_idx++;
      artificial_idx++;
    }
    else if (strcmp(op, "=") == 0)
    {
      model->constraints_symbols[i] = 'E';
      int artif_col = artificial_col + artificial_idx;
      model->columns[i][artif_col] = 1.0;
      model->equalities_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;
      // Big M will be set later
      artificial_idx++;
    }
    else
  {
      fprintf(stderr, "Error: Invalid operator '%s' in constraint %d\n", op, i);
    }
  }

  if (model->num_constraints > model->num_vars)
  {
    printf("Warning: More constraints than variables! The revised simplex algorithm works best when there are more variables than constraints!\n");
  }

  model->solver_iterations = 1;

  // Set objective coefficients
  for (int i = 0; i < model->num_vars; i++)
  {
    model->coeffs[i] = objective_coeffs[i];
  }

  // Applying Big M to artificial variables
  int artificial_start = model->num_vars + model->inequalities_count;
  double BIG_M = 9999.0;
  for (int i = 0; i < model->equalities_count; i++)
  {
    if (strcmp(model->objective, "MAXIMIZE") == 0)
    {
      model->coeffs[artificial_start + i] = -BIG_M; // Penalty for MAX
    }
    else
  {
      model->coeffs[artificial_start + i] = BIG_M; // Penalty for MIN
    }
  }

  // Converting MINIMIZE to MAXIMIZE by negating ALL coefficients
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

  model->objective_function = 0.0;

  free(objective_coeffs);
  return model;
}

void PrintColumns(Model *model)
{
  int total_cols = model->num_vars + model->inequalities_count + model->equalities_count;

  printf("Constraints in canonical form:\n");
  for (int i = 0; i < model->num_constraints; i++)
  {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < total_cols; j++)
    {
      printf("%5.1f ", model->columns[i][j]);
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

void InvertMatrix(double **matrix, int n)
{
  double **augmented = (double **)malloc(n * sizeof(double *));
  for (int i = 0; i < n; i++)
  {
    augmented[i] = (double *)malloc(2 * n * sizeof(double));
    for (int j = 0; j < n; j++)
    {
      augmented[i][j] = matrix[i][j];
    }
    for (int j = n; j < 2 * n; j++)
    {
      augmented[i][j] = (j - n == i) ? 1.0 : 0.0;
    }
  }

  // Gauss-Jordan with partial pivoting
  for (int i = 0; i < n; i++)
  {
    // Find row with largest absolute value in column i
    int max_row = i;
    double max_val = fabs(augmented[i][i]);
    for (int k = i + 1; k < n; k++)
    {
      if (fabs(augmented[k][i]) > max_val)
      {
        max_val = fabs(augmented[k][i]);
        max_row = k;
      }
    }

    // Swap rows if needed
    if (max_row != i)
    {
      double *temp = augmented[i];
      augmented[i] = augmented[max_row];
      augmented[max_row] = temp;
    }

    double pivot = augmented[i][i];
    if (fabs(pivot) < 1e-10)
    {
      printf("ERROR: Singular matrix - basis is not invertible!\n");
      for (int k = 0; k < n; k++)
      {
        free(augmented[k]);
      }
      free(augmented);
      exit(1);
    }

    // Scale pivot row
    for (int j = 0; j < 2 * n; j++)
    {
      augmented[i][j] /= pivot;
    }

    // Eliminate column
    for (int k = 0; k < n; k++)
    {
      if (k != i)
      {
        double factor = augmented[k][i];
        for (int j = 0; j < 2 * n; j++)
        {
          augmented[k][j] -= factor * augmented[i][j];
        }
      }
    }
  }

  // Copy inverse back
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      matrix[i][j] = augmented[i][n + j];
    }
  }

  for (int i = 0; i < n; i++)
  {
    free(augmented[i]);
  }
  free(augmented);
}

double **Get_BasicsMatrix(Model *model)
{
  double **Basics_matrix = (double **)malloc(model->num_constraints * sizeof(double *));

  for (size_t i = 0; i < model->num_constraints; i++)
  {
    Basics_matrix[i] = (double *)malloc(model->num_constraints * sizeof(double));
    for (size_t j = 0; j < model->num_constraints; j++)
    {
      int basics_index = model->basics_vector[j];
      Basics_matrix[i][j] = model->columns[i][basics_index];
    }
  }

  return Basics_matrix;
}

void RevisedSimplex(Model *model)
{
  int termination = 0;
  int n = model->num_constraints;
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1; // In theory, the worst case is n*m where m is the number of constraints and n the number of variables in the model

  while (termination != 1)
  {
    if (model->solver_iterations == MAX_ITERATIONS)
    {
      printf("Max iterations reached. Terminating!\n");
      exit(0);
    }

    int feasibility_check = 0;
    printf("Beginning solver iteration %i ... \n", model->solver_iterations);

    double **B = Get_BasicsMatrix(model);
    double original_RHS[n];
    for (size_t i = 0; i < n; i++)
    {
      original_RHS[i] = model->rhs_vector[i];
    }

    if (model->solver_iterations > 1)
    {
      InvertMatrix(B, model->num_constraints);
    }
    double *Simplex_multiplier = Get_SimplexMultiplier(model, B);

    int entering_var_idx = 0;
    int entering_var = 0;
    UpdateRhs(model, original_RHS, B);

    double best_reduced_cost = -DBL_MAX;

    // Getting the entering variable
    for (size_t i = 0; i < model->non_basics_count; i++)
    {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B, non_basic_idx, Simplex_multiplier);

      if (reduced_cost > best_reduced_cost)
      {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0)
      {
        feasibility_check++;
      }
    }

    // If all reduced costs are <= 0, optimal solution found
    if (feasibility_check == model->non_basics_count)
    {
      printf("Solver loop terminated!\n");
      termination++;
      Get_ObjectiveFunction(model, original_RHS);
      for (size_t i = 0; i < model->num_constraints; i++)
      {
        free(B[i]);
      }
      free(B);
      free(Simplex_multiplier);

      break;
    }

    double *Pivot = Get_pivot_column(B, model, entering_var);

    // Getting the exiting variable
    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = DBL_MAX;

    for (size_t i = 0; i < model->num_constraints; i++)
    {
      if (Pivot[i] <= 1e-6)
      {
        continue;
      }
      double ratio = original_RHS[i] / Pivot[i];
      if (ratio < best_ratio && ratio >= 0)
      {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }

    if (exiting_var_idx == -1)
    {
      printf("LP is unbounded! Terminating!\n");
      exit(0);
    }

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;
    model->solver_iterations++;

    free(Pivot);
    for (size_t i = 0; i < model->num_constraints; i++)
    {
      free(B[i]);
    }
    free(B);
    free(Simplex_multiplier);
  }
}

double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector)
{
  int n = model->num_constraints;
  double dot_product = 0.0;

  for (int i = 0; i < n; i++)
  {
    dot_product += multiplier_vector[i] * model->columns[i][var_col];
  }

  double reduced_cost = model->coeffs[var_col] - dot_product;
  return reduced_cost;
}

double *Get_SimplexMultiplier(Model *model, double **B_inv)
{
  int n = model->num_constraints;
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
  int n = model->num_constraints;
  double *Pivot = (double *)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++)
  {
    double sum = 0.0;
    for (int j = 0; j < n; j++)
    {
      sum += B_inv[i][j] * model->columns[j][best_cost_idx];
    }
    Pivot[i] = sum;
  }

  return Pivot;
}

void UpdateRhs(Model *model, double *rhs_vector_copy, double **B)
{
  int n = model->num_constraints;
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
  int n = model->num_constraints;
  double obj_value = 0.0;

  for (int i = 0; i < n; i++)
  {
    int basic_idx = model->basics_vector[i];

    // Infeasibility check!
    if (fabs(model->coeffs[basic_idx]) > 9000.0 && rhs_vector[i] > 1e-6)
    {
      printf("Model Infeasible. Artificial basic variable has a positive RHS value. Terminating!\n");
      exit(0);
    }

    obj_value += model->coeffs[basic_idx] * rhs_vector[i];

    if (fabs(model->coeffs[basic_idx]) > 1e-6 &&
      fabs(model->coeffs[basic_idx]) < 9000.0 &&
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
    obj_value *= -1;
  }

  model->objective_function = obj_value;
  printf("Optimal solution found! Objective value: %f\n", obj_value);
}

void RevisedSimplex_Debug(Model *model)
{
  int termination = 0;
  int n = model->num_constraints;
  printf("\n");
  printf("===================================================================\n");

  printf("Initial Non basic variable indices:\n");
  for (size_t i = 0; i < model->non_basics_count; i++)
  {
    printf(" %i ", model->non_basics[i]);
  }
  printf("\n");
  printf("Initial Basic variable indices:\n");

  for (size_t i = 0; i < model->num_constraints; i++)
  {
    printf(" %i ", model->basics_vector[i]);
  }
  printf("\n");
  printf("===================================================================\n");
  int max_iterations = (model->num_vars * model->num_constraints) + 1; // the worst case number of iterations in theory is n*m, where m is the number of constraints and n is the number  of variables in the model

  while (termination != 1)
  {
    int feasibility_check = 0;
    if (model->solver_iterations == max_iterations)
    {
      printf("Max iterations reached. Terminating!\n");
      exit(0);
    }

    printf("Beginning solver iteration %i ... \n", model->solver_iterations);
    printf("Getting the B matrix:\n");
    double **B = Get_BasicsMatrix(model);
    printf("===================================================================\n");

    for (size_t i = 0; i < n; i++)
    {
      for (size_t j = 0; j < n; j++)
      {
        printf(" %f ", B[i][j]);
      }
      printf("\n");
    }

    printf("\n");
    printf("===================================================================\n");

    double original_RHS[n];
    for (size_t i = 0; i < n; i++)
    {
      original_RHS[i] = model->rhs_vector[i];
    }
    if (model->solver_iterations > 1)
    {
      printf("Inverting the B matrix now:\n");
      InvertMatrix(B, model->num_constraints);
      printf("===================================================================\n");
      for (size_t i = 0; i < n; i++)
      {
        for (size_t j = 0; j < n; j++)
        {
          printf(" %f ", B[i][j]);
        }
        printf("\n");
      }

      printf("\n");
      printf("===================================================================\n");
    }
    else
  {
      printf("Solver iteration 1, skipping matrix inversion sequence\n");
    }

    printf("Getting the simplex multiplier P vector:\n");
    double *Simplex_multiplier = Get_SimplexMultiplier(model, B);
    printf("P vector elements:\n");

    for (size_t i = 0; i < n; i++)
    {
      printf(" %f ", Simplex_multiplier[i]);
    }
    printf("\n");
    printf("===================================================================\n");

    int entering_var_idx = 0;
    int entering_var = 0;
    UpdateRhs(model, original_RHS, B);

    printf("Updating the Right hand side vector\n");

    for (size_t i = 0; i < n; i++)
    {
      printf(" %f ", model->rhs_vector[i]);
    }
    printf("\n");

    double best_reduced_cost = -DBL_MAX;
    printf("MAXIMIZATION problem (after conversion), choosing the highest reduced cost value\n");

    printf("===================================================================\n");

    printf("Getting the entering variable:\n");
    printf("===================================================================\n");

    // Getting the entering variable
    for (size_t i = 0; i < model->non_basics_count; i++)
    {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B, non_basic_idx, Simplex_multiplier);

      printf("The reduced cost of non basic variable %i is %f \n", non_basic_idx, reduced_cost);

      if (reduced_cost > best_reduced_cost)
      {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0)
      {
        feasibility_check++;
      }
    }

    printf("Entering variable index: %i, cost: %f \n", entering_var, best_reduced_cost);
    printf("===================================================================\n");

    if (feasibility_check == model->non_basics_count)
    {
      printf("Solver loop terminated!\n");
      termination++;
      Get_ObjectiveFunction(model, original_RHS);
      for (size_t i = 0; i < model->num_constraints; i++)
      {
        free(B[i]);
      }
      free(B);
      free(Simplex_multiplier);


      break;

    }

    double *Pivot = Get_pivot_column(B, model, entering_var);

    printf("Getting the exiting variable:\n");
    printf("===================================================================\n");

    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = DBL_MAX;

    for (size_t i = 0; i < model->num_constraints; i++)
    {
      if (Pivot[i] <= 1e-6)
      {
        continue;
      }
      double ratio = original_RHS[i] / Pivot[i];
      printf("Ratio of variable %i is %f which has a pivot value of %f \n", model->basics_vector[i], ratio, Pivot[i]);
      if (ratio < best_ratio && ratio >= 0)
      {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }

    if (exiting_var_idx == -1)
    {
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
    for (size_t i = 0; i < model->non_basics_count; i++)
    {
      printf(" %i ", model->non_basics[i]);
    }
    printf("\n");
    printf("Basics vector:\n");

    for (size_t i = 0; i < model->num_constraints; i++)
    {
      printf(" %i ", model->basics_vector[i]);
    }

    model->solver_iterations++;
    free(Pivot);
    for (size_t i = 0; i < model->num_constraints; i++)
    {
      free(B[i]);
    }
    free(B);
    free(Simplex_multiplier);

    printf("\n");
    printf("Iteration %i ends here, press any key to continue!\n", model->solver_iterations - 1);
    printf("===================================================================\n");

    getchar();
  }
}

void FreeModel(Model *model)
{
  size_t size = 0;

  size += sizeof(char) * 8; // objective string
  size += sizeof(int) * model->num_constraints * 2;
  size += sizeof(int) * model->num_vars;

  int total_cols = model->num_vars + model->inequalities_count + model->equalities_count;

  for (int i = 0; i < model->num_constraints; i++)
  {
    free(model->columns[i]);
    size += sizeof(double) * total_cols;
  }

  size += sizeof(double) * total_cols;             // Coefficients
  size += sizeof(double) * model->num_constraints; // RHS vector
  size += sizeof(int) * model->num_constraints;    // Basics vector
  size += sizeof(double);                          // Objective function
  size += sizeof(char) * model->num_constraints;   // Constraints symbols
  size += sizeof(int) * 2;                         // Inequalities and equalities count
  size += sizeof(int) * model->equalities_count;   // Equalities vector
  size += sizeof(int);                             // Solver iterations
  size += sizeof(int) * model->non_basics_count;

  free(model->columns);
  free(model->objective);
  free(model->coeffs);
  free(model->rhs_vector);
  free(model->basics_vector);
  free(model->constraints_symbols);
  free(model->equalities_vector);
  free(model->non_basics);
  free(model);
  printf("Model's estimated memory usage: %zu bytes\n", size);
  printf("-------------------------------\n");
  printf("Model free'd from the heap\n");
}
  
void ValidateModelPointers(Model *model)
{
  if (!model)
  {
    fprintf(stderr, "Fatal Error: Model pointer is NULL.\n");
    exit(0);
  }

  if (!model->columns ||
    !model->objective ||
    !model->coeffs ||
    !model->rhs_vector ||
    !model->basics_vector ||
    !model->constraints_symbols ||
    !model->equalities_vector ||
    !model->non_basics)
  {
    fprintf(stderr, "Fatal Error: One or more model pointers are NULL.\n");
    exit(0);
  }

  for (int i = 0; i < model->num_constraints; i++)
  {
    if (!model->columns[i])
    {
      fprintf(stderr, "Fatal Error: columns[%d] is NULL.\n", i);
      exit(0);
    }
  }
}
