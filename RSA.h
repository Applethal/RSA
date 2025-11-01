#ifndef RSA_H
#define RSA_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> // Do I need this to set the Big Ms? lol
typedef struct {
  char *objective;   // MINIMIZE or MAXIMIZE
  int num_constraints; // Number of constraints
  int num_vars;       // Number of variables
  double **columns;   // Constraints x variables coeffs matrix  
  double *coeffs;     // Variable coefficients
  double *rhs_vector; // Right hand side vector of the constraints 
  int *basics_vector; // Vector where I keep track of the basic variables, will contain only the inequalities' indices in the columns matrix
  double objective_function; // Final objective function 
  char *constraints_symbols; // Tracks the constraints' symbols, Used for debugging only
  int inequalities_count;// Counts the number of inequalities in the model, will be used for debugging only
  int equalities_count; // Counts the number of equality constraints in the model, used for debugging only 
  int *equalities_vector; // Contains the indices of the equality constraints in the model
  int solver_iterations; // Counts how many iterations it took for the solver to finalize the processing
  int *non_basics; // Contains indices of non-basic variables
  int non_basics_count; // Counts the number of non-basic variables 
} Model;

Model* ReadCsv(FILE *textfile);

void FreeModel(Model *model);

void InvertMatrix(double **matrix, int n);

void RevisedSimplex(Model* model);

double** Get_BasicsMatrix(Model* model);

void PrintColumns(Model* model);

double Get_ReducedPrice(Model* model, double **B_inv, int var_col, double *multiplier_vector);

double* Get_SimplexMultiplier(Model* model, double **B_inv);

double* Get_pivot_column(double** B_inv, Model* model, int best_cost_idx);

void UpdateRhs(Model* model, double* rhs_vector_copy, double** B);

void Get_ObjectiveFunction(Model* model, double *rhs_vector); 

void RevisedSimplex_Debug(Model* model);

Model* ReadCsv(FILE *csvfile) {
  Model *model = (Model*)malloc(sizeof(Model));
  char line[1024];

  // Line 1: Read objective type (MINIMIZE or MAXIMIZE)
  if (!fgets(line, sizeof(line), csvfile)) {
    fprintf(stderr, "Error: Could not read objective type\n");
    free(model);
    return NULL;
  }

  // Remove newline
  line[strcspn(line, "\r\n")] = 0;

  model->objective = (char*)malloc(strlen(line) + 1);
  strcpy(model->objective, line);

  // Line 2: Read objective coefficients and count variables
  if (!fgets(line, sizeof(line), csvfile)) {
    fprintf(stderr, "Error: Could not read objective coefficients\n");
    free(model->objective);
    free(model);
    return NULL;
  }

  // Count variables by counting commas + 1
  model->num_vars = 1;
  for (char *p = line; *p; p++) {
    if (*p == ',') model->num_vars++;
  }

  // Parse objective coefficients
  double *objective_coeffs = (double*)malloc(model->num_vars * sizeof(double));
  char *token = strtok(line, ",");
  int idx = 0;
  while (token != NULL && idx < model->num_vars) {
    objective_coeffs[idx++] = atof(token);
    token = strtok(NULL, ",");
  }

  // First pass: count constraints and their types
  long constraint_start = ftell(csvfile);
  model->num_constraints = 0;
  int num_le = 0, num_ge = 0, num_eq = 0;

  while (fgets(line, sizeof(line), csvfile)) {
    if (strstr(line, "<=")) num_le++;
    else if (strstr(line, ">=")) num_ge++;
    else if (strstr(line, "=")) num_eq++;
    model->num_constraints++;
  }

  int num_slack_surplus = num_le + num_ge;
  int num_artificial = num_eq + num_ge;
  int total_cols = model->num_vars + num_slack_surplus + num_artificial;

  // Allocate model arrays
  model->columns = (double**)malloc(model->num_constraints * sizeof(double*));
  for (int i = 0; i < model->num_constraints; i++) {
    model->columns[i] = (double*)malloc(total_cols * sizeof(double));
    for (int j = 0; j < total_cols; j++) {
      model->columns[i][j] = 0.0;
    }
  }

  model->rhs_vector = (double*)malloc(model->num_constraints * sizeof(double));
  model->constraints_symbols = (char*)malloc(model->num_constraints * sizeof(char));
  model->equalities_vector = (int*)malloc(num_artificial * sizeof(int));
  model->basics_vector = (int*)malloc(model->num_constraints * sizeof(int));
  model->coeffs = (double*)malloc(total_cols * sizeof(double));

  model->inequalities_count = num_slack_surplus;
  model->equalities_count = num_artificial;

  // Second pass: read constraints
  fseek(csvfile, constraint_start, SEEK_SET);

  int slack_col = model->num_vars;
  int artificial_col = model->num_vars + num_slack_surplus;
  int slack_idx = 0;
  int artificial_idx = 0;

  for (int i = 0; i < model->num_constraints; i++) {
    if (!fgets(line, sizeof(line), csvfile)) {
      fprintf(stderr, "Error: Unexpected end of file at constraint %d\n", i);
      // TODO: Cleanup and return NULL
      break;
    }

    // Parse constraint coefficients
    token = strtok(line, ",");
    int col_idx = 0;

    // Read variable coefficients
    while (token != NULL && col_idx < model->num_vars) {
      model->columns[i][col_idx] = atof(token);
      col_idx++;
      token = strtok(NULL, ",");
    }

    // Read operator (<=, >=, =)
    if (token == NULL) {
      fprintf(stderr, "Error: Missing operator in constraint %d\n", i);
      continue;
    }

    // Remove spaces from operator token
    char op[4] = "";
    int op_idx = 0;
    for (char *p = token; *p && op_idx < 3; p++) {
      if (*p != ' ' && *p != '\t') {
        op[op_idx++] = *p;
      }
    }
    op[op_idx] = '\0';

    // Read RHS
    token = strtok(NULL, ",");
    if (token == NULL) {
      fprintf(stderr, "Error: Missing RHS in constraint %d\n", i);
      continue;
    }
    model->rhs_vector[i] = atof(token);

    // Process constraint type
    if (strcmp(op, "<=") == 0) {
      model->constraints_symbols[i] = 'L';
      int col = slack_col + slack_idx;
      model->columns[i][col] = 1.0;
      model->basics_vector[i] = col;
      model->coeffs[col] = 0.0;
      slack_idx++;

    } else if (strcmp(op, ">=") == 0) {
      model->constraints_symbols[i] = 'G';
      int surplus_col = slack_col + slack_idx;
      int artif_col = artificial_col + artificial_idx;

      model->columns[i][surplus_col] = -1.0;
      model->columns[i][artif_col] = 1.0;
      model->equalities_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;
      model->coeffs[artif_col] = 9999.00;
      slack_idx++;
      artificial_idx++;

    } else if (strcmp(op, "=") == 0) {
      model->constraints_symbols[i] = 'E';
      int artif_col = artificial_col + artificial_idx;
      model->columns[i][artif_col] = 1.0;
      model->equalities_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;
      model->coeffs[artif_col] = 9999.00;
      artificial_idx++;

    } else {
      fprintf(stderr, "Error: Invalid operator '%s' in constraint %d\n", op, i);
    }
  }

  // printf("Number of slack/surplus variables: %d\n", model->inequalities_count);
  // printf("Number of artificial variables: %d\n", model->equalities_count);

  if (model->num_constraints > model->num_vars) {
    printf("Warning: More constraints than variables! The Revised Simplex Algorithm works best when the number of variables exceeds the number of constraints\n");
  }

  model->solver_iterations = 1;

  // Set objective coefficients
  for (int i = 0; i < model->num_vars; i++) {
    model->coeffs[i] = objective_coeffs[i];
  }

  // Build non-basics list
  model->non_basics = (int*)malloc(total_cols * sizeof(int));
  int non_basic_idx = 0;
  for (int col = 0; col < total_cols; col++) {
    int is_basic = 0;
    for (int j = 0; j < model->num_constraints; j++) {
      if (model->basics_vector[j] == col) {
        is_basic = 1;
        break;
      }
    }
    if (!is_basic) {
      model->non_basics[non_basic_idx++] = col;
    }
  }
  model->non_basics_count = non_basic_idx;

  model->objective_function = 0.0;
  printf("Objective function mode: %s\n", model->objective);
  printf("Number of variables: %d\n", model->num_vars);
  printf("Number of constraints: %d\n", model->num_constraints);
  printf("Objective coefficients (Slack and Artificial coeffs included):");
  for (int i = 0; i < model->num_vars + model->equalities_count + model->inequalities_count; i++) {
    printf(" %.1f", model->coeffs[i]);  
  }
  printf("\n\n");



  // Converting the problem to its opposite objective direction
  if (strcmp(model->objective, "MINIMIZE") == 0) {
    for (int i = 0; i < model->num_vars; i++) {
      model->coeffs[i] *= -1;
    }
    strcpy(model->objective, "MAXIMIZE");
  } else {
    for (int i = 0; i < model->num_vars; i++) {
      model->coeffs[i] *= -1;
    }
    strcpy(model->objective, "MINIMIZE");
  }

  // Apply Big M to artificial variables
  int artificial_start = model->num_vars + model->inequalities_count;
  double BIG_M = 9999.0;
  for (int i = 0; i < model->equalities_count; i++) {
    model->coeffs[artificial_start + i] = -BIG_M;
  }

  free(objective_coeffs);
  return model;
}

void PrintColumns(Model* model){
  int total_cols = model->num_vars + model->inequalities_count + model->equalities_count;

  printf("Constraints in canonical form:\n");
  for (int i = 0; i < model->num_constraints; i++) {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < total_cols; j++) {
      printf("%5.1f ", model->columns[i][j]);  
    }
    printf("| RHS: %5.1f", model->rhs_vector[i]);
    printf("\n");
    // printf(" | Basic var: col %d\n", model->basics_vector[i]);
  }
  printf("\n");

}


void InvertMatrix(double **matrix, int n) {
  double **augmented = (double**)malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++) {
    augmented[i] = (double*)malloc(2 * n * sizeof(double));
    for (int j = 0; j < n; j++) {
      augmented[i][j] = matrix[i][j];
    }
    for (int j = n; j < 2 * n; j++) {
      augmented[i][j] = (j - n == i) ? 1.0 : 0.0;
    }
  }

  // Gauss-Jordan with partial pivoting
  for (int i = 0; i < n; i++) {
    // Find row with largest absolute value in column i
    int max_row = i;
    double max_val = fabs(augmented[i][i]);
    for (int k = i + 1; k < n; k++) {
      if (fabs(augmented[k][i]) > max_val) {
        max_val = fabs(augmented[k][i]);
        max_row = k;
      }
    }

    // Swap rows if needed
    if (max_row != i) {
      double *temp = augmented[i];
      augmented[i] = augmented[max_row];
      augmented[max_row] = temp;
    }

    double pivot = augmented[i][i];
    if (fabs(pivot) < 1e-10) {
      printf("ERROR: Singular matrix - basis is not invertible!\n");
      printf("Check your basis selection at iteration %d\n", i);
      for (int k = 0; k < n; k++) {
        free(augmented[k]);
      }
      free(augmented);
      exit(1);
    }

    // Scale pivot row
    for (int j = 0; j < 2 * n; j++) {
      augmented[i][j] /= pivot;
    }

    // Eliminate column
    for (int k = 0; k < n; k++) {
      if (k != i) {
        double factor = augmented[k][i];
        for (int j = 0; j < 2 * n; j++) {
          augmented[k][j] -= factor * augmented[i][j];
        }
      }
    }
  }

  // Copy inverse back
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      matrix[i][j] = augmented[i][n + j];
    }
  }

  for (int i = 0; i < n; i++) {
    free(augmented[i]);
  }
  free(augmented);
}

double** Get_BasicsMatrix(Model* model) {

  double **Basics_matrix = (double**)malloc(model->num_constraints* sizeof(double));



  // printf("Basics matrix:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    Basics_matrix[i] = (double*)malloc(model->num_constraints * sizeof(double));
    for (size_t j = 0; j < model->num_constraints; j++) {
      int basics_index = model->basics_vector[j];
      Basics_matrix[i][j] = model->columns[i][basics_index]; 
      //    printf(" %f ", Basics_matrix[i][j]);
    }
    //    printf("\n");
  }

  return Basics_matrix;
}

void Print_Basics_matrix(double **B, int size) {
  printf("Basis Matrix (B):\n");
  for (int i = 0; i < size; i++) {

    for (int j = 0; j < size; j++) {
      printf("%6.2f ", B[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}
void RevisedSimplex_Debug(Model* model){
  int termination = 0;
  int n = model->num_constraints;
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

  while (termination != 1) {
    int feasibility_check = 0;
    printf("Beginning solver iteration %i ... \n", model->solver_iterations);
    printf("Getting the B matrix:\n");
    double **B = Get_BasicsMatrix(model);
    printf("===================================================================\n");

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        printf(" %f ", B[i][j]); 
      } 
      printf("\n");
    }

    printf("\n");
    printf("===================================================================\n");

    double original_RHS[n];
    for (size_t i = 0; i < n; i++) {
      original_RHS[i] = model->rhs_vector[i]; 
    }
    printf("Inverting the B matrix now:\n");
    InvertMatrix(B, model->num_constraints);
    printf("===================================================================\n");

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        printf(" %f ", B[i][j]); 
      } 
      printf("\n");
    }

    printf("\n");
    printf("===================================================================\n");





    printf("Getting the simplex multiplier P vector:\n");
    double *Simplex_multiplier = Get_SimplexMultiplier(model, B);
    printf("P vector elements:\n");

    for (size_t i = 0; i < n; i++) {
      printf(" %f ", Simplex_multiplier[i]); 
    }
    printf("\n");
    printf("===================================================================\n");

    int entering_var_idx = 0;
    int entering_var = 0;
    UpdateRhs(model, original_RHS, B); 

    printf("Updating the Right hand side vector\n");

    for (size_t i = 0; i < n; i++) {
      printf(" %f ",model->rhs_vector[i]);
    }
    printf("\n");

    double best_reduced_cost;
    if (strcmp(model->objective, "MINIMIZE") == 0) {
      best_reduced_cost = DBL_MAX;
      printf("MINIMIZATION problem detected, choosing the lowest reduced cost value\n");
    }
    else {
      best_reduced_cost = -DBL_MAX;
      printf("MAXIMIZATION problem detected, choosing the highest reduced cost value\n");
    }
    printf("===================================================================\n");

    printf("Getting the entering variable:\n");
    printf("===================================================================\n");
    // Getting the entering variable
    for (size_t i = 0; i < model->non_basics_count; i++) {

      int non_basic_idx = model->non_basics[i]; 


      double reduced_cost = Get_ReducedPrice(model, B, non_basic_idx, Simplex_multiplier ) ;

      printf("The reduced cost of non basic variable %i is %f \n", non_basic_idx, reduced_cost);

      if (strcmp(model->objective, "MAXIMIZE") == 0) {

        if (reduced_cost > best_reduced_cost) {
          best_reduced_cost = reduced_cost;
          entering_var_idx = i;
          entering_var = non_basic_idx;
        }
        if (reduced_cost <= 0) {
          feasibility_check++;
        }
      } else if (strcmp(model->objective, "MINIMIZE") == 0) {
        if (reduced_cost < best_reduced_cost) {
          best_reduced_cost = reduced_cost;
          entering_var_idx = i;
          entering_var = non_basic_idx;
        }
        if (reduced_cost >= 0) {
          feasibility_check++;
        }

      }

    }
    printf("Entering variable index: %i, cost: %f \n", entering_var, best_reduced_cost);
    // If all the reduced costs are positive in a minimization problem / negative in a maximization problem => optimal solution found
    printf("===================================================================\n");

    if (feasibility_check == model->non_basics_count) {
      printf("Optimal solution found, terminating!\n");
      termination++;
      Get_ObjectiveFunction(model, original_RHS);
      break;
    }
    // printf("Best reduced cost is %f for the entering variable %i \n", best_reduced_cost, entering_var);
    double *Pivot = Get_pivot_column(B, model, entering_var);
    //UpdateRhs(model, original_RHS, B); 
    // Getting the exiting variable   
    printf("Getting the exiting variable:\n");
    printf("===================================================================\n");

    int exiting_var_idx = 0;

    int exiting_var = 0;
    double best_ratio = DBL_MAX;
    //
    for (size_t i = 0; i < model->num_constraints; i++) {
      if (Pivot[i] == 0.0) {
        continue;
      }
      double ratio = original_RHS[i] / Pivot[i];
      printf("Ratio of variable %i is %f which has a pivot value of %f \n", model->basics_vector[i], ratio, Pivot[i]);
      if (ratio < best_ratio && ratio > 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];

      }

    }
    if (best_ratio < 0) {
      printf("The lowest ratio is %f, LP is unbounded! Terminating!\n", best_ratio);
      exit(0); 
    }
    printf("The exiting variable is %i with ratio of %f \n", exiting_var, best_ratio, exiting_var_idx );


    printf("===================================================================\n");

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;
    // model->basics_vector[entering_var_idx] = entering_var;
    // model->non_basics[exiting_var_idx] = exiting_var;
    // //
    printf("Entering variable:%i at position %i \n", entering_var, entering_var_idx);
    printf("Exiting variable:%i at position: %i \n",exiting_var,exiting_var_idx);
    printf("Non basics updated:\n");
    for (size_t i = 0; i < model->non_basics_count; i++) {
      printf(" %i ", model->non_basics[i]); 
    }
    printf("\n");
    printf("Basics vector:\n");

    for (size_t i = 0; i < model->num_constraints; i++) {
      printf(" %i ", model->basics_vector[i]); 
    }
    // termination++;
    model->solver_iterations++;
    free(Pivot);
    for (size_t i = 0; i < model->num_constraints; i++) {
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



void RevisedSimplex(Model* model){
  int termination = 0;
  int n = model->num_constraints;
  // printf("\n");
  //   printf("Non basics start:\n");
  //   for (size_t i = 0; i < model->non_basics_count; i++) {
  //     printf(" %i ", model->non_basics[i] + 1); 
  //   }
  //   printf("\n");
  //   printf("Basics vector:\n");
  //
  //   for (size_t i = 0; i < model->num_constraints; i++) {
  //     printf(" %i ", model->basics_vector[i] + 1); 
  //   }
  //   printf("\n");
  while (termination != 1) {
    int feasibility_check = 0;
    printf("Beginning solver iteration %i ... \n", model->solver_iterations);
    double **B = Get_BasicsMatrix(model);
    double original_RHS[n];
    for (size_t i = 0; i < n; i++) {
      original_RHS[i] = model->rhs_vector[i]; 
    }

    InvertMatrix(B, model->num_constraints);

    double *Simplex_multiplier = Get_SimplexMultiplier(model, B);

    int entering_var_idx = 0;
    int entering_var = 0;
    UpdateRhs(model, original_RHS, B); 

    double best_reduced_cost;
    if (strcmp(model->objective, "MINIMIZE") == 0) {
      best_reduced_cost = DBL_MAX;
      // printf("MINIMIZATION problem detected, choosing the lowest reduced cost value\n");
    }
    else {
      best_reduced_cost = -DBL_MAX;
      // printf("MAXIMIZATION problem detected, choosing the highest reduced cost value\n");
    }



    // Getting the entering variable
    for (size_t i = 0; i < model->non_basics_count; i++) {

      int non_basic_idx = model->non_basics[i]; 


      double reduced_cost = Get_ReducedPrice(model, B, non_basic_idx, Simplex_multiplier ) ;

      // printf("reduced cost of non basic variable %i is %f \n", non_basic_idx, reduced_cost);

      if (strcmp(model->objective, "MAXIMIZE") == 0) {

        if (reduced_cost > best_reduced_cost) {
          best_reduced_cost = reduced_cost;
          entering_var_idx = i;
          entering_var = non_basic_idx;
        }
        if (reduced_cost <= 0) {
          feasibility_check++;
        }
      } else if (strcmp(model->objective, "MINIMIZE") == 0) {
        if (reduced_cost < best_reduced_cost) {
          best_reduced_cost = reduced_cost;
          entering_var_idx = i;
          entering_var = non_basic_idx;
        }
        if (reduced_cost >= 0) {
          feasibility_check++;
        }

      }

    }
    // If all the reduced costs are positive in a minimization problem / negative in a maximization problem => optimal solution found
    if (feasibility_check == model->non_basics_count) {
      printf("Optimal solution found, terminating!\n");
      termination++;
      Get_ObjectiveFunction(model, original_RHS);
      break;
    }
    // printf("Best reduced cost is %f for the entering variable %i \n", best_reduced_cost, entering_var);
    double *Pivot = Get_pivot_column(B, model, entering_var);
    //UpdateRhs(model, original_RHS, B); 
    // Getting the exiting variable   
    int exiting_var_idx = 0;
    int exiting_var = 0;
    double best_ratio = DBL_MAX;
    //
    for (size_t i = 0; i < model->num_constraints; i++) {
      if (Pivot[i] == 0.0) {
        continue;
      }
      double ratio = original_RHS[i] / Pivot[i];
      // printf("Ratio of variable %i is %f which has a pivot value of %f \n", model->basics_vector[i], ratio, Pivot[i]);
      if (ratio < best_ratio && ratio > 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];

      }

    }
    if (best_ratio < 0) {
      printf("LP is unbounded! Terminating!\n");
      exit(0); 
    }
    // printf("Best exiting variable is %i with ratio of %f at position %i \n", exiting_var, best_ratio, exiting_var_idx );
    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;
    // model->basics_vector[entering_var_idx] = entering_var;
    // model->non_basics[exiting_var_idx] = exiting_var;
    // //
    // printf("Entering variable:%i at position %i \n", entering_var + 1, entering_var_idx);
    // printf("Exiting variable:%i at position: %i \n",exiting_var + 1,exiting_var_idx);
    // printf("Non basics updated:\n");
    // for (size_t i = 0; i < model->non_basics_count; i++) {
    //   printf(" %i ", model->non_basics[i] + 1); 
    // }
    // printf("\n");
    // printf("Basics vector:\n");
    //
    // for (size_t i = 0; i < model->num_constraints; i++) {
    //   printf(" %i ", model->basics_vector[i] + 1); 
    // }
    // termination++;
    model->solver_iterations++;
    free(Pivot);
    for (size_t i = 0; i < model->num_constraints; i++) {
      free(B[i]); 
    }
    free(B);
    free(Simplex_multiplier);


  }
}

double Get_ReducedPrice(Model* model, double **B_inv, int var_col, double *multiplier_vector) {
  int n = model->num_constraints;


  double dot_product = 0.0;


  // printf("column vector for the variable %i :\n", var_col);
  // printf("(");
  for (int i = 0; i < n; i++) {
    // printf(" %f ", model->columns[i][var_col]);
    dot_product += multiplier_vector[i] * model->columns[i][var_col];
  }
  // printf(")\n");
  // printf("Dot product: %f \n", dot_product);
  // printf("Objective coeff of variable %i is %f \n", var_col, model->coeffs[var_col]);
  double reduced_cost = model->coeffs[var_col] - dot_product; 
  // printf("Reduced cost of variable %i is %f \n", var_col, reduced_cost);
  return reduced_cost;
}
double* Get_SimplexMultiplier(Model* model, double **B_inv){
  int n = model->num_constraints;
  // THIS IS CORRECT AS Of THE VIDEO

  double *multiplier_vector = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      int basic_col_idx = model->basics_vector[j];
      sum += model->coeffs[basic_col_idx] * B_inv[j][i];
    }
    multiplier_vector[i] = sum;
  }

  // printf("multiplier vector:\n");
  // printf("(");
  // for (size_t i = 0; i < n; i++) {
  //   printf(" %f ", multiplier_vector[i]);
  //
  // }
  //
  // printf(")\n");
  //
  return multiplier_vector;

}

double* Get_pivot_column(double** B_inv, Model* model, int best_cost_idx) {
  int n = model->num_constraints;
  double* Pivot = (double*) malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      sum += B_inv[i][j] * model->columns[j][best_cost_idx];
    }
    Pivot[i] = sum;

  }
  return Pivot;
}

void UpdateRhs(Model* model, double *rhs_vector_copy, double** B) {
  int n = model->num_constraints;
  double* temp = (double*)malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    temp[i] = 0.0;
    for (int j = 0; j < n; j++) {
      temp[i] += B[i][j] * rhs_vector_copy[j];


    }
  }


  for (int i = 0; i < n; i++) {
    if (temp[i] < 0) {
      printf("Negative RHS value detected, LP model is infeasible. Terminating!\n");
      exit(0);
    }

    rhs_vector_copy[i] = temp[i];
  }

  free(temp);
  //
  // printf("Updated RHS vector:\n");
  // for (int i = 0; i < n; i++) {
  //     printf("RHS[%d] = %f\n", i, rhs_vector[i]);
  // }
}
void Get_ObjectiveFunction(Model* model, double *rhs_vector) {
  int n = model->num_constraints;
  for (int i = 0; i < n; i++) {
    int basic_idx = model->basics_vector[i];
    if (model->coeffs[basic_idx] > 9000.0 && rhs_vector[i] > 0) {
      printf("Model Infeasible. Terminating!\n");
      exit(0);

    }


    if (model->coeffs[basic_idx] == 0 ) {
      continue;

    }

    double value = (model->coeffs[basic_idx] * rhs_vector[i]) * -1;
    model->objective_function += value;

    if (model->objective_function > 9000 * 0.5) {
      printf("Model Infeasible. Artificial variable x%i has a positive RHS value Terminating!\n", basic_idx);
      exit(0);
    }
    if ( rhs_vector[i] == 0 ) {
      continue;
    }
    printf("Value of the variable x%i is %f with a coefficient of %f \n", basic_idx +1, rhs_vector[i], model->coeffs[basic_idx] * -1);



  }

}
void FreeModel(Model* model){
  for (int i = 0; i < model->num_constraints; i++) {
    free(model->columns[i]);
  }
  free(model->columns);

  free(model->coeffs);
  free(model->rhs_vector);
  free(model->basics_vector);
  free(model->constraints_symbols);
  free(model->equalities_vector);
  free(model->non_basics);

  printf("-------------------------------\n");
  printf("Model free'd from the heap\n"); 
}

#endif 


