#ifndef RSA_H
#define RSA_H
#include <string.h>
#include <math.h>
#include <float.h> // Do I need this to set the Big Ms? lol
typedef struct {
  char *objective;   // MINIMIZE or MAXIMIZE
  int num_constraints; // Number of constraints
  int num_vars;       // Number of variables
  double **columns;   // Constraints x variables coeffs matrix  
  double *coeffs;     // Variable coefficients
  double *RHS_vector; // Right hand side vector of the constraints 
  int *Basics_vector; // Vector where I keep track of the basic variables, will contain only the inequalities' indices in the columns matrix
  double objective_function; // Final objective function 
  char *constraints_symbols; // Tracks the constraints' symbols, Used for debugging only
  int inequalities_count;// Counts the number of inequalities in the model, will be used for debugging only
  int equalities_count; // Counts the number of equality constraints in the model, used for debugging only 
  int *equalities_vector; // Contains the indices of the equality constraints in the model
  int solver_iterations; // Counts how many iterations it took for the solver to finalize the processing
  int *Non_basics; // Contains indices of non-basic variables
  int Non_basics_count; // Counts the number of non-basic variables 
} Model;

Model* read(FILE *textfile, char *objective);

void free_model(Model *model);

void invert_matrix(double **matrix, int n);

void RevisedSimplex(Model* model);

double** Get_Basics_Matrix(Model* model);

void Print_columns(Model* model);

double Get_reduced_price(Model* model, double **B_inv, int var_col, double *multiplier_vector);

double* Get_simplex_multiplier(Model* model, double **B_inv);

double* Get_pivot_column(double** B_inv, Model* model, int best_cost_idx);

void update_RHS(Model* model, double* RHS_vector_copy, double** B);

void Get_Objective_Function(Model* model, double *RHS_vector); 

Model* read(FILE *textfile, char *objective) {
  Model *model = (Model*)malloc(sizeof(Model));

  model->objective = objective;
  // Count the number of coefficients on the MINIMIZE line
  long pos = ftell(textfile);
  model->num_vars = 0;
  char ch;
  while ((ch = fgetc(textfile)) != ';' && ch != '\n' && ch != EOF) {
    if (ch == ' ' || ch == '\t') {
      continue;
    } else {
      ungetc(ch, textfile);
      double temp;
      if (fscanf(textfile, "%lf", &temp) == 1) {
        model->num_vars++;
      }
    }
  }

  // Go back and read the coefficients
  fseek(textfile, pos, SEEK_SET);
  double *objective_coeffs = (double*)malloc(model->num_vars * sizeof(double));

  for (int i = 0; i < model->num_vars; i++) {

    fscanf(textfile, "%lf", &objective_coeffs[i]);
  }

  // Skip the rest of the first line (semicolon and newline)
  char buffer[256];
  fgets(buffer, sizeof(buffer), textfile);

  pos = ftell(textfile);
  model->num_constraints = 0;
  int num_le = 0;  // <=
  int num_ge = 0;  // >=
  int num_eq = 0;  // =

  char line[256];
  while (fgets(line, sizeof(line), textfile)) {
    if (strchr(line, '<')) {
      num_le++;
      model->num_constraints++;
    } else if (strchr(line, '>')) {
      num_ge++;
      model->num_constraints++;
    } else if (strchr(line, '=')) {
      // Check if it's just = and not <= or >=
      if (!strchr(line, '<') && !strchr(line, '>')) {
        num_eq++;
        model->num_constraints++;
      }
    }
  }

  int num_slack_surplus = num_le + num_ge;  
  int num_artificial = num_eq + num_ge;      
  int total_cols = model->num_vars + num_slack_surplus + num_artificial;

  printf("Original variables: %d\n", model->num_vars);
  printf("<= constraints: %d (need %d slack vars)\n", num_le, num_le);
  printf(">= constraints: %d (need %d surplus + %d artificial)\n", num_ge, num_ge, num_ge);
  printf("=  constraints: %d (need %d artificial)\n", num_eq, num_eq);
  printf("Total matrix columns: %d\n", total_cols);
  printf("========================\n\n");

  model->columns = (double**)malloc(model->num_constraints * sizeof(double*));
  for (int i = 0; i < model->num_constraints; i++) {
    model->columns[i] = (double*)malloc(total_cols * sizeof(double));
    for (int j = 0; j < total_cols; j++) {
      model->columns[i][j] = 0.0;
    }
  }

  model->RHS_vector = (double*)malloc(model->num_constraints * sizeof(double));
  model->constraints_symbols = (char*)malloc(model->num_constraints * sizeof(char));

  // Allocate tracking vectors
  model->equalities_vector = (int*)malloc(num_artificial * sizeof(int));
  model->Basics_vector = (int*)malloc(model->num_constraints * sizeof(int));

  model->inequalities_count = num_slack_surplus;
  model->equalities_count = num_artificial;

  fseek(textfile, pos, SEEK_SET);

  // Temporary storage for constraint data
  typedef struct {
    double coeffs[256];
    double rhs;
    char type;  // 'L', 'G', 'E'
  } TempConstraint;

  TempConstraint *temp_constraints = (TempConstraint*)malloc(model->num_constraints * sizeof(TempConstraint));

  int constraint_idx = 0;
  while (fgets(line, sizeof(line), textfile) && constraint_idx < model->num_constraints) {
    char *token = strtok(line, " \t\n");
    int col_idx = 0;
    char inequality[3] = "";

    while (token != NULL) {
      if (strcmp(token, "<=") == 0 || strcmp(token, ">=") == 0 || strcmp(token, "=") == 0) {
        strcpy(inequality, token);
        token = strtok(NULL, " \t\n");
        break;
      }
      temp_constraints[constraint_idx].coeffs[col_idx++] = atof(token);
      token = strtok(NULL, " \t\n");
    }

    if (strcmp(inequality, "<=") == 0) {
      temp_constraints[constraint_idx].type = 'L';
    } else if (strcmp(inequality, "=") == 0) {
      temp_constraints[constraint_idx].type = 'E';
    } else if (strcmp(inequality, ">=") == 0) {
      temp_constraints[constraint_idx].type = 'G';
    }

    if (token != NULL) {
      temp_constraints[constraint_idx].rhs = atof(token);
    }

    constraint_idx++;
  }


  model->coeffs = (double*)malloc(sizeof(double) * (model->num_vars + model->equalities_count + model->inequalities_count));

  int slack_col = model->num_vars;
  int artificial_col = model->num_vars + num_slack_surplus;

  int slack_idx = 0;
  int artificial_idx = 0;

  for (int i = 0; i < model->num_constraints; i++) {
    // Copy original variable coefficients
    for (int j = 0; j < model->num_vars; j++) {
      model->columns[i][j] = temp_constraints[i].coeffs[j];
    }

    model->RHS_vector[i] = temp_constraints[i].rhs;
    model->constraints_symbols[i] = temp_constraints[i].type;

    if (temp_constraints[i].type == 'L') {
      // <= constraint: slack is basic
      int col = slack_col + slack_idx;
      model->columns[i][col] = 1.0;
      model->Basics_vector[i] = col; 
      model->coeffs[col] = 0.0;
      slack_idx++;

    } else if (temp_constraints[i].type == 'G') {
      // >= constraint: artificial is basic (not surplus!)
      int surplus_col = slack_col + slack_idx;
      int artif_col = artificial_col + artificial_idx;

      model->columns[i][surplus_col] = -1.0;
      model->columns[i][artif_col] = 1.0;

      model->equalities_vector[artificial_idx] = artif_col;
      model->Basics_vector[i] = artif_col;  

      model->coeffs[artif_col] = 9999.00; 
      slack_idx++;
      artificial_idx++;

    } else if (temp_constraints[i].type == 'E') {
      // = constraint: artificial is basic
      int artif_col = artificial_col + artificial_idx;
      model->columns[i][artif_col] = 1.0;
      model->equalities_vector[artificial_idx] = artif_col;
      model->Basics_vector[i] = artif_col;  
      model->coeffs[artif_col] = 9999.00;
      artificial_idx++;
    }
  }

  free(temp_constraints);

  printf("Number of slack/surplus variables: %d\n", model->inequalities_count);
  printf("Number of artificial variables: %d\n", model->equalities_count);

  if (model->num_constraints > model->num_vars) {

    printf("Warning: More constraints than variables! The Revised Simplex Algorithm works best when the number of variables exceeds the number of constraints\n");

  }


  model->solver_iterations = 1;

  for (size_t i = 0; i < model->num_vars; i++) {

    model->coeffs[i] = objective_coeffs[i];

  }

  model->Non_basics = (int*)malloc(sizeof(int) * total_cols);

  int non_basic_idx = 0;
  for (int col = 0; col < total_cols; col++) {
    int is_basic = 0;
    for (int j = 0; j < model->num_constraints; j++) {
      if (model->Basics_vector[j] == col) {
        is_basic = 1;
        break;
      }
    }
    if (!is_basic) {
      model->Non_basics[non_basic_idx++] = col;
    }
  }

  model->Non_basics_count = non_basic_idx;

  model->objective_function = 0.0;
  // Convert the problem to a maximization problem if it's a MIN problem, if it's a MIN problem then we convert to MAX

  if (strcmp(model->objective, "MINIMIZE") == 0) {
    for (size_t i = 0; i < model->num_vars; i++) {
      model->coeffs[i] *= -1; 
    }
    strcpy(model->objective, "MAXIMIZE");
  }else {
    for (size_t i = 0; i < model->num_vars; i++) {
      model->coeffs[i] *= -1; 
    }
    strcpy(model->objective, "MINIMIZE");

  }

  int artificial_start = model->num_vars + model->inequalities_count;
  double BIG_M = 9999.0;
  for (size_t i = 0; i < model->equalities_count; i++) {
    model->coeffs[artificial_start + i] = -BIG_M; 
  }


  return model;
}


void Print_columns(Model* model){
  int total_cols = model->num_vars + model->inequalities_count + model->equalities_count;

  printf("Constraint matrix:\n");
  for (int i = 0; i < model->num_constraints; i++) {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < total_cols; j++) {
      printf("%5.1f ", model->columns[i][j]);  
    }
    printf("| RHS: %5.1f", model->RHS_vector[i]);
    printf(" | Basic var: col %d\n", model->Basics_vector[i]);
  }
  printf("\n");

}


void invert_matrix(double **matrix, int n) {
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

double** Get_Basics_Matrix(Model* model) {

  double **Basics_matrix = (double**)malloc(model->num_constraints* sizeof(double));



  // printf("Basics matrix:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    Basics_matrix[i] = (double*)malloc(model->num_constraints * sizeof(double));
    for (size_t j = 0; j < model->num_constraints; j++) {
      int basics_index = model->Basics_vector[j];
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


void RevisedSimplex(Model* model){
  int termination = 0;
  int n = model->num_constraints;
  // printf("\n");
  //   printf("Non basics start:\n");
  //   for (size_t i = 0; i < model->Non_basics_count; i++) {
  //     printf(" %i ", model->Non_basics[i] + 1); 
  //   }
  //   printf("\n");
  //   printf("Basics vector:\n");
  //
  //   for (size_t i = 0; i < model->num_constraints; i++) {
  //     printf(" %i ", model->Basics_vector[i] + 1); 
  //   }
  //   printf("\n");
  while (termination != 1) {
    int feasibility_check = 0;
    printf("Beginning solver iteration %i ... \n", model->solver_iterations);
    double **B = Get_Basics_Matrix(model);
    double original_RHS[n];
    for (size_t i = 0; i < n; i++) {
      original_RHS[i] = model->RHS_vector[i]; 
    }

    invert_matrix(B, model->num_constraints);

    double *Simplex_multiplier = Get_simplex_multiplier(model, B);

    int entering_var_idx = 0;
    int entering_var = 0;
    update_RHS(model, original_RHS, B); 

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
    for (size_t i = 0; i < model->Non_basics_count; i++) {

      int non_basic_idx = model->Non_basics[i]; 


      double reduced_cost = Get_reduced_price(model, B, non_basic_idx, Simplex_multiplier ) ;

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
    if (feasibility_check == model->Non_basics_count) {
      printf("Optimal solution found, terminating!\n");
      termination++;
      Get_Objective_Function(model, original_RHS);
      break;
    }
    // printf("Best reduced cost is %f for the entering variable %i \n", best_reduced_cost, entering_var);
    double *Pivot = Get_pivot_column(B, model, entering_var);
    //update_RHS(model, original_RHS, B); 
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
      // printf("Ratio of variable %i is %f which has a pivot value of %f \n", model->Basics_vector[i], ratio, Pivot[i]);
      if (ratio < best_ratio && ratio > 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->Basics_vector[i];

      }

    }
    if (best_ratio < 0) {
      printf("LP is unbounded! Terminating!\n");
      exit(0); 
    }
    // printf("Best exiting variable is %i with ratio of %f at position %i \n", exiting_var, best_ratio, exiting_var_idx );
    model->Non_basics[entering_var_idx] = exiting_var;
    model->Basics_vector[exiting_var_idx] = entering_var;
    // model->Basics_vector[entering_var_idx] = entering_var;
    // model->Non_basics[exiting_var_idx] = exiting_var;
    // //
    // printf("Entering variable:%i at position %i \n", entering_var + 1, entering_var_idx);
    // printf("Exiting variable:%i at position: %i \n",exiting_var + 1,exiting_var_idx);
    // printf("Non basics updated:\n");
    // for (size_t i = 0; i < model->Non_basics_count; i++) {
    //   printf(" %i ", model->Non_basics[i] + 1); 
    // }
    // printf("\n");
    // printf("Basics vector:\n");
    //
    // for (size_t i = 0; i < model->num_constraints; i++) {
    //   printf(" %i ", model->Basics_vector[i] + 1); 
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

double Get_reduced_price(Model* model, double **B_inv, int var_col, double *multiplier_vector) {
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
double* Get_simplex_multiplier(Model* model, double **B_inv){
  int n = model->num_constraints;
  // THIS IS CORRECT AS Of THE VIDEO

  double *multiplier_vector = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      int basic_col_idx = model->Basics_vector[j];
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

void update_RHS(Model* model, double *RHS_vector_copy, double** B) {
  int n = model->num_constraints;
  double* temp = (double*)malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    temp[i] = 0.0;
    for (int j = 0; j < n; j++) {
      temp[i] += B[i][j] * RHS_vector_copy[j];
    
    }
  }


  for (int i = 0; i < n; i++) {
    RHS_vector_copy[i] = temp[i];
  }

  free(temp);
  //
  // printf("Updated RHS vector:\n");
  // for (int i = 0; i < n; i++) {
  //     printf("RHS[%d] = %f\n", i, RHS_vector[i]);
  // }
}
void Get_Objective_Function(Model* model, double *RHS_vector) {
  int n = model->num_constraints;
  for (int i = 0; i < n; i++) {
    int basic_idx = model->Basics_vector[i];
    if (model->coeffs[basic_idx] == 0 || RHS_vector[i] == 0  ) {
      continue;

    }
    double value = (model->coeffs[basic_idx] * RHS_vector[i]) * -1;
    model->objective_function += value;
    printf("Value of the variable x%i is %f with a coefficient of %f \n", basic_idx +1, RHS_vector[i], model->coeffs[basic_idx] * -1);

  }

}
void free_model(Model* model){
  for (int i = 0; i < model->num_constraints; i++) {
    free(model->columns[i]);
  }
  free(model->columns);

  free(model->coeffs);
  free(model->RHS_vector);
  free(model->Basics_vector);
  free(model->constraints_symbols);
  free(model->equalities_vector);
  free(model->Non_basics);

  printf("-------------------------------\n");
  printf("Model free'd from the heap\n"); 
}

#endif 


