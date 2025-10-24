#ifndef RSA_H
#define RSA_H
#include <string.h>
#include <math.h>
// #include <float.h> Do I need this to set the Big Ms? lol
typedef struct {
  char objective[10];   // MINIMIZE or MAXIMIZE
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

Model* read(FILE *textfile);

void Get_Objective_Function(Model* model, double **B); 

void free_model(Model *model);

void invert_matrix(double **matrix, int n);

void RevisedSimplex(Model* model);

double** Get_Basics_Matrix(Model* model);

void Print_columns(Model* model);

double Get_reduced_price(Model* model, double **B_inv, int var_col, double *multiplier_vector);

double* Get_simplex_multiplier(Model* model, double **B_inv);

double* Get_pivot_column(double** B_inv, Model* model, int best_cost_idx);



Model* read(FILE *textfile) {
  Model *model = (Model*)malloc(sizeof(Model));

  fscanf(textfile, "%s", model->objective);

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
  printf("Total columns in tableau: %d\n", total_cols);
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


  model->solver_iterations = 0;

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

  for (size_t i = 0; i < model->num_vars + model->inequalities_count + model->equalities_count; i++) {
    model->coeffs[i] *= -1; 
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

  printf("Basic variables (column indices): ");
  for (int i = 0; i < model->num_constraints; i++) {
    printf("%d ", model->Basics_vector[i]);
  }
  printf("\n\n");
}


void invert_matrix(double **matrix, int n) {
  // Check if matrix is square (you'll need to pass n as parameter)
  // Create augmented matrix [A|I]
  double **augmented = (double**)malloc(n * sizeof(double*));
  for (int i = 0; i < n; i++) {
    augmented[i] = (double*)malloc(2 * n * sizeof(double));
    // Copy original matrix
    for (int j = 0; j < n; j++) {
      augmented[i][j] = matrix[i][j];
    }
    // Add identity matrix
    for (int j = n; j < 2 * n; j++) {
      augmented[i][j] = (j - n == i) ? 1.0 : 0.0;
    }
  }

  // Gauss-Jordan elimination
  for (int i = 0; i < n; i++) {
    // Find pivot
    double pivot = augmented[i][i];

    if (fabs(pivot) < 1e-10) {
      printf("WARNING: Degenerate matrix!\n");
      // Free memory
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

  // Copy inverse back to original matrix
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      matrix[i][j] = augmented[i][n + j];
    }
  }

  // Free augmented matrix
  for (int i = 0; i < n; i++) {
    free(augmented[i]);
  }
  free(augmented);
}

double** Get_Basics_Matrix(Model* model) {

  double **Basics_matrix = (double**)malloc(model->num_constraints* sizeof(double));




  for (size_t i = 0; i < model->num_constraints; i++) {
    Basics_matrix[i] = (double*)malloc(model->num_constraints * sizeof(double));
    for (size_t j = 0; j < model->num_constraints; j++) {
      int basics_index = model->Basics_vector[j];
      Basics_matrix[i][j] = model->columns[i][basics_index]; 
    } 
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

  // Save original RHS values - ESSENTIAL!
  double *original_RHS = malloc(model->num_constraints * sizeof(double));
  for (size_t i = 0; i < model->num_constraints; i++) {
    original_RHS[i] = model->RHS_vector[i];
  }

  // Determine actual solving direction
  int is_maximization = (strcmp(model->objective, "MAXIMIZE") == 0) || 
                        (strcmp(model->objective, "MINIMIZE") == 0);  // Both become MAX after negation

  while (termination != 1) {
    printf("Beginning solver iteration %i \n", model->solver_iterations);
    double **B = Get_Basics_Matrix(model);
    invert_matrix(B, model->num_constraints);

    // UPDATE RHS: x_B = B^(-1) * b (original RHS)
    for (size_t i = 0; i < model->num_constraints; i++) {
      model->RHS_vector[i] = 0.0;
      for (size_t j = 0; j < model->num_constraints; j++) {
        model->RHS_vector[i] += B[i][j] * original_RHS[j];
      }
    }

    double *Simplex_multiplier = Get_simplex_multiplier(model, B);
    int entering_var_idx = -1;
    int entering_var = -1;
    double best_reduced_cost = 0.0;
    
    printf("Problem solving as MAXIMIZATION (coefficients already negated if MIN)\n");
    
    // Getting the entering variable - ALWAYS use MAX logic since coeffs are negated
    for (size_t i = 0; i < model->Non_basics_count; i++) {
      int non_basic_idx = model->Non_basics[i]; 
      
      // Skip if non_basic_idx is out of bounds
      if (non_basic_idx >= model->num_vars) {
        continue;
      }
      
      double reduced_cost = Get_reduced_price(model, B, non_basic_idx, Simplex_multiplier);
      printf("Variable %i reduced cost: %f\n", non_basic_idx, reduced_cost);
      
      // For MAXIMIZATION: choose most negative reduced cost
      if (reduced_cost < best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
    }
    
    // Check for optimality - MAX logic (all reduced costs >= 0)
    if (best_reduced_cost >= 0) {
      printf("Optimal solution found at iteration %i! All reduced costs >= 0\n", model->solver_iterations);
      Get_Objective_Function(model, B);
      termination = 1;
      free(Simplex_multiplier);
      for (int i = 0; i < model->num_constraints; i++) {
        free(B[i]);
      }
      free(B);
      break;
    }
    
    // If no entering variable found
    if (entering_var_idx == -1) {
      printf("No entering variable found - optimal solution reached\n");
      Get_Objective_Function(model, B);
      termination = 1;
      free(Simplex_multiplier);
      for (int i = 0; i < model->num_constraints; i++) {
        free(B[i]);
      }
      free(B);
      break;
    }
    
    printf("Best reduced cost is %f for entering variable %i\n", best_reduced_cost, entering_var);
    double *Pivot = Get_pivot_column(B, model, entering_var);
    
    // Getting the exiting variable   
    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = 1e20;
    
    int positive_pivot_found = 0;
    for (size_t i = 0; i < model->num_constraints; i++) {
      if (Pivot[i] > 1e-10) {
        positive_pivot_found = 1;
        double ratio = model->RHS_vector[i] / Pivot[i];
        if (ratio >= 0 && ratio < best_ratio) {
          best_ratio = ratio;
          exiting_var_idx = i;
          exiting_var = model->Basics_vector[i];
        }
      }
    }
    
    // Check for unboundedness
    if (!positive_pivot_found || exiting_var_idx == -1) {
      printf("Problem is unbounded! No positive pivot element found.\n");
      termination = 1;
      free(Pivot);
      free(Simplex_multiplier);
      for (int i = 0; i < model->num_constraints; i++) {
        free(B[i]);
      }
      free(B);
      break;
    }
    
    printf("Best exiting variable is %i with ratio of %f\n", exiting_var, best_ratio);
    
    // Update basis and non-basics
    model->Basics_vector[exiting_var_idx] = entering_var;
    model->Non_basics[entering_var_idx] = exiting_var;
    
    printf("Entering variable: %i at position %i\n", entering_var, entering_var_idx);
    printf("Exiting variable: %i at position: %i\n", exiting_var, exiting_var_idx);
    printf("Non basics updated: ");
    for (size_t i = 0; i < model->Non_basics_count; i++) {
      printf("%i ", model->Non_basics[i]); 
    }
    printf("\nBasics vector: ");
    for (size_t i = 0; i < model->num_constraints; i++) {
      printf("%i ", model->Basics_vector[i]); 
    }
    printf("\n\n");
    
    model->solver_iterations++;
    
    // Free allocated memory
    for (int i = 0; i < model->num_constraints; i++) {
      free(B[i]);
    }
    free(B);
    free(Simplex_multiplier);
    free(Pivot);
  }
  
  free(original_RHS);
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
        // printf("Pivot[%d] = %f\n", i, Pivot[i]);
    }

    return Pivot;
}
void Get_Objective_Function(Model* model, double **B) {
    double objective_value = 0.0;
    
    // Compute x_B = B^(-1) * b
    // Since B is already inverted, x_B is in model->RHS_vector
    // Or you need to multiply B_inv with original RHS
    if (model->num_vars < model->num_constraints) {
      for (size_t i = 0; i < model->num_vars; i++) {
        int basic_var = model->Basics_vector[i];
        double basic_var_value = model->RHS_vector[i];
        double coefficient = model->coeffs[basic_var] * -1;
        printf("Final value for the x%i variable: %f \n", basic_var, basic_var_value); 
        model->objective_function += coefficient * basic_var_value;
    }

      
    }
    if (model->num_constraints < model->num_vars) {
      for (size_t i = 0; i < model->num_constraints; i++) {
        int basic_var = model->Basics_vector[i];
        double basic_var_value = model->RHS_vector[i];
        double coefficient = model->coeffs[basic_var] * -1;
        printf("Final value for the x%i variable: %f \n", basic_var, basic_var_value); 
        model->objective_function += coefficient * basic_var_value;
    }

      
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
  free(model);
  printf("-------------------------------\n");
  printf("Model free'd from the heap\n"); 
}

#endif 


