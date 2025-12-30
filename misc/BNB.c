#include <string.h>
#include "BNB.h"
#include <stdio.h>
#include <stdlib.h>
#include "RSA.h"
#include <math.h>



void BNB(Model *model, Subproblem *subproblem) {

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
      SubproblemCheck(model, original_RHS, subproblem);

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

void SubproblemCheck(Model *model, double *rhs_vector, Subproblem* subproblem)
{
  size_t n = model->num_constraints;



  for (int i = 0; i < n; i++)
  {
    int basic_idx = model->basics_vector[i];

    // Infeasibility check!
    if (fabs(model->coeffs[basic_idx]) == model->BIG_M  && rhs_vector[i] > 1e-6)
    {
      printf("Model Infeasible. Artificial basic variable has a positive RHS value. Terminating!\n");
      subproblem->status = INFEASIBLE;  
      return;
    }

    model->objective_function += model->coeffs[basic_idx] * rhs_vector[i];

    if (fabs(model->coeffs[basic_idx]) > 1e-6 &&
      fabs(model->coeffs[basic_idx]) < model->BIG_M &&
      rhs_vector[i] > 1e-6)
    {

      
      // I am aware of 4 branching techniques, I guess I will branch whichever fractional variable I get first  
      if (fabs(rhs_vector[i] - fabs(rhs_vector[i])) > 1e-9){
        subproblem->variable = basic_idx;
        subproblem->status = BRANCH; 
        return;
      } 

      

     }

  }
  // if the solution is feasible and we have nothing to branch -> solution is integer 
  subproblem->status = INTEGER;
}


