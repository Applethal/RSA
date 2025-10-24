#include <stdio.h>
#include <stdlib.h>
#include "RSA.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("Insert the file, exiting\n");
    exit(0);
  }

  FILE *file = fopen(argv[1], "r");
  if (file == NULL) {
    printf("Cannot open file\n");
    exit(1);
  }

  Model *model = read(file);

  printf("Objective: %s\n", model->objective);
  printf("Number of variables: %d\n", model->num_vars);
  printf("Number of constraints: %d\n", model->num_constraints);

  printf("Objective coefficients:");
  for (int i = 0; i < model->num_vars + model->equalities_count + model->inequalities_count; i++) {
    printf(" %.1f", model->coeffs[i]);  
  }
  printf("\n\n");


  printf("Constraint matrix (With slack variables):\n");

  Print_columns(model);

  printf("Starting solving\n"); 


  printf("\n");

  printf("Symbols:");

  for (int i = 0;i < model->num_constraints;i++) {
    printf(" %c ", model->constraints_symbols[i]);
  }
  printf("\n");
  
  printf("Non basic variables:\n");
  for (int i = 0;i < model->Non_basics_count ; i++) {
    printf(" %i ", model->Non_basics[i]);

  }

  printf("\n");

  RevisedSimplex(model); 
  

  printf("Objective function: %f \n", model->objective_function);


  fclose(file);
  free_model(model); 
  return 0;
}
