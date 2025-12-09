#include <stdio.h>
#include <stdlib.h>
#include "RSA.h"
#include <time.h>
int main(int argc, char *argv[]) {
  clock_t begin = clock(); 
  int Debug = 0;

  if (argc >= 3 && strcmp(argv[2], "-Debug") == 0) {
    Debug = 1;
  }




  if (argc < 2) {
    printf("Insert the file, exiting\n");
    exit(0);
  }

  FILE *file = fopen(argv[1], "r");
  if (file == NULL) {
    printf("Cannot open file\n");
    exit(1);
  }
  // char  *objective= argv[2];

  Model *model = ReadCsv(file);





  printf("Starting solver, to enable iterative debugging add the flag '-Debug' as an argument after the file path. \n"); 

  printf("\n");



  if (Debug == 1) {

    printf("Debug Mode: On\n");
    printf("Objective function mode: %s\n", model->objective);
    printf("Number of variables: %d\n", model->num_vars);
    printf("Number of constraints: %d\n", model->num_constraints);
    printf("Objective coefficients (Slack and Artificial coeffs included):");
    for (int i = 0; i < model->num_vars + model->equalities_count + model->inequalities_count; i++) {
      printf(" %.1f", model->coeffs[i]);  
    }
    printf("\n\n");
    PrintColumns(model);

    RevisedSimplex_Debug(model);
    
  } else {
    RevisedSimplex(model); 

  } 





  printf("\n");
  printf("Objective function: %f \n", model->objective_function);
  clock_t end = clock();
  
  printf("Solving time: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);

  fclose(file);
  FreeModel(model);

  return 0;
}
