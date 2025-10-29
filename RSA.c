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
  // char  *objective= argv[2];

  Model *model = read_csv(file);

 
  
  Print_columns(model);

  printf("Starting solving\n"); 

  printf("\n");

 
  RevisedSimplex(model); 
  
  printf("\n");
  printf("Objective function: %f \n", model->objective_function);


  fclose(file);
  free_model(model); 
  return 0;
}
