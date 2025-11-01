#include <stdio.h>
#include <stdlib.h>
#include "RSA.h"

int main(int argc, char *argv[]) {
  
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

 
  
  PrintColumns(model);

  printf("Starting solving\n"); 

  printf("\n");

  

 if (Debug == 1) {

    printf("Debug Mode: On\n");
    RevisedSimplex_Debug(model);

 } else {
    RevisedSimplex(model); 

  } 

 


  
  printf("\n");
  printf("Objective function: %f \n", model->objective_function);


  fclose(file);
  FreeModel(model); 
  return 0;
}
