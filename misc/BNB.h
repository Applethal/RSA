#ifndef BNB_H 
#define BNB_H 
#include <RSA.h>


#define INFEASIBLE 2
#define INTEGER 1
#define BRANCH 0 
#define VISITED 3


typedef struct Subproblem{

  int id;
  int variable; // contains the idx of the variable to be fixed
  int status; // 1 = integer solution, 0 = branch, 2 = infeasible, 3 = already visited

} Subproblem;


void BNB(Model* model, Subproblem* subproblem);


void SubproblemCheck(Model* model,double* original_RHS,Subproblem* subproblem);




#endif // !BNB_H 
