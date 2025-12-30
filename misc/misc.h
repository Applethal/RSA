#ifdef MISC_H 
#define MISC_H 
#include "RSA.h"


void AddConstraint(Model* model, double* lhs, size_t lhs_count, char symbol, double rhs);
void AddVar(Model* model, double var_coeff);



#endif // MISC_H 

