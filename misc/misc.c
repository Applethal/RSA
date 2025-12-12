#include "misc.h"
#include <stdio.h>
#include <stdlib.h>


void AddConstraint(Model* model, double* lhs, size_t lhs_count, char symbol, double rhs) 
{
     if (lhs_count > model->num_vars) {
        fprintf(stderr,
                "Error: LHS has %zu coefficients, but the model only has %zu variables.\n",
                lhs_count, model->num_vars);
        exit(1);
    }

    if (symbol != 'L' && symbol != 'R' && symbol != 'E') {
        fprintf(stderr,
                "Error: Invalid constraint symbol '%c'. Allowed: L, R, E.\n",
                symbol);
        exit(1);
    }

    size_t newCount = model->num_constraints + 1;

 
    double **tmpColumns = realloc(model->columns,
                                  newCount * sizeof(double*));

    if (tmpColumns == NULL) {
        fprintf(stderr, "Error: realloc failed for model->columns\n");
        exit(1);
    }
    model->columns = tmpColumns;

    model->columns[model->num_constraints] =
        malloc(model->num_vars * sizeof(double));
    if (model->columns[model->num_constraints] == NULL) {
        fprintf(stderr, "Error: malloc failed for new row in columns\n");
        exit(1);
    }

    for (size_t j = 0; j < lhs_count; j++) {
        model->columns[model->num_constraints][j] = lhs[j];
    }

    for (size_t j = lhs_count; j < model->num_vars; j++) {
        model->columns[model->num_constraints][j] = 0.0;
    }

    double *tmpRhs = realloc(model->rhs_vector,
                             newCount * sizeof(double));
    if (tmpRhs == NULL) {
        fprintf(stderr, "Error: realloc failed for rhs_vector\n");
        exit(1);
    }
    model->rhs_vector = tmpRhs;

    model->rhs_vector[model->num_constraints] = rhs;

  
    char *tmpSymbols = realloc(model->constraints_symbols,
                               newCount * sizeof(char));
    if (tmpSymbols == NULL) {
        fprintf(stderr, "Error: realloc failed for constraints_symbols\n");
        exit(1);
    }
    model->constraints_symbols = tmpSymbols;

    model->constraints_symbols[model->num_constraints] = symbol;

    model->num_constraints++;
}

void AddVar(Model* model, double var_coeff)
{
    size_t newVarCount = model->num_vars + 1;

    double *tmpCoeffs =
        realloc(model->coeffs, sizeof(double) * newVarCount);

    if (tmpCoeffs == NULL) {
        fprintf(stderr, "Error: realloc failed for coeffs\n");
        exit(1);
    }

    model->coeffs = tmpCoeffs;
    model->coeffs[model->num_vars] = var_coeff;

    for (size_t i = 0; i < model->num_constraints; i++) {

        double *tmpRow =
            realloc(model->columns[i], sizeof(double) * newVarCount);

        if (tmpRow == NULL) {
            fprintf(stderr,
                    "Error: realloc failed for columns row %zu\n", i);
            exit(1);
        }

        model->columns[i] = tmpRow;

        model->columns[i][model->num_vars] = 0.0;
    }


    model->num_vars++;
}
