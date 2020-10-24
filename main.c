#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
int main(void) {
    matrix_t* A = matrix_id(2);
    matrix_set(A, 0, 1, 0.5);
    matrix_set(A, 1, 0, -0.1);
    printf("A = \n\n");
    matrix_print(A);
    matrix_t* b = matrix_create(2, 1);
    matrix_set(b, 1, 0, 1.0);
    
    matrix_t* x = matrix_create(2, 1);
    matrix_set(x, 0, 0, 1.0);
    matrix_set(x, 1, 0, 1.0);
    double err = matrix_jacobi(A, b, &x, 0.01);
    printf("Ax = (%llu x %llu) x (%llu x %llu) \n\n", A->m, A->n, x->m, x->n);
     
    matrix_t* B = matrix_mul(A, x);
    matrix_print(B);
}

