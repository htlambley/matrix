#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
    
int main(void) {
    matrix_t* A = matrix_id(2);
    matrix_set(A, 0, 1, 5.0);
    matrix_t* b = matrix_create(2, 1);
    matrix_set(b, 0, 0, 1.0);
    matrix_t* x = matrix_gepp(A, b);
    printf("x = \n\n");
    matrix_print(x);
    matrix_destroy(A);
    matrix_destroy(x);
    matrix_destroy(b);
}
