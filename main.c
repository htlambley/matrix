#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main(void) {
    matrix_t* id = matrix_id(2);
    double a_11 = matrix_get(id, 1, 1);
    printf("%f\n", a_11);
    return 0;
}

