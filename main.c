#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int main(void) {
    matrix_t* A = matrix_id(250);
    for (int i = 0; i < 100; i++) {
	matrix_t* U = matrix_clone(A);
	matrix_t* L = matrix_lu(U);	
	matrix_destroy(L);
	matrix_destroy(U);
    }
    return 0;
}

