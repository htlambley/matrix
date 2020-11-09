# Matrix computation library
A C library for basic matrix algebra including:
- Forward and backward substitution;
- **Decompositions**: LU (including with pivoting), QR;
- **Solvers for systems of linear equations**: Gaussian elimination, Jacobi method;
- **Solvers for least squares problems**: using QR factorisation.

## Compilation
Compile using `make`:
```
make all
```
In order to compile, `gcc` and `make` are required. These should be easily obtained using your 
computer's package manager.

## Usage
The following example demonstrates how to use the library to solve a system Ax = b using the 
Gaussian elimination algorithm. This is suitable for any invertible matrix A. 
```
#include "matrix.h"
int main(void) {
    // Create a 3 x 3 identity matrix.
    matrix_t* A = matrix_id(3);
    // Set the element in the second row, first column to 2.
    // Note that indexing of matrices begins at zero.
    matrix_set(A, 1, 0, 2.0);
    // Create a column vector b = (1, 2, 3)^T.
    matrix_t* b = matrix_create(3, 1);
    matrix_set(b, 0, 0, 1.0);
    matrix_set(b, 1, 0, 2.0);
    matrix_set(b, 2, 0, 3.0);
    // Use Gaussian elimination to solve the system Ax = b for the vector x.
    matrix_t* x = matrix_gepp(A, b);
    printf("Solution:\n");
    matrix_print(x);
    // Destroy matrices after use.
    matrix_destroy(A);
    matrix_destroy(x);
    matrix_destroy(b);
}
```
