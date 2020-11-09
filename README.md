# Matrix computation library
A C library for basic matrix algebra including:
- Forward and backward substitution;
- **Decompositions**: LU (including with pivoting), QR;
- **Solvers for systems of linear equations**: Gaussian elimination, Jacobi method;
- **Solvers for least squares problems**: using QR factorisation.

## Compilation
Compile using `make`, and if desired, install the library to your system (root may be required):
```shell
make package
cp libmatrix.a /usr/local/lib/libmatrix.a
cp matrix.h /usr/local/include/matrix.h
```
In order to compile, `gcc` and `make` are required. These should be easily obtained using your 
computer's package manager.

If the library has been installed as above, compile your program with `-lmatrix` in order to
use the library.

## Usage
The following example demonstrates how to use the library to solve a system Ax = b using the 
Gaussian elimination algorithm. This is suitable for any invertible matrix A. 
```c
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

