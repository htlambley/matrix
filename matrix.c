#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

#define EPSILON 1e-6

bool approx_equal(double x, double y, double eps) {
    return fabs(x - y) < eps;
}

bool matrix_upper_triangular(matrix_t* A) {
    if (A->m != A->n) {
	return false;
    }

    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < i; j++) {
	    if (!approx_equal(matrix_get(A, i, j), 0, EPSILON)) {
		return false;
	    }
	}	
    }
    return true;
}

bool matrix_lower_triangular(matrix_t* A) {
    if (A->m != A->n) {
	return false;
    }

    for (uint64_t i = 0; i < A->m - 1; i++) {
	for (uint64_t j = i + 1; j < A->n; j++) {
	    if (!approx_equal(matrix_get(A, i, j), 0, EPSILON)) {
		return false;
	    }
	}
    }
    return true;
}

bool matrix_diagonal(matrix_t* A) {
    if (A->m != A->n) {
	return false;
    }
    
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    if (i == j) {
		continue;
	    } else {
		if (!approx_equal(matrix_get(A, i, j), 0, EPSILON)) {
		    return false;
		}
	    }
	}
    }
    return true;
}

/*
 * Performs forward substitution on a lower triangular matrix A to
 * solve the equation Ax = b for x.
 * Takes an n x n lower triangular matrix A and a n x 1 vector b.
 * You must check that A is lower triangular yourself: for example,
 * use matrix_lower_triangular(A).
 */
matrix_t* matrix_fs(matrix_t* A, matrix_t* b) {
    matrix_t* x = matrix_create(A->n, 1);
    for (uint64_t i = 0; i < A->n; i++) {
	double b_i = matrix_get(b, i, 0);
	double sum = 0.;
	for (uint64_t j = 0; j < i; j++) {
	    sum += matrix_get(A, i, j) * matrix_get(x, j, 0);
	}
	double a_ii = matrix_get(A, i, i);
	matrix_set(x, i, 0, (b_i - sum) / a_ii);
    }
    return x;
}

/* 
 * Performs backward substitution on an upper triangular matrix A. 
 * See matrix_fs.
 */
matrix_t* matrix_bs(matrix_t* A, matrix_t* b) {
    matrix_t* x = matrix_create(A->n, 1);
    for (uint64_t i = A->n; i > 0; i--) {
	double i_ = i - 1;
	double b_i = matrix_get(b, i_, 0);
	double sum = 0.;
	for (uint64_t j = A->n - 1; j > i_; j--) {
	    sum += matrix_get(A, i_, j) * matrix_get(x, j, 0);
	}
	double a_ii = matrix_get(A, i_, i_);
	matrix_set(x, i_, 0, (b_i - sum) / a_ii);
    }
    return x;
}

/* Solves a linear system Ax = b for any square matrix A
 * by Gaussian elimination (without pivoting).
 * 
 * Takes an n x n matrix A and a n x 1 vector b.
 * Consumes the matrix A (so do not use A after this).
 * Returns an n x 1 vector x satisfiying Ax = b.
 */
matrix_t* matrix_ge(matrix_t* A, matrix_t* b) {
    matrix_t* L = matrix_lu(A); 
    // Solve Ly = b.
    matrix_t* y = matrix_fs(L, b);
    // Solve Ux = y.
    matrix_t* x = matrix_bs(A, y);
    return x;
}

matrix_t* matrix_gepp(matrix_t* A, matrix_t* b) {
    matrix_lup_t* decomp = matrix_lupp(A);
    // Solve PAx = Pb by evaluating
    // LUx = Pb. So, evaluate:
    // Ly = Pb by FS and Ux = y by BS.
    matrix_t* Pb = matrix_mul(decomp->P, b);
    matrix_t* y = matrix_fs(decomp->L, Pb);
    matrix_t* x = matrix_bs(A, y);
    return x;
}

/*
 * Allocates an m x n matrix and returns a pointer to it.
 * Each element is guaranteed to be zero.
 */
matrix_t* matrix_create(uint64_t m, uint64_t n) {
    double* A = calloc(m * n, sizeof(double));
    struct matrix *mat = (struct matrix*)malloc(sizeof(struct matrix));
    mat->A = A;
    mat->m = m;
    mat->n = n;
    return mat;
}

void matrix_destroy(matrix_t* A) {
    free(A->A);
    free(A);
}

__attribute__((always_inline))
inline double matrix_get(struct matrix* m, uint64_t i, uint64_t j) {
    return m->A[i * m->n + j]; 
}

__attribute__((always_inline))
inline void matrix_set(struct matrix* m, uint64_t i, uint64_t j, double value) {
    m->A[i * m->n + j] = value;
}

/*
 * Adds two matrices, storing the value in the first matrix.
 */
void matrix_add(struct matrix* A, struct matrix* B) {
    assert(A->m == B->m && A->n == B->n);
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    double result = matrix_get(A, i, j) + matrix_get(B, i, j);
	    matrix_set(A, i, j, result);
	}
    }
}

// Creates the matrix I_n.
struct matrix* matrix_id(uint64_t n) {
    struct matrix* A = matrix_create(n, n);
    for (uint64_t i = 0; i < A->m; i++) {
	matrix_set(A, i, i, 1);
    }
    return A;
}

/* 
 * Computes the matrix product AB, storing the result in a new matrix of 
 * the appropriate dimension.
 */
matrix_t* matrix_mul(matrix_t* A, matrix_t* B) {
    assert(A->n == B->m);
    matrix_t* C = matrix_create(A->m, B->n);
    for (uint64_t i = 0; i < C->m; i++) {
	for (uint64_t j = 0; j < C->n; j++) {
	    double element = 0.;
	    for (uint64_t k = 0; k < A->n; k++) {
		element += matrix_get(A, i, k) * matrix_get(B, k, j);
	    }
	    matrix_set(C, i, j, element);
	}
    }
    return C;
}

/*
 * Computes the LU decomposition of a matrix. An LU decomposition is a 
 * factorisation of a matrix A into a unit lower triangular matrix L
 * and an upper triangular matrix U such that A = LU.
 *
 * This function takes the given matrix and returns the matrix L. The given
 * matrix will be altered in place to give U.
 *
 * The given matrix must be square.
 */
matrix_t* matrix_lu(matrix_t* U) {
    assert(U->m == U->n);
    matrix_t* L = matrix_id(U->m);
    for (uint64_t j = 0; j < U->n - 1; j++) {
	for (uint64_t i = j + 1; i < U->n; i++) {
	    double l_ij = matrix_get(U, i, j) / matrix_get(U, j, j);
	    matrix_set(L, i, j, l_ij);
	    matrix_set(U, i, j, 0.0);
	    for (uint64_t k = j + 1; k < U->n; k++) {
		double val = matrix_get(U, i, k) - l_ij * matrix_get(U, j, k);
		matrix_set(U, i, k, val);
	    }
	}
    }
    return L;
}

/* Swaps the rows k and l, from column a (inclusive) to column b (exclusive).
 * Example: matrix_row_swap_partial(A, 0, 1, 1, A->n) swaps all but the first column.
 */
void matrix_row_swap_partial(matrix_t* A, uint64_t k, uint64_t l, uint64_t a, uint64_t b) {
    for (uint64_t j = a; j < b; j++) {
	double tmp = matrix_get(A, k, j);
	matrix_set(A, k, j, matrix_get(A, l, j));
	matrix_set(A, l, j, tmp);
    }
}

void matrix_row_swap(matrix_t* A, uint64_t k, uint64_t l) {
    matrix_row_swap_partial(A, k, l, 0, A->n);
}

matrix_lup_t* matrix_lupp(matrix_t* U) {
    assert(U->m == U->n);
    matrix_t* L = matrix_id(U->m);
    matrix_t* P = matrix_id(U->m);
    for (uint64_t j = 0; j < U->n - 1; j++) {
	// Find the largest absolute value in the column.
	uint64_t max_i = 0;
	double max_abs = 0.0;
	for (uint64_t i = 0; i < U->n; i++) {
	    double u_ij = fabs(matrix_get(U, i, j));
	    if (u_ij > max_abs) {
		max_abs = u_ij;
		max_i = i;	
	    } 
	}

	// Swap (u_ij, ..., u_in) with (u_jj, ..., u_jn),
	//      (l_i1, ..., l_i,j-1) with (l_j1, ..., l_j,j-1)
	// Swap row i with k in P.
	matrix_row_swap_partial(U, max_i, j, j, U->n);
	matrix_row_swap_partial(L, max_i, j, 0, j);
	matrix_row_swap(P, max_i, j);
	for (uint64_t i = j + 1; i < U->m; i++) {
	    double l_ij = matrix_get(U, i, j) / matrix_get(U, j, j);
	    matrix_set(L, i, j, l_ij);
	    matrix_set(U, i, j, 0.0);
	    for (uint64_t k = j + 1; k < U->n; k++) {
		double u_ik = matrix_get(U, i, k);
		u_ik -= l_ij * matrix_get(U, j, k);
		matrix_set(U, i, k, u_ik);
	    }
	}
    } 
    matrix_lup_t* ret = (matrix_lup_t*)malloc(sizeof(matrix_lup_t));
    ret->L = L;
    ret->P = P;
    return ret;
}


void matrix_print(matrix_t* A) {
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    printf("%+f   ", matrix_get(A, i, j));
	}
	printf("\n\n");
    }
}

/* Computes the Hadamard product (the elementwise product 
 * of two matrices) of A and B.
 * The result is stored in A.
 */
void matrix_hadamard(matrix_t* A, matrix_t* B) {
    assert(A->m == B->m && A->n == B->n);
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    double a_ij = matrix_get(A, i, j);
	    double b_ij = matrix_get(B, i, j);
	    matrix_set(A, i, j, a_ij * b_ij);
	}
    }
}

double matrix_dot(matrix_t* A, matrix_t* B) {
    assert(A->n == 1 && B->n == 1 && A->m == B->m);
    double sum = 0.0;
    for (uint64_t i = 0; i < A->m; i++) {
	sum += matrix_get(A, i, 0) * matrix_get(B, i, 0);	
    }
    return sum;
}

double matrix_norm_frob(matrix_t* A) {
    double norm = 0.0;
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    double a_ij = matrix_get(A, i, j);
	    norm += a_ij * a_ij;
	}
    }
    return sqrt(norm);
}

void matrix_scalar_mul(matrix_t* A, double k) {
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    double a_ij = matrix_get(A, i, j);
	    matrix_set(A, i, j, k * a_ij);
	}
    }
}

matrix_t* matrix_get_col(matrix_t* A, uint64_t j) {
    matrix_t* col = matrix_create(A->m, 1);
    for (uint64_t i = 0; i < A->m; i++) {
	matrix_set(col, i, 0, matrix_get(A, i, j));
    }
    return col;
}

void matrix_sub_smul(matrix_t* A, matrix_t* B, double k) {
    assert(A->m == B->m && A->n == B->n);
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    matrix_set(A, i, j, matrix_get(A, i, j) - k * matrix_get(B, i, j));
	}
    }
}

void matrix_set_col(matrix_t* A, uint64_t j, matrix_t* col) {
    assert(A->m == col->m); 
    for (uint64_t i = 0; i < A->m; i++) {
	matrix_set(A, i, j, matrix_get(col, i, 0));
    } 
}

/* Computes the reduced QR factorisation of an m x n
 * matrix A. 
 * This is a decomposition A = QR, where Q is an m x n
 * orthogonal matrix and R is an upper triangular n x n
 * matrix.
 * We require m >= n for this factorisation to exist.
 * Returns a matrix_qr_t containing members Q and R
 * which are newly allocated matrices of the two
 * required forms.
 */
matrix_qr_t* matrix_reduced_qr(matrix_t* A) {
    assert(A->m >= A->n);
    matrix_t* Q = matrix_create(A->m, A->n);
    matrix_t* R = matrix_create(A->n, A->n);
    matrix_t* q[A->n];
    matrix_t* a[A->n];
    for (uint64_t j = 0; j < A->n; j++) {
	q[j] = matrix_get_col(A, j);
	a[j] = matrix_get_col(A, j);

	// Remove the projections so that
	// q_j = a_j - \sum_{i = 1}^{n - 1} <q_i, a_j> q_i.
	for (uint64_t i = 0; i < j; i++) {
	    double dot = matrix_dot(q[i], q[j]);
	    matrix_sub_smul(q[j], q[i], dot);	    
	}

	// Normalise u_k.
	double norm = matrix_norm_frob(q[j]);
	matrix_scalar_mul(q[j], 1.0/norm);
	matrix_set_col(Q, j, q[j]);
    }
    
    for (uint64_t j = 0; j < A->n; j++) {
	// Fill out the R matrix.
	for (uint64_t i = 0; i < A->n; i++) {
	    double dot = matrix_dot(q[i], a[j]);
	    matrix_set(R, i, j, dot);
	}
    }
    for (uint64_t j = 0; j < A->n; j++) {
	matrix_destroy(q[j]);
	matrix_destroy(a[j]);
    }
    matrix_qr_t* ret = (matrix_qr_t*)malloc(sizeof(matrix_qr_t));
    ret->Q = Q;
    ret->R = R;
    return ret;
}

matrix_t* matrix_transpose(matrix_t* A) {
    matrix_t* A_T = matrix_create(A->n, A->m);
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    matrix_set(A_T, j, i, matrix_get(A, i, j));
	}
    }
    return A_T;
}

/* Solves the least squares (LSQ) problem 
 * to minimise the function
 * g(x) = 1/2 |Ax - b|^2.
 * Takes an m x n matrix A and vector b in R^m
 * and returns a vector x in R^n.
 */
matrix_t* matrix_lsq(matrix_t* A, matrix_t* b) {
    matrix_qr_t* decomp = matrix_reduced_qr(A);
    matrix_t* Q_T = matrix_transpose(decomp->Q);
    matrix_t* y = matrix_mul(Q_T, b);
    matrix_t* x = matrix_bs(decomp->R, y);

    matrix_destroy(Q_T);
    matrix_destroy(y);
    matrix_destroy(decomp->Q);
    matrix_destroy(decomp->R);
    free(decomp);

    return x;
}

matrix_t* matrix_clone(matrix_t* A) {
    matrix_t* B = matrix_create(A->m, A->n);
    for (uint64_t i = 0; i < A->m; i++) {
	for (uint64_t j = 0; j < A->n; j++) {
	    matrix_set(B, i, j, matrix_get(A, i, j));
	}
    }
    return B;
}

/* Solves the system Ax = b, where A is an n x n matrix and
 * b is an n-dimensional vector using the iterative Jacobi method,
 * terminating when the residual error is less than the tolerance tol.
 * Returns the residual error, and the solution is stored in x.
 * An initial guess x_0 is required at the solution which will be 
 * iterated upon.
 * The matrix A will be altered in the process.
 */
double matrix_jacobi(matrix_t* A, matrix_t* b, matrix_t** x, double tol) {
    assert(A->m == A->n);
    matrix_t* D_inv = matrix_create(A->m, A->n);
    for (uint64_t i = 0; i < A->m; i++) {
	double a_ii = matrix_get(A, i, i);
	matrix_set(D_inv, i, i, 1.0 / a_ii);
	matrix_set(A, i, i, 0.0);
    } 

    double err = 0.0;
    do {
	matrix_t* x_old = *x;
	// At each iteration, x_k = D^{-1} (b - (L + U)x_{k - 1}).
	matrix_t* inner = matrix_mul(A, x_old);
	matrix_sub_smul(inner, b, 1.0);
	matrix_scalar_mul(inner, -1.0);
	// inner now contains b - (L + U) x_k.
	// We then multiply by the inverse of D.
	matrix_t* x_k = matrix_mul(D_inv, inner);
	matrix_destroy(inner);

	// To calculate the residual we need to find
	// |b - Ax_k|, so we set up the necessary matrices.
	matrix_t* r_k = matrix_clone(b);
	matrix_t* Ax_k = matrix_mul(A, x_k);
	matrix_sub_smul(r_k, Ax_k, 1.0); 
	matrix_destroy(Ax_k);
	err = matrix_norm_frob(r_k);
	matrix_destroy(x_old);
	matrix_destroy(r_k);
	*x = x_k;
	printf("%f", err);
    } while (err > tol);
    matrix_destroy(D_inv);
    return err;
}

