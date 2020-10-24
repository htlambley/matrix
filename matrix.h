#include <stdbool.h>
#include <stdint.h>

/**
 * @brief The main matrix type.
 *
 * This is the main matrix structure used in the library. Most of the the time,
 * you should expect to be given a pointer to a heap-allocated version of this
 * structure. 
 *
 * The alias <tt>matrix_t</tt> for this structure is defined for convienience.
 * 
 * To create a matrix, see the <tt>matrix_create</tt> function. 
 */
struct matrix {
    double* A;
    uint64_t m;
    uint64_t n;
};

typedef struct matrix matrix_t;

struct matrix_lup {
    matrix_t* L;
    matrix_t* P;
};

struct matrix_qr {
    matrix_t* Q;
    matrix_t* R;
};

typedef struct matrix_lup matrix_lup_t;
typedef struct matrix_qr matrix_qr_t;

bool matrix_lower_triangular(matrix_t* A); 
bool matrix_upper_triangular(matrix_t* A);
bool matrix_diagonal(matrix_t* A);
matrix_t* matrix_fs(matrix_t* A, matrix_t* b);
matrix_t* matrix_bs(matrix_t* A, matrix_t* b);
matrix_t* matrix_ge(matrix_t* A, matrix_t* b);
matrix_t* matrix_gepp(matrix_t* A, matrix_t* b);
matrix_t* matrix_create(uint64_t m, uint64_t n);
double matrix_get(matrix_t* m, uint64_t i, uint64_t j);
void matrix_set(matrix_t* m, uint64_t i, uint64_t j, double value);
void matrix_add(matrix_t* A, matrix_t* B);
matrix_t* matrix_id(uint64_t n);
matrix_t* matrix_mul(matrix_t* A, matrix_t* B);
matrix_t* matrix_lu(matrix_t* U);
matrix_lup_t* matrix_lupp(matrix_t* U);
void matrix_row_swap(matrix_t* A, uint64_t k, uint64_t l);
void matrix_print(matrix_t* A);
matrix_qr_t* matrix_reduced_qr(matrix_t* A);
matrix_t* matrix_transpose(matrix_t* A);
void matrix_destroy(matrix_t* A);
matrix_t* matrix_lsq(matrix_t* A, matrix_t* b);
double matrix_jacobi(matrix_t* A, matrix_t* b, matrix_t** x, double tol);

