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
    double* A; /** The internal array storing the elements of the matrix.
                 * Do not access directly: instead use <tt>matrix_get</tt>.
                 */
    uint64_t m; /** The number of rows in the matrix. */
    uint64_t n; /** The number of columns in the matrix */
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
/** 
 * @brief Solves a lower triangular system of equations using the forward 
 * substitution algorithm.
 *
 * This function can be used to solve systems of linear equations where the
 * matrix \f$L\f$ is a (square) lower triangular matrix and we wish to find
 * \f$x\f$ such that \f$Lx = b\f$.
 *
 * Takes a lower triangular \f$n \times n\f$ matrix @p A and a vector @b
 * (represented as a \f$ n \times 1 \f$ matrix) and returns a newly allocated
 * \f$n \times 1\f$ matrix <tt>x</tt> such that \f$Ax = b\f$.
 *
 * Example: suppose we wish to solve the system \f$\begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}x = (1, 0)^T\f$. The following code will use forward substitution to solve this system.
 * @code
 * matrix_t* A = matrix_id(2);
 * matrix_set(A, 1, 0, 2.0);
 * matrix_t* b = matrix_create(2, 1);
 * matrix_set(b, 0, 0, 1.0);
 * matrix_t* x = matrix_fs(A, b);
 * @endcode
 * 
 */
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

