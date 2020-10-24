#include <stdbool.h>
#include <stdint.h>

/**
 * @file
 * @brief Matrix creation, algebra and decompositions.
 */

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
    /** The internal array storing the elements of the matrix.
     * Do not access directly: instead use <tt>matrix_get</tt>.
     */
    double* A; 
    /** The number of rows in the matrix. */
    uint64_t m;
    /** The number of columns in the matrix */
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
/**
 * @brief Solves an upper triangular system of equations using the backward
 * substitution algorithm.
 * 
 * This function operates analogously to <tt>matrix_fs</tt>, except this algorithm
 * can be used to solve systems \f$Ax = b\f$ where \f$A\f$ is upper triangular.
 */
matrix_t* matrix_bs(matrix_t* A, matrix_t* b);
/**
 * @brief Solves a system of linear equations \f$Ax = b\f$ using Gaussian 
 * elimination (without pivoting).
 * 
 * Given an \f$n \times n\f$ matrix @p A and a \f$n \times 1\f$ matrix @p b
 * representing a vector, this function uses Gaussian elimination to find a
 * vector \f$x\f$ such that \f$Ax = b\f$. 
 * 
 * The Gaussian elimination algorithm consists of taking the LU decomposition
 * of A, then using backward and forward substitution because we then obtain
 * lower and upper triangular matrices that we can solve to obtain x.
 *
 * @warning As this matrix uses the LU decomposition algorithm (without
 * pivoting), we require all the <i>principal submatrices</i>, i.e. the 
 * submatrix consisting of the first j rows and columns, to be non-singular
 * for j = 1 to n.
 * Failing to verify this will result in the LU decomposition algorithm
 * producing errors. This function does <b>not</b> check that the condition
 * is satisfied.
 */
matrix_t* matrix_ge(matrix_t* A, matrix_t* b);
/**
 * @brief Solves a system of linear equations \f$Ax = b\f$ using Gaussian
 * elimination (with pivoting). 
 *
 * For an explanation of the Gaussian elimination algorithm, see
 * <tt>matrix_ge</tt>. Instead of using LU factorisation without pivoting,
 * this method uses the LU factorisation with "partial pivoting". This algorithm
 * is designed to deal with small or zero diagonal elements (called pivots). 
 * The algorithm must divide by the pivots, so small or zero pivots lead to
 * significant numerical instabiliy.
 *  
 * The pivoting looks for a permutation matrix \f$P\f$ to maximise the pivot 
 * element by permuting rows, and then calculates the LU decomposition of 
 * \f$PA\f$.
 */
matrix_t* matrix_gepp(matrix_t* A, matrix_t* b);
/** 
 * @brief Initialises a new matrix with @p m rows and @p n columns.
 *
 * This function initialises an \f$m \times n\f$ zero matrix, allocating
 * space for it on the heap and returning a pointer to the newly created
 * matrix.
 * 
 * @warning When you no longer require the matrix, you must use 
 * <tt>matrix_destroy</tt> in order to free the memory used by the matrix.
 */
matrix_t* matrix_create(uint64_t m, uint64_t n);
/**
 * @brief Accesses the element (@p i, @p j) from the matrix @p A.
 */
inline double matrix_get(matrix_t* A, uint64_t i, uint64_t j);
/**
 * @brief Sets the element (@p i, @p j) in the matrix @p A to @p value.
 */
inline void matrix_set(matrix_t* A, uint64_t i, uint64_t j, double value);
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

