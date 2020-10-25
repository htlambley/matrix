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

/**
 * @brief Tests whether a matrix @p A is lower triangular. 
 */
bool matrix_lower_triangular(matrix_t* A); 
/** 
 * @brief Tests whether a matrix @p A is upper triangular.
 */
bool matrix_upper_triangular(matrix_t* A);
/**
 * @brief Tests whether a matrix @p A is diagonal.
 */
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
 * Example: suppose we wish to solve the system \f$Ax = (1, 0)^T\f$, where A = [1, 0; 2, 1]. The following code will use forward substitution to solve this system.
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
 * @brief Accesses the element (@p i, @p j) from the matrix @p A,
 * where @p i accesses the (i + 1)th row and @p j accesses the (j + 1)th
 * column.
 * 
 * @warning In contrast to the mathematical notation \f$a_{11}\f$ as the
 * first element of the matrix, we denote the top left element by the index
 * (0, 0), and the bottom right by (m - 1, n - 1) in an m x n matrix. 
 */
static inline double matrix_get(matrix_t* A, uint64_t i, uint64_t j);
/**
 * @brief Sets the element (@p i, @p j) in the matrix @p A to @p value.
 */
static inline void matrix_set(matrix_t* A, uint64_t i, uint64_t j, double value);
/**
 * @brief Adds the matrices @p A and @p B, storing the result in the matrix
 * @p A. 
 *
 * @warning Note carefully that this will overwrite what is stored in the
 * matrix @p A, and leave @p B unchanged. This is to avoid allocating another
 * matrix unnecessarily. If you need to keep the original @p A, allocate a
 * copy using <tt>matrix_clone</tt>.
 */
void matrix_add(matrix_t* A, matrix_t* B);
/**
 * @brief Creates a new identity matrix of dimension @p n.
 *
 * This function allocates a new matrix that is equal to the identity
 * matrix \f$I_n\f$ with 1 on the diagonal and zero elsewhere. As
 * with <tt>matrix_create</tt>, you must destroy this manually using
 * <tt>matrix_destroy</tt> after usage.
 */
matrix_t* matrix_id(uint64_t n);
/**
 * @brief Creates a new matrix equal to the product of @p A and @p B.
 *
 * The matrices must be of compatible dimension, so if \f$A\f$ is an
 * \f$m \times n\f$ matrix, then \f$B\f$ must be a \f$n \times p\f$ matrix
 * and the function will return a newly allocated \f$m \times p\f$ matrix.
 * 
 * Unlike the <tt>matrix_add</tt> function, this function allocates a new
 * matrix because, unless @p A and @p B are square, the product matrix will
 * not be of the same dimension as @p A or @p B.
 */
matrix_t* matrix_mul(matrix_t* A, matrix_t* B);
/**
 * @brief Computes the LU decomposition of the given matrix.
 *
 * See also <tt>matrix_ge</tt> for discussion of this function. This
 * performs the LU decomposition algorithm without any pivoting. 
 *
 * @returns A unit lower triangular matrix (i.e. a matrix with 1 on the
 * diagonal and zero above). The upper triangular matrix will be stored in the
 * matrix passed to the function.
 *
 * @warning This algorithm requires that @p U has all <i>principal 
 * submatrices</i> non-singular, as described in <tt>matrix_ge</tt>.
 * You must be sure of this before running the algorithm.
 */
matrix_t* matrix_lu(matrix_t* U);
/**
 * @brief Computes the LU decomposition algorithm (with partial pivoting of
 * rows).
 *
 * Analogously to <tt>matrix_lu</tt>, this function can be used to find matrices
 * \f$L\f$ and \f$U\f$, where \f$L\f$ is unit lower triangular and \f$U\f$ is
 * upper triangular.
 *
 * In contrast to <tt>matrix_lu</tt>, this function also finds a permutation 
 * matrix \f$P\f$ such that \f$PA = LU\f$. This approach has the advantage of
 * being able to work with matrices not satisfying the non-singular submatrix
 * criterion. As discussed in <tt>matrix_gepp</tt>, pivoting the rows of the
 * matrix prevents numerical instability for small (or zero) pivot elements.
 * @returns A structure containing the matrices L and P. The upper triangular
 * matrix is stored in the argument passed to the function.
 */
matrix_lup_t* matrix_lupp(matrix_t* U);
/**
 * @brief Swaps the rows @p k and @p l in the matrix @p A.
 */
void matrix_row_swap(matrix_t* A, uint64_t k, uint64_t l);
/** 
 * @brief Prints a compact representation of the matrix to stdout.
 */
void matrix_print(matrix_t* A);
/**
 * @brief Computes a reduced QR factorisation of the matrix @p A using
 * the Gram--Schmidt orthonormalisation procedure.
 */
matrix_qr_t* matrix_reduced_qr(matrix_t* A);
/**
 * @brief Computes the Householder reflection with normal vector @p v.
 * 
 * This function creates a new matrix \f$ H = I_n - 2 v v^T \f$ which is
 * called a <i>Householder reflection</i>. This is a reflection by the 
 * hyperplane normal to @p v, which should be given as an \f$ n \times 1\f$
 * matrix.
 */
matrix_t* matrix_householder(matrix_t* v);
/**
 * @brief Computes a QR decomposition of @p A using the Householder
 * reflections method.
 *
 * This is an alternative to computing the QR decomposition with 
 * <tt>matrix_reduced_qr</tt>, which uses the Gram--Schmidt orthonormalisation
 * procedure. The Householder reflections method tends to be numerically more
 * stable than the Gram--Schmidt method.
 *
 * @returns A heap-allocated structure storing newly allocated matrices Q and R
 * as calculated by the algorithm. The original matrix @p A is left unchanged.
 */
matrix_qr_t* matrix_qr_hh(matrix_t* A);
/**
 * @brief Creates a new matrix equal to the transpose of @p A.
 *
 * The transpose of an \f$ m \times n \f$ matrix is an \f$ n \times m \f$
 * matrix where \f$ a_{ij} = b_{ji} \f$ where \f$ b_{ij} \f$ are the elements
 * of the transpose matrix.
 */
matrix_t* matrix_transpose(matrix_t* A);
/**
 * @brief Frees the memory used by the matrix @p A.
 *
 * This function frees the array of doubles used to store the elements
 * of the matrix and then frees the matrix structure itself.
 */
void matrix_destroy(matrix_t* A);
/** 
 * @brief Finds a vector x minimising the function 
 * \f$\frac12 \lVert Ax - b \rVert^2\f$.
 *
 * This function can be used to obtain a minimiser of the function \f$g(x) =
 * \frac12 \lVert Ax - b \rVert^2\f$ given an \f$m \times n\f$ matrix @p A
 * and a \f$m \times 1\f$ matrix @p b. These problems are known as <b>least
 * squares problems</b>.
 * 
 * Internally, this function calculates the reduced QR factorisation of @p A
 * and then solves the system using backward substitution.
 *
 * @warning We require @p A to be of full rank or else the solution is not
 * uniquely determined to the least squares problem. 
 */
matrix_t* matrix_lsq(matrix_t* A, matrix_t* b);
/**
 * @brief Iteratively solves the system \f$Ax = b\f$ where @p A satisfies 
 * certain convergence criteria.
 *
 * In contrast to <tt>matrix_ge</tt> and related functions, this is an 
 * iterative solver to find the solution to a system of linear equations.
 *
 * As this is an iterative solver, a tolerance @p tol must be passed which
 * sets when to stop: when the residual \f$\lVert Ax - b \rVert \f$ is less 
 * than @p tol, we terminate.
 *
 * This function will never iterate more than MAX_ITER times as defined in 
 * <tt>matrix.c</tt> to prevent infinite looping if the process does not 
 * converge.
 *
 * @returns The 2-norm of the residual \f$Ax - b\f$ for the vector @p x which
 * stores the solution found.
 *
 * @warning The iterative method will only converge if certain criteria are
 * satisified by @p A. If these are not satisfied, the iterations may never
 * converge to the true value. The following are sufficient conditions for the
 * convergence of the Jacobi method:
 * - the 'strong row sum' criterion, 
 *   \f$|a_{ii}| > \sum_{j \neq i} |a_{ij}|\f$ for all \f$i\f$
 * - the 'strong column sum' criterion, formulated analogously to above,
 * - irreducibility and the 'weak row sum' criterion
 * \f$|a_{ii}| \geq \sum_{j \neq i} |a_{ij}| for all i, and
 * \f$|a_{kk}| > \sum_{j \neq k} |a_{kj}| for some row k.
 */
double matrix_jacobi(matrix_t* A, matrix_t* b, matrix_t** x, double tol);


