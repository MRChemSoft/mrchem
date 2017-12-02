#include <cmath>
#include <cstdio>
#include <fstream>

#include "MathUtils.h"
#include "TelePrompter.h"
#include "constants.h"
#include "parallel.h"
#include "eigen_disable_warnings.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

/** Compute the norm of a matrix given as a vector
 *
 * The norm of the matrix is computed by iterating the following operation:
 *	\f$ x_n = M^t \cdot M \cdot x_{n-1} \f$
 *
 *	The norm of the matrix is obtained as:
 *	 \f$ ||M|| \lim_{n \rightarrow \infty} ||x_n||/||x_{n-1}||\f$
 */
double MathUtils::matrixNorm2(const VectorXd &vector) {
    int n = (int) sqrt(vector.size());
    MatrixXd matrix = MatrixXd::Map(vector.data(), n, n);
    return matrixNorm2(matrix);
}

/** Compute the norm of a matrix */
double MathUtils::matrixNorm2(const MatrixXd &M) {
    const int jMax = 100;
    const int maxTry = 10;
    int size = M.cols();

    const double tol = 1.0e-2;
    double newNorm;
    double newOperNorm;
    double oldOperNorm;
    double ratio;

    Eigen::VectorXd x(size);
    Eigen::VectorXd y(size);

    std::srand(1);
    for (int i = 0; i < maxTry; i++) {
        y = Eigen::VectorXd::Random(size);
        newNorm = sqrt(y.squaredNorm());

        ratio = 0.0;
        oldOperNorm = 0.0;
        newOperNorm = 0.0;

        for (int j = 0; j <= jMax; j++) {
            x = y / newNorm;
            y = M.transpose() * M * x;
            newNorm = sqrt(y.squaredNorm());
            oldOperNorm = newOperNorm;
            newOperNorm = sqrt(newNorm);
            if (oldOperNorm == 0.0) {
                continue;
            }
            ratio = newOperNorm/oldOperNorm;
            if (fabs(ratio - 1.0) <= tol) {
                return newOperNorm;
            }
        }
        MSG_WARN("Random guess did not converge, starting over!");
    }
    MSG_FATAL("Could not calculate matrix norm!");
    return -1.0;
}

/** Compute the norm of a matrix given as a vector.
 *
 * The norm of the matrix is obtained by taking the column with the
 * largest norm.
 */
double MathUtils::matrixNorm1(const VectorXd &vector) {
    int n = (int) sqrt(vector.size());
    MatrixXd matrix = MatrixXd::Map(vector.data(), n, n);
    return matrixNorm1(matrix);
}

/** Compute the norm of a matrix.
*/
double MathUtils::matrixNorm1(const MatrixXd &M) {
    int nRows = M.rows();
    int nCols = M.rows();
    double maxNorm = 0.0;
    for (int i = 0; i < nCols; i++) {
        double colNorm = 0.0;
        for (int j = 0; j < nRows; j++) {
            colNorm += fabs(M(j,i));
        }
        if (colNorm > maxNorm) {
            maxNorm = colNorm;
        }
    }
    return maxNorm;
}


/** Compute the infinity norm of a matrix given as a vector.
 * The norm of the matrix is obtained by taking the row with the largest norm.
 */
double MathUtils::matrixNormInfinity(const VectorXd &vector) {
    int n = (int) sqrt(vector.size());
    MatrixXd matrix = MatrixXd::Map(vector.data(), n, n);
    return matrixNormInfinity(matrix);
}

/** Compute the infinity norm of a matrix.
 * The norm of the matrix is obtained by taking the row with the largest norm.
 */
double MathUtils::matrixNormInfinity(const MatrixXd &M) {
    int nRows = M.rows();
    int nCols = M.rows();
    double maxNorm = 0.0;
    for (int i = 0; i < nRows; i++) {
        double rowNorm = 0.0;
        for (int j = 0; j < nCols; j++) {
            rowNorm += fabs(M(i,j));
        }
        if (rowNorm > maxNorm) {maxNorm = rowNorm;}
    }
    return maxNorm;
}

/** Compute the binomial coefficient n!/((n-j)! j!) */
double MathUtils::binomialCoeff(int n, int j) {
    double binomial_n_j = 1.0;
    int k = 0;

    if (n < 0 || j < 0 || j > n) {
        MSG_ERROR("Negative argument or j > n is not defined.");
    } else {
        k = 0;
        while (k < j) {
            binomial_n_j *= (double) (n - k);
            k += 1;
        }
        binomial_n_j /= factorial(j);
    }
    return binomial_n_j;
}

/** Compute k! = GAMMA(k+1) for integer argument k */
double MathUtils::factorial(int n) {
    int k = 1;
    double fac_n = 1.0;

    if (n < 0) {
        MSG_FATAL("Negative argument is not defined.");
    } else if (n == 0 || n == 1)
        return 1.0;
    else if (n > 1) {
        while (k <= n) {
            fac_n *= (double) k;
            k += 1;
        }
    }
    return fac_n;
}

/** Compute the tensor product of two matrices */
MatrixXd MathUtils::tensorproduct(const MatrixXd &A,
                                  const MatrixXd &B) {
    int Ar = A.rows();
    int Ac = A.cols();
    int Br = B.rows();
    int Bc = B.cols();
    MatrixXd tprod(Ar * Br, Ac * Bc);
    for (int i = 0; i < Ar; i++) {
        for (int j = 0; j < Ac; j++) {
            tprod.block(i * Br, j * Bc, Br, Bc) = A(i, j) * B;
        }
    }
    return tprod;
}

/** Compute the tensor product of a matrix and a vector */
MatrixXd MathUtils::tensorproduct(const MatrixXd &A,
                                  const VectorXd &B) {
    int Ar = A.rows();
    int Ac = A.cols();
    int Br = B.rows();
    MatrixXd tprod(Ar * Br, Ac);
    for (int i = 0; i < Br; i++) {
        tprod.block(i * Br, 0, Ar, Ac) = A * B(i);
    }
    return tprod;
}

/** Compute the tensor product of a matrix and a vector */
MatrixXd MathUtils::tensorproduct(const VectorXd &A,
                                  const MatrixXd &B) {
    int Ar = A.rows();
    int Br = B.rows();
    int Bc = B.cols();
    MatrixXd tprod(Ar * Br, Ar);
    for (int i = 0; i < Ar; i++) {
        tprod.block(i * Br, 0, Br, Bc) = A(i) * B;
    }
    return tprod;
}

/** Compute the tensor product of a column vector and a row vector */
MatrixXd MathUtils::tensorproduct(const VectorXd &A,
                                  const VectorXd &B) {
    int Ar = A.rows();
    int Br = B.rows();
    MatrixXd tprod(Ar, Br);
    for (int i = 0; i < Ar; i++) {
        tprod.block(i, 0, 1, Br) = A(i) * B.transpose();
    }
    return tprod;
}

/** Compute the tensor product of a vector and itself */
void MathUtils::tensorSelfProduct(const VectorXd &A, VectorXd &tprod) {
    int Ar = A.rows();
    for (int i = 0; i < Ar; i++) {
        tprod.segment(i*Ar, Ar) = A(i) * A;
    }
}

/** Compute the tensor product of a vector and itself */
void MathUtils::tensorSelfProduct(const VectorXd &A, MatrixXd &tprod) {
    int Ar = A.rows();
    for (int i = 0; i < Ar; i++) {
        tprod.block(i, 0, 1, Ar) = A(i) * A;
    }
}

// for debugging
void MathUtils::print_vector(int n, const double *vec) {
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(vec[i]) < 1.0e-10) {
            printf("%e\n", 0.0);
        } else {
            printf("%e\n", vec[i]);
        }
    }
    printf("\n");
}

VectorXd MathUtils::getBinomialCoefs(unsigned int order) {
    VectorXd coefs = VectorXd::Ones(order + 1);
    for (int k = 0; k <= order; k++) {
        coefs[k] = MathUtils::binomialCoeff(order, k);
    }
    return coefs;
}

void MathUtils::applyFilter(double *out, double *in,
                            const MatrixXd &filter,
                            int kp1, int kp1_dm1, double fac) {
#ifdef HAVE_BLAS
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                kp1_dm1, kp1, kp1, 1.0, in, kp1, filter.data(),
                kp1, fac, out, kp1_dm1);
#else
    Eigen::Map<MatrixXd> f(in, kp1, kp1_dm1);
    Eigen::Map<MatrixXd> g(out, kp1_dm1, kp1);
    if (fac < MachineZero) {
        g = f.transpose() * filter;
    } else {
        g += f.transpose() * filter;
    }
#endif
}

/** Make a nD-representation from 1D-representations of separable functions.
 *
 * This method uses the "output" vector as initial input, in order to
 * avoid the use of temporaries.
 */
void MathUtils::tensorExpandCoefs(int dim, int dir, int kp1, int kp1_d,
                                  const MatrixXd &primitive, VectorXd &expanded) {
    if (dir < dim - 1) {
        int idx = MathUtils::ipow(kp1, dir + 1);
        int nelem = idx * kp1;
        int pos = kp1_d - nelem;
        int inpos = kp1_d - idx;
        for (int i = 0; i < kp1; i++) {
            expanded.segment(pos + i * idx, idx) =
                    expanded.segment(inpos, idx) * primitive.col(dir + 1)(i);
        }
        tensorExpandCoefs(dim, dir + 1, kp1, kp1_d, primitive, expanded);
    }
}

void MathUtils::tensorExpandCoords_2D(int kp1, const MatrixXd &primitive, MatrixXd &expanded) {
    int n = 0;
    for (int i = 0; i < kp1; i++) {
        for (int j = 0; j < kp1; j++) {
            expanded(0,n) = primitive(0,j);
            expanded(1,n) = primitive(1,i);
            n++;
        }
    }
}

void MathUtils::tensorExpandCoords_3D(int kp1, const MatrixXd &primitive, MatrixXd &expanded) {
    int n = 0;
    for (int i = 0; i < kp1; i++) {
        for (int j = 0; j < kp1; j++) {
            for (int k = 0; k < kp1; k++) {
                expanded(0,n) = primitive(0,k);
                expanded(1,n) = primitive(1,j);
                expanded(2,n) = primitive(2,i);
                n++;
            }
        }
    }
}

/** Calculate the distance between two points in n-dimensions */
double MathUtils::calcDistance(int D, const double *a, const double *b) {
    assert(a != 0 and b != 0 and D >= 0);
    double r = 0.0;
    for (int i = 0; i < D; i++) {
        r += pow(a[i] - b[i], 2.0);
    }
    return sqrt(r);
}

void MathUtils::swapCols(Matrix<double, Dynamic, Dynamic, RowMajor> &mat, int i, int j) {
    if (i == j) {
        return;
    }
    VectorXd tmpVec = mat.col(j);
    mat.col(j) = mat.col(i);
    mat.col(i) = tmpVec;
}

void MathUtils::swapRows(Matrix<double, Dynamic, Dynamic, RowMajor> &mat, int i, int j) {
    if (i == j) {
        return;
    }
    VectorXd tmpVec = mat.row(j);
    mat.row(j) = mat.row(i);
    mat.row(i) = tmpVec;
}

Eigen::MatrixXd MathUtils::hermitianMatrixPow(Eigen::MatrixXd A, double b) {
    MatrixXd U = diagonalizeHermitianMatrix(A);
    for (int i = 0; i < A.cols(); i++) {
        A(i,i) = pow(A(i,i), b);
    }
    return U * A * U.transpose();
}

Eigen::MatrixXd MathUtils::diagonalizeHermitianMatrix(Eigen::MatrixXd &A) {
    SelfAdjointEigenSolver<MatrixXd> es(A.cols());
    es.compute(A);
    A = es.eigenvalues().asDiagonal();
    return es.eigenvectors();
}


/** Compute the exponential of minus a skew (=antisymmetric) real matrix \f$A\f$
 *
 * The result is a unitary real matrix.
 *	\f$ U=\exp(-A)=\exp(i(iA))\f$ ; \f$ iA=\f$ Hermitian Matrix
 *
 *      \f$ \exp(i(iA))=V\exp(id)V^\dagger \f$ ;
 *      \f$ d \f$: eigenvalues of \f$iA\f$,  \f$V\f$ eigenvectors of
 *      \f$iA\f$ (V unitary complex matrix)
 */
Eigen::MatrixXd MathUtils::SkewMatrixExp(Eigen::MatrixXd &A) {
    //calculates U=exp(-A)=exp(i(iA)) iA=HermitianMatrix
    //skew=antisymmetric real
    complex<double> im(0.0,1.0);
    MatrixXcd Aim(A.cols(),A.cols()) ;
    VectorXd diag(A.cols());
    MatrixXcd diagim(A.cols(),A.cols()) ;
    Aim = im*A;
    MatrixXcd U(A.cols(),A.cols()) ;

    //NB: eigenvalues are real, but eigenvectors are complex
    U = diagonalizeHermitianMatrix(Aim,diag);

    for (int j = 0; j < A.cols(); j++) {
        for (int i = 0; i < A.cols(); i++) {
            diagim(i,j) = 0.0;
        }
        diagim(j,j) = exp(im*diag(j));
    }

    Aim = U * diagim * U.adjoint();
    return Aim.real() ;  //imaginary part is zero
}

Eigen::MatrixXcd MathUtils::diagonalizeHermitianMatrix(Eigen::MatrixXcd &A,Eigen::VectorXd &diag) {
    SelfAdjointEigenSolver<MatrixXcd> es(A.cols());
    es.compute(A);
    diag = es.eigenvalues();//real
    return es.eigenvectors();//complex
}

MatrixXd MathUtils::readMatrixFile(const string &file) {
    int nTerms;
    fstream ifs;
    ifs.open(file.c_str());
    if (not ifs) {
        MSG_ERROR("Failed to open file: " << file);
    }
    string line;
    getline(ifs, line);
    istringstream iss(line);
    iss >> nTerms;

    MatrixXd matrix = MatrixXd::Zero(nTerms, nTerms);
    for (int i = 0; i < nTerms; i++) {
        for (int j = 0; j < nTerms; j++) {
            getline(ifs, line);
            istringstream iss(line);
            iss >> matrix(j,i);
        }
    }
    return matrix;
}

void MathUtils::diagonalizeBlock(Eigen::MatrixXd &M , Eigen::MatrixXd &U, int nstart, int nsize) {
    SelfAdjointEigenSolver<MatrixXd> es(nsize);
    es.compute(M.block(nstart, nstart, nsize, nsize));
    MatrixXd tmp = es.eigenvectors();
    U.block(nstart, nstart, nsize, nsize) = tmp.transpose();
    M.block(nstart, nstart, nsize, nsize) = es.eigenvalues().asDiagonal();
}
