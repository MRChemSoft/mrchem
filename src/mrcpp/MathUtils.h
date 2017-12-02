/* \file MathUtils.h
 *
 * \breif Collection of misc math funcs.
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#pragma GCC system_header
#include <Eigen/Eigenvalues>

#include "TelePrompter.h"

class MathUtils {
public:
    static double binomialCoeff(int n, int j);
    static double factorial(int n);
    /** Calculate \f$ m^e\f$ for integers. */
    static int ipow(int m, int e) {
        if (e < 0) {
            MSG_FATAL("Exponent cannot be negative: " << e)
        }
        int result = 1;
        for (int i=0; i < e; i++) {
            result *= m;
        }
        return result;
    }

    static Eigen::MatrixXd hermitianMatrixPow(Eigen::MatrixXd A, double b);
    static Eigen::MatrixXd diagonalizeHermitianMatrix(Eigen::MatrixXd &A);

    static Eigen::MatrixXd tensorproduct(const Eigen::MatrixXd &A,
                                         const Eigen::MatrixXd &B);
    static Eigen::MatrixXd tensorproduct(const Eigen::MatrixXd &A,
                                         const Eigen::VectorXd &B);
    static Eigen::MatrixXd tensorproduct(const Eigen::VectorXd &A,
                                         const Eigen::MatrixXd &B);
    static Eigen::MatrixXd tensorproduct(const Eigen::VectorXd &A,
                                         const Eigen::VectorXd &B);

    static void tensorSelfProduct(const Eigen::VectorXd &A, Eigen::VectorXd &B);
    static void tensorSelfProduct(const Eigen::VectorXd &A, Eigen::MatrixXd &B);

    static double matrixNormInfinity(const Eigen::MatrixXd &M);
    static double matrixNorm2(const Eigen::MatrixXd &M);
    static double matrixNorm1(const Eigen::MatrixXd &M);

    static double matrixNormInfinity(const Eigen::VectorXd &M);
    static double matrixNorm2(const Eigen::VectorXd &M);
    static double matrixNorm1(const Eigen::VectorXd &M);

    static Eigen::VectorXd getBinomialCoefs(unsigned int order);

    static void applyFilter(double *out, double *in, const Eigen::MatrixXd &filter,
                            int kp1, int kp1_dm1, double fac);

    static void tensorExpandCoefs(int dim, int dir, int kp1, int kp1_d,
                                  const Eigen::MatrixXd &primitive, Eigen::VectorXd &expanded);

    static void tensorExpandCoords_2D(int kp1,
                                      const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);
    static void tensorExpandCoords_3D(int kp1,
                                      const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);

    static double calcDistance(int D, const double *a, const double *b);

    static void swapRows(Eigen::Matrix<double, Eigen::Dynamic,
                         Eigen::Dynamic, Eigen::RowMajor> &mat, int i, int j);
    static void swapCols(Eigen::Matrix<double, Eigen::Dynamic,
                         Eigen::Dynamic, Eigen::RowMajor> &mat, int i, int j);

    static Eigen::MatrixXd SkewMatrixExp(Eigen::MatrixXd &A);
    static Eigen::MatrixXcd diagonalizeHermitianMatrix(Eigen::MatrixXcd &A,Eigen::VectorXd &diag);

    static Eigen::MatrixXd readMatrixFile(const std::string &file);

    static void diagonalizeBlock(Eigen::MatrixXd &M , Eigen::MatrixXd &U, int nstart, int nsize);

private:
    static void print_vector(int n, const double *vec);
};

