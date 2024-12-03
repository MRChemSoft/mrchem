/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <fstream>

#include "MRCPP/Printer"

#include "math_utils.h"

namespace mrchem {
namespace math_utils {

DoubleVector init_nan(int I) {
    DoubleVector V(I);
    for (int i = 0; i < I; i++) V(i) = std::numeric_limits<double>::quiet_NaN();
    return V;
}

DoubleMatrix init_nan(int I, int J) {
    DoubleMatrix M(I, J);
    for (int i = 0; i < I; i++) {
        for (int j = 0; j < J; j++) M(i, j) = std::numeric_limits<double>::quiet_NaN();
    }
    return M;
}

/** @brief Calculate the distance between two points in three dimensions */
double calc_distance(const mrcpp::Coord<3> &a, const mrcpp::Coord<3> &b) {
    double r_x = a[0] - b[0];
    double r_y = a[1] - b[1];
    double r_z = a[2] - b[2];
    return std::sqrt(r_x * r_x + r_y * r_y + r_z * r_z);
}

/** @brief Read Eigen matrix from file
 *
 * @param file: file name
 *
 * Format of the file:
 * First entry is the size of the matrix (assumed square).
 * After this all entries of the matrix are listed, columns
 * concatenated into a long vector, single entry per line.
 */
DoubleMatrix read_matrix_file(const std::string &file) {
    int nTerms;
    std::ifstream ifs(file.c_str());
    if (not ifs) MSG_ERROR("Failed to open file: " << file);

    std::string line;
    std::getline(ifs, line);
    std::istringstream iss(line);
    iss >> nTerms;

    DoubleMatrix M = DoubleMatrix::Zero(nTerms, nTerms);
    for (int i = 0; i < nTerms; i++) {
        for (int j = 0; j < nTerms; j++) {
            std::getline(ifs, line);
            std::istringstream iss(line);
            iss >> M(j, i);
        }
    }

    return M;
}

/** @brief Compute the exponential of minus a skew (=antisymmetric) real matrix
 *
 * @param A: matrix to exponentiate
 *
 * The result is a unitary real matrix.
 *	\f$ U=\exp(-A)=\exp(i(iA))\f$ ; \f$ iA=\f$ Hermitian Matrix
 *
 *      \f$ \exp(i(iA))=V\exp(id)V^\dagger \f$ ;
 *      \f$ d \f$: eigenvalues of \f$iA\f$,  \f$V\f$ eigenvectors of
 *      \f$iA\f$ (V unitary complex matrix)
 */
DoubleMatrix skew_matrix_exp(const DoubleMatrix &A) {
    // calculates U=exp(-A)=exp(i(iA)) iA=HermitianMatrix
    // skew=antisymmetric real
    ComplexDouble im(0.0, 1.0);
    ComplexMatrix Aim = im * A;

    // NB: eigenvalues are real, but eigenvectors are complex
    DoubleVector diag;
    ComplexMatrix U = diagonalize_hermitian_matrix(Aim, diag);

    ComplexMatrix diagim = ComplexMatrix::Zero(A.cols(), A.cols());
    for (int j = 0; j < A.cols(); j++) { diagim(j, j) = std::exp(im * diag(j)); }

    Aim = U * diagim * U.adjoint();
    return Aim.real(); // imaginary part is zero
}

/** @brief Compute the power of a Hermitian matrix
 *
 * @param A: matrix
 * @param b: exponent
 *
 * The matrix is first diagonalized, then the diagonal elements are raised
 * to the given power, and the diagonalization is reversed. Sanity check for
 * eigenvalues close to zero, necessary for negative exponents in combination
 * with slightly negative eigenvalues.
 */
ComplexMatrix hermitian_matrix_pow(const ComplexMatrix &A, double b) {
    DoubleVector diag;
    ComplexMatrix U = diagonalize_hermitian_matrix(A, diag);

    DoubleMatrix B = DoubleMatrix::Zero(A.rows(), A.cols());
    for (int i = 0; i < diag.size(); i++) {
        if (std::abs(diag(i)) < mrcpp::MachineZero) {
            B(i, i) = 0.0;
        } else {
            B(i, i) = std::pow(diag(i), b);
        }
    }
    return U * B * U.adjoint();
}

/** @brief Compute the eigenvalues and eigenvectors of a Hermitian matrix
 *
 * @param A: matrix to diagonalize (not modified)
 * @param b: vector to store eigenvalues
 *
 * Returns the matrix of eigenvectors and stores the eigenvalues in the input vector.
 */
ComplexMatrix diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag) {
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(A.cols());
    es.compute(A);
    diag = es.eigenvalues();  // real
    return es.eigenvectors(); // complex
}

/** @brief Compute the eigenvalues and eigenvectors of a Hermitian matrix block
 *
 * @param A: matrix to diagonalize (updated in place)
 * @param U: matrix of eigenvectors
 * @param nstart: upper left corner of block
 * @param nsize: size of block
 *
 * Assumes that the given block is a proper Hermitian sub matrix.
 */
void diagonalize_block(ComplexMatrix &A, ComplexMatrix &U, int nstart, int nsize) {
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(nsize);
    es.compute(A.block(nstart, nstart, nsize, nsize));
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexVector ei_val = es.eigenvalues().cast<ComplexDouble>();
    U.block(nstart, nstart, nsize, nsize) = ei_vec;
    A.block(nstart, nstart, nsize, nsize) = ei_val.asDiagonal();
}

/**
 * @brief Compute the log of the sum of exponentials of the elements of a vector.
 * Numerically stable for large and small x.
 * @param x Vector of values where the logsumexp is to be computed.
 */
double logsumexp(const Eigen::VectorXd &x){
    double max = x.maxCoeff();
    return max + std::log((x.array() - max).exp().sum());
}

} // namespace math_utils
} // namespace mrchem
