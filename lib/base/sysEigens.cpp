#include "sysEigens.hpp"

double find_max_evals(Scalar2D& raw_matrix){
    const int N = raw_matrix.size();

    Eigen::MatrixXd matrix;
    matrix.resize(N, N);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            matrix(i, j) = raw_matrix[i][j];
        }
    }

    cout << matrix << endl;

    Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    return eigenvalues.cwiseAbs().maxCoeff();
}
