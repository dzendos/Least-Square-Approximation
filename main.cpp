#include <iostream>
#include <fstream>
#include <iomanip>

#include "matrix.h"

using namespace std;

void input(istream &is, Matrix<double> &matrix_a, Matrix<double> &matrix_b) {
    int n, m;

    is >> n >> m;

    matrix_a.resize(m, n + 1); // m by n+1 matrix
    matrix_b.resize(m, 1);     // m by 1   vector

    for (int i = 0; i < m; i++) {
        matrix_a[i][0] = 1;

        for (int j = 0; j < n; j++)
            is >> matrix_a[i][j + 1];

        is >> matrix_b[i][0];
    }
}

void output(ostream &os, const Matrix<double> &matrix_a, const Matrix<double> &matrix_b,
            const Matrix<double> &matrix_a_transposed_m_a, const Matrix<double> &inversed_a_t_mult_a,
            const Matrix<double> &inversed_a_t_mult_a__mult_a_transposed, const Matrix<double> &matrix_x)
{
    os << std::fixed << setprecision(2);

    os << "A:"              << std::endl << matrix_a                               << std::endl;
    os << "b:"              << std::endl << matrix_b                               << std::endl;
    os << "A_T*A:"          << std::endl << matrix_a_transposed_m_a                << std::endl;
    os << "(A_T*A)_-1:"     << std::endl << inversed_a_t_mult_a                    << std::endl;
    os << "(A_T*A)_-1*A_T:" << std::endl << inversed_a_t_mult_a__mult_a_transposed << std::endl;
    os << "x:"              << std::endl << matrix_x                               << std::endl;
}

void run() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    Matrix<double> matrix_a;
    Matrix<double> matrix_b;
    Matrix<double> matrix_a_transposed;

    // getting matrices a and b
    input(fin, matrix_a, matrix_b);

    // getting matrix a_transposed and a_transposed multiplied by a
    matrix_a.transposed(matrix_a_transposed);
    Matrix<double> matrix_a_transposed_m_a = matrix_a_transposed * matrix_a;

    // getting previous matrix inversed
    Matrix<double> inversed_a_t_mult_a(matrix_a_transposed_m_a.getNumberOfRows(),
                                       matrix_a_transposed_m_a.getNumberOfColumns());
    matrix_a_transposed_m_a.inversed(inversed_a_t_mult_a);

    // getting previous matrix multiplied by a_transposed
    Matrix<double> inversed_a_t_mult_a__mult_a_transposed = inversed_a_t_mult_a * matrix_a_transposed;

    // getting the answer
    Matrix<double> matrix_x = inversed_a_t_mult_a__mult_a_transposed * matrix_b;

    output(fout, matrix_a, matrix_b, matrix_a_transposed_m_a, inversed_a_t_mult_a, inversed_a_t_mult_a__mult_a_transposed, matrix_x);

    fin.close();
    fout.close();
}

int main() {
    run();

    return 0;
}
