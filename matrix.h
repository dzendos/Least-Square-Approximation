#pragma once

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

template <typename T>
class Matrix {
private:
    size_t number_of_rows;
    size_t number_of_columns;

    std::vector<std::vector<T>> values;

    /**
     * @brief resizeMatrix - changes the size of container 'values'
     * @param number_of_rows - new number of rows
     * @param number_of_columns - new number of columns
     * @see values
     */
    void resizeMatrix(size_t number_of_rows, size_t number_of_columns)  {
        values.resize(number_of_rows);

        for (auto &row : values)
            row.resize(number_of_columns);
    }

public:
    Matrix() = default;

    Matrix(size_t number_of_rows, size_t number_of_columns) :
        number_of_rows(number_of_rows),
        number_of_columns(number_of_columns)
    {
        resizeMatrix(number_of_rows, number_of_columns);
    }

    /**
     * @brief  operator [] returns given row of a matrix
     * @param  position of row in the matrix
     * @return row of a matrix
     */
    std::vector<T>& operator [] (const size_t &position) {
        return values[position];
    }

    const std::vector<T>& operator [] (const size_t &position) const {
        return values[position];
    }

    /**
     * @brief operator << prints values of the matrix
     * @param os - output stream
     * @return output stream
     */
    friend std::ostream& operator << (std::ostream &os, const Matrix<T> &matrix) {
        for (int i = 0; i < matrix.number_of_rows; i++) {
            for (int j = 0; j < matrix.number_of_columns; j++) {
                os << matrix.values[i][j] << " ";
            }
            os << std::endl;
        }

        return os;
    }

    /**
     * @brief operator * returns the result of multiplication of 2 matrices
     * @param first - first matrix
     * @param second - second matrix
     * @return the result of multiplication of 2 matrices
     * @throws std::invalid_argument when it is impossible to multiply 2 matrices
     */
    friend Matrix<T> operator * (const Matrix<T> &first, const Matrix<T> &second) {
        if (first.number_of_columns != second.number_of_rows)
            throw std::invalid_argument("cannot multiply 2 matrices");

        Matrix<T> result(first.number_of_rows, second.number_of_columns);

        for (int i = 0; i < result.number_of_rows; i++) {
            for (int j = 0; j < result.number_of_columns; j++)
                for (int k = 0; k < first.number_of_columns; k++)
                    result[i][j] += first[i][k] * second[k][j];

        }

        return result;
    }

    /**
     * @brief resize - changes the size of the matrix
     * @param number_of_rows - new number of rows
     * @param number_of_columns - new number of columns
     */
    void resize(const size_t &number_of_rows, const size_t &number_of_columns)  {
        this->number_of_rows    = number_of_rows;
        this->number_of_columns = number_of_columns;
        resizeMatrix(number_of_rows, number_of_columns);
    }

    /**
     * @brief transposed - puts transpose of the current matrix in parameter matrix
     * @param transposed_matrix - matrix that contains result of transposing
     */
    void transposed(Matrix<T> &transposed_matrix) {
        transposed_matrix.resize(number_of_columns, number_of_rows);

        for (int i = 0; i < number_of_rows; i++)
            for (int j = 0; j < number_of_columns; j++)
                transposed_matrix[j][i] = values[i][j];
    }

    void inversed(Matrix<T> &inversed_matrix) {
        if (inversed_matrix.number_of_rows != inversed_matrix.number_of_columns ||
                number_of_rows    != inversed_matrix.number_of_rows             ||
                number_of_columns != inversed_matrix.number_of_columns)
        {
            throw std::invalid_argument("cannot inverse matrix");
        }

        Matrix<double> matrix_copy = *this;

        // making identity matrix
        for (int i = 0; i < inversed_matrix.number_of_rows; i++) {
            for (int j = 0; j < inversed_matrix.number_of_columns; j++)
                inversed_matrix[i][j] = (i == j ? 1 : 0);
        }

        // getting inversed matrix
        for (int i = 0; i < inversed_matrix.number_of_rows; i++) {
            for (int j = 0; j < inversed_matrix.number_of_rows; j++) {
                if (i == j) // the same row
                    continue;

                double coefficient = matrix_copy[j][i] / matrix_copy[i][i];
                for (int k = 0; k < matrix_copy.number_of_columns; k++) {
                    matrix_copy[j][k]     -= matrix_copy[i][k]     * coefficient;
                    inversed_matrix[j][k] -= inversed_matrix[i][k] * coefficient;
                }

                matrix_copy[j][i] = 0;
            }
        }

        for (int i = 0; i < inversed_matrix.number_of_rows; i++) {
            double coefficient = matrix_copy[i][i];
            for (int k = 0; k < inversed_matrix.number_of_columns; k++) {
                inversed_matrix[i][k] /= coefficient;
            }
        }
    }

    size_t getNumberOfRows() {
        return number_of_rows;
    }

    size_t getNumberOfColumns() {
        return number_of_columns;
    }
};


















