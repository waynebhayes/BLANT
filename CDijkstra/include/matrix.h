#ifndef CDIJKSTRA_MATRIX_HEADER
#define CDIJKSTRA_MATRIX_HEADER

namespace cdijkstra {

template<typename T>
class Matrix {
public:
    Matrix(unsigned int rows, unsigned int cols, T default_val);
    Matrix(const Matrix<T>& m);
    Matrix& operator =(const Matrix<T>& m);
    ~Matrix();

    unsigned int rows() const;
    unsigned int cols() const;

    const T& at(unsigned int n, unsigned int m) const;
    T& at(unsigned int n, unsigned int m);

private:
    unsigned int rows_;
    unsigned int cols_;
    T** mat;
};

template<typename T>
Matrix<T>::Matrix(unsigned int r, unsigned int c, T default_val)
 : rows_{c}, cols_{c}, mat{new T*[rows()]} {
    for (unsigned int i = 0; i < rows(); ++i) {
        mat[i] = new T[cols()];

        for (unsigned int j = 0; j < cols(); ++j) {
            mat[i][j] = default_val;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& m)
 : rows_{m.rows()}, cols_{m.cols()}, mat{new T*[rows()]} {
    for (unsigned int i = 0; i < rows(); ++i) {
        mat[i] = new T[cols()];

        for (unsigned int j = 0; j < cols(); ++j) {
            mat[i][j] = m.at(i, j);
        }
    }
}

template<typename T>
Matrix<T>& Matrix<T>::operator =(const Matrix<T>& m) {
    // TODO: code reuse could be better here
    if (this != &m) {
        for (unsigned int i = 0; i < rows(); ++i) {
            delete[] mat[i];
        }

        delete[] mat;

        rows_ = m.rows();
        cols_ = m.cols();
        mat = new T*[rows()];

        for (unsigned int i = 0; i < rows(); ++i) {
            mat[i] = new T[cols()];

            for (unsigned int j = 0; j < cols(); ++j) {
                mat[i][j] = m.at(i, j);
            }
        }
    }

    return *this;
}

template<typename T>
Matrix<T>::~Matrix() {
    for (unsigned int i = 0; i < rows(); ++i) {
        delete[] mat[i];
    }

    delete[] mat;
}

template<typename T>
inline unsigned int Matrix<T>::rows() const {
    return rows_;
}

template<typename T>
inline unsigned int Matrix<T>::cols() const {
    return cols_;
}

template<typename T>
inline const T& Matrix<T>::at(unsigned int i, unsigned int j) const {
    return mat[i][j];
}

template<typename T>
inline T& Matrix<T>::at(unsigned int i, unsigned int j) {
    return mat[i][j];
}

}

#endif
