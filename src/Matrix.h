#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <iostream>

namespace matrix {
using OrdinalType = std::size_t;

enum class MatCtorType { ZEROS, IDENTITY };

template <typename ScalarType>
class Matrix {
 protected:
  std::size_t nrows_;
  std::vector<ScalarType> mat_;

 public:
  Matrix() : nrows_(0) {}
  Matrix(const std::size_t nrows, const MatCtorType c = MatCtorType::ZEROS)
      : nrows_(nrows), mat_(nrows * nrows, 0) {
    if (c == MatCtorType::IDENTITY) {
      for (std::size_t i = 0; i < nrows_; ++i) {
        mat_[i * (1 + nrows_)] = 1.0;
      }
    }
  }

  using FieldType = ScalarType;
  using VectorType = std::vector<ScalarType>;
  using VectorPtrType = std::vector<ScalarType*>;

  std::size_t nrows() const { return nrows_; }
  std::vector<ScalarType> mat() const { return mat_; }
  std::vector<ScalarType> row(const OrdinalType rowIdx) {
    auto rowBegin = mat_.begin() + rowIdx * nrows_;
    auto rowEnd = mat_.begin() + (rowIdx + 1) * nrows_;
    return std::vector<ScalarType>(rowBegin, rowEnd);
  }

  std::vector<ScalarType*> get_row_ptrs(const OrdinalType rowIdx) {
    std::vector<ScalarType*> ret;
    ret.reserve(nrows_);
    for (OrdinalType colIdx = 0; colIdx < nrows_; ++colIdx) {
      ret.push_back(&((*this)(rowIdx, colIdx)));
    }
    return ret;
  }

  ScalarType& operator()(const std::size_t flatIdx) { return mat_[flatIdx]; }
  ScalarType& operator()(const std::size_t rowIdx, const std::size_t colIdx) {
    return (*this)(rowIdx * nrows_ + colIdx);
  }
  const ScalarType& operator()(const std::size_t flatIdx) const {
    return mat_[flatIdx];
  }
  const ScalarType& operator()(const std::size_t rowIdx,
                               const std::size_t colIdx) const {
    return (*this)(rowIdx * nrows_ + colIdx);
  }

  void print() const {
    std::cout << "mat(" << nrows_ << "x" << nrows_ << ") = [\n";
    for (std::size_t i = 0; i < nrows_; ++i) {
      for (std::size_t j = 0; j < nrows_; ++j) {
        std::cout << (*this)(i, j) << " ";
      }
      std::cout << ((i + 1 != nrows_) ? ',' : ']') << '\n';
    }
  }
};

template <typename ScalarType>
Matrix<ScalarType> operator+(const Matrix<ScalarType>& left,
                             const Matrix<ScalarType>& right) {
  const std::size_t n = left.nrows();
  Matrix<ScalarType> tmp(n, false);
  for (std::size_t i = 0; i < n * n; ++i) {
    tmp(i) = left(i) + right(i);
  }
  return tmp;
}

template <typename ScalarType>
Matrix<ScalarType> operator-(const Matrix<ScalarType>& left,
                             const Matrix<ScalarType>& right) {
  const std::size_t n = left.nrows();

  assert(n == right.nrows());

  Matrix<ScalarType> tmp(n, false);
  for (std::size_t i = 0; i < n * n; ++i) {
    tmp(i) = left(i) - right(i);
  }
  return tmp;
}

template <typename ScalarType>
Matrix<ScalarType> operator*(const Matrix<ScalarType>& left,
                             const Matrix<ScalarType>& right) {
  const std::size_t n = left.nrows();

  assert(n == right.nrows());

  Matrix<ScalarType> tmp(n, false);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        tmp(i, j) += left(i, k) * right(k, j);
      }
    }
  }
  return tmp;
}
}

#endif /* MATRIX_H_ */
