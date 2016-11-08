#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <iostream>
#include <cassert>

namespace matrix {
using OrdinalType = std::size_t;

/*
 * @enum MatCtorType
 * @brief matrix construction enumeration, to faciliate building all-zero or
 * identity matrices
 */
enum class MatCtorType { ZEROS, IDENTITY };

/*
 * @class Matrix
 *
 * This class represents a 2-D square array and provides indexing methods for a
 * row-major storage in a flat vector.
 */
template <typename ScalarType>
class Matrix {
 protected:
  OrdinalType nrows_;
  std::vector<ScalarType> mat_;

 public:
  /*
   * types for classes that use matrices
   */
  using FieldType = ScalarType;
  using VectorType = std::vector<ScalarType>;
  using VectorPtrType = std::vector<ScalarType*>;

  /*
   * @brief construct an empty matrix
   */
  Matrix() : nrows_(0) {}

  /*
   * @brief construct a matrix
   * @param nrows number of rows and columns in the matrix
   * @param c style of constructor, default is MatCtorType::ZEROS
   */
  Matrix(const OrdinalType nrows, const MatCtorType c = MatCtorType::ZEROS)
      : nrows_(nrows), mat_(nrows * nrows, 0) {
    if (c == MatCtorType::IDENTITY) {
      for (std::size_t i = 0; i < nrows_; ++i) {
        mat_[i * (1 + nrows_)] = 1.0;
      }
    }
  }

  /*
   * @brief get the number of rows in this matrix
   */
  OrdinalType nrows() const { return nrows_; }

  /*
   * @brief get the std::vector object that stores the dense matrix
   */
  std::vector<ScalarType> mat() const { return mat_; }

  /*
   * @brief get the row of a certain index
   * @param rowIdx the index of the row to be obtained
   */
  std::vector<ScalarType> row(const OrdinalType rowIdx) {
    auto rowBegin = mat_.begin() + rowIdx * nrows_;
    auto rowEnd = mat_.begin() + (rowIdx + 1) * nrows_;
    return std::vector<ScalarType>(rowBegin, rowEnd);
  }

  /*
   * @brief get a vector of pointers to fields of a certain row
   * @param rowIdx the index of the row whose pointers are returned
   */
  std::vector<ScalarType*> get_row_ptrs(const OrdinalType rowIdx) {
    std::vector<ScalarType*> ret;
    ret.reserve(nrows_);
    for (OrdinalType colIdx = 0; colIdx < nrows_; ++colIdx) {
      ret.push_back(&((*this)(rowIdx, colIdx)));
    }
    return ret;
  }

  /*
   * @brief non-const flat-indexing (row-major) of this matrix through the
   * parentheses operator
   * @param flatIdx the row-major index
   */
  ScalarType& operator()(const OrdinalType flatIdx) { return mat_[flatIdx]; }

  /*
   * @brief non-const 2-D indexing of this matrix through the parentheses
   * operator
   * @param rowIdx the row index
   * @param colIdx the column index
   */
  ScalarType& operator()(const OrdinalType rowIdx, const OrdinalType colIdx) {
    return (*this)(rowIdx * nrows_ + colIdx);
  }

  /*
   * @brief const flat-indexing (row-major) of this matrix through the
   * parentheses operator
   * @param flatIdx the row-major index
   */
  const ScalarType& operator()(const OrdinalType flatIdx) const {
    return mat_[flatIdx];
  }

  /*
   * @brief const 2-D indexing of this matrix through the parentheses
   * operator
   * @param rowIdx the row index
   * @param colIdx the column index
   */
  const ScalarType& operator()(const OrdinalType rowIdx,
                               const OrdinalType colIdx) const {
    return (*this)(rowIdx * nrows_ + colIdx);
  }

  /*
   * @brief matrix equality operator
   * @param mat matrix to compare
   */
  bool operator==(const Matrix& mat) {
    if (mat.nrows() != nrows_) {
      return false;
    }
    assert(mat.nrows() == nrows_);
    for (OrdinalType rowIdx = 0; rowIdx < nrows_; ++rowIdx) {
      for (OrdinalType colIdx = 0; colIdx < nrows_; ++colIdx) {
        if (mat(rowIdx, colIdx) != (*this)(rowIdx, colIdx)) {
          return false;
        }
      }
    }
    return true;
  }

  /*
   * @brief print the matrix to std::cout
   */
  void print() const {
    std::cout << "mat(" << nrows_ << "x" << nrows_ << ") = [\n";
    for (OrdinalType i = 0; i < nrows_; ++i) {
      for (OrdinalType j = 0; j < nrows_; ++j) {
        std::cout << (*this)(i, j) << " ";
      }
      std::cout << ((i + 1 != nrows_) ? ',' : ']') << '\n';
    }
  }
};
}

#endif /* MATRIX_H_ */
