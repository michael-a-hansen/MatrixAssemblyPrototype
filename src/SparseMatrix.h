#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#include "AssemblerBase.h"
#include "MapUtil.h"

namespace ast {

/*
 * @class SpareMatrix
 *
 * This class is a generic sparse matrix that is capable of all four operations:
 * emplacement, addition-in, subtraction-from, and right multiplication.
 *
 * A SparseMatrix is constructed without arguments, and elements should be added
 * with the parentheses operator, as M(i,j) = aTag; where i is the row index, j
 * is the column index, and aTag is the field tag that goes in element (i,j).
 *
 * The is_completely_entered() function MUST BE CALLED before a SparseMatrix can
 * be used in assembly.
 *
 * An example construction and filling procedure is:
 *
 * \code
 * SparseMatrix<FieldT> M;
 * M(0,0) = aTag;
 * M(1,0) = bTag;
 * M(1,1) = cTag;
 * M.is_completely_entered();
 * \endcode
 */
template <typename FieldT>
class SparseMatrix : public AssemblerBase<FieldT> {
 public:
  using OrdinalType = typename AssemblerBase<FieldT>::OrdinalType;

 protected:
  using ElementMapType = maputil::OrdinalPairTagMapType<matrix::OrdinalType>;
  using OrdinalMapType = maputil::OrdinalTagMapType<matrix::OrdinalType>;
  using RowColMapType = std::map<OrdinalType, OrdinalMapType>;

  bool isCompletelyEntered_ = false;
  ElementMapType elementMap_;
  RowColMapType activeRowMaps_;
  RowColMapType activeColMaps_;

 public:
  /*
   * types from AssemblerBase
   */
  using FieldType = typename AssemblerBase<FieldT>::FieldType;
  using MatrixType = typename AssemblerBase<FieldT>::MatrixType;
  using VectorType = typename AssemblerBase<FieldT>::VectorType;
  using VectorPtrType = typename AssemblerBase<FieldT>::VectorPtrType;
  using TagFieldMapType = typename AssemblerBase<FieldT>::TagFieldMapType;

  /*
   * @brief construct an empty SparseMatrix
   */
  SparseMatrix() : AssemblerBase<FieldT>() {}

  /*
   * @brief Obtain a reference to the field tag associated with a certain
   * element, identified by row and column indices. Use this function to insert
   * tags into the matrix, as M(i,j) = aTag. If an i,j pair is already
   * associated with a tag, it will be overwritten by assignment.
   * @param rowIdx row index
   * @param colIdx column index
   */
  Tag& operator()(const OrdinalType rowIdx, const OrdinalType colIdx) {
    if (isCompletelyEntered_) {
      throw std::runtime_error(
          "SparseMatrix is fully entered! Can't add to it!");
    }
    std::pair<OrdinalType, OrdinalType> ij = std::make_pair(rowIdx, colIdx);
    auto iter = elementMap_.find(ij);
    if (iter == elementMap_.end()) {
      elementMap_[ij] = Tag();
    }
    return elementMap_[ij];
  }

  /*
   * @brief finalize the matrix and prevent further entries from being added,
   * build faster column and row mappings, and set required tags
   */
  void is_completely_entered() {
    isCompletelyEntered_ = true;
    activeRowMaps_ = maputil::rowify_element_map(elementMap_);
    activeColMaps_ = maputil::colify_element_map(elementMap_);
    AssemblerBase<FieldT>::set_required_tags(
        maputil::extract_values(elementMap_));
  }

  /*
   * @brief assemble this sparse matrix
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a Place() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const Place) const override {
    assert(rowIdx < row.size());
    if (!isCompletelyEntered_) {
      throw std::runtime_error("SparseMatrix was not fully entered!");
    }
    if (activeRowMaps_.find(rowIdx) != activeRowMaps_.end()) {
      OrdinalMapType columnMap = activeRowMaps_.at(rowIdx);
      for (auto c : columnMap) {
        *row[c.first] = fieldReqs.at(c.second);
      }
    }
  }

  /*
   * @brief add this sparse matrix
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an AddIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const AddIn) const {
    assert(rowIdx < row.size());
    if (!isCompletelyEntered_) {
      throw std::runtime_error("SparseMatrix was not fully entered!");
    }
    if (activeRowMaps_.find(rowIdx) != activeRowMaps_.end()) {
      OrdinalMapType columnMap = activeRowMaps_.at(rowIdx);
      for (auto c : columnMap) {
        *row[c.first] += fieldReqs.at(c.second);
      }
    }
  }

  /*
   * @brief subtract this sparse matrix
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a SubIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const SubIn) const {
    assert(rowIdx < row.size());
    if (!isCompletelyEntered_) {
      throw std::runtime_error("SparseMatrix was not fully entered!");
    }
    if (activeRowMaps_.find(rowIdx) != activeRowMaps_.end()) {
      OrdinalMapType columnMap = activeRowMaps_.at(rowIdx);
      for (auto c : columnMap) {
        *row[c.first] -= fieldReqs.at(c.second);
      }
    }
  }

  /*
   * @brief multiply by this sparse matrix
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an RMult() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const RMult) const {
    assert(rowIdx < row.size());
    if (!isCompletelyEntered_) {
      throw std::runtime_error("SparseMatrix was not fully entered!");
    }
    const OrdinalType nrows = row.size();
    VectorType leftMatRow(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      leftMatRow[colIdx] = *row[colIdx];
    }
    for (OrdinalType colIdxToAssign = 0; colIdxToAssign < nrows;
         colIdxToAssign++) {
      FieldT temp = 0.0;
      if (activeColMaps_.find(colIdxToAssign) != activeColMaps_.end()) {
        OrdinalMapType rowMap = activeColMaps_.at(colIdxToAssign);
        for (auto r : rowMap) {
          temp += leftMatRow[r.first] * fieldReqs.at(r.second);
        }
      }
      *row[colIdxToAssign] = temp;
    }
  }
};
}

#endif /* SPARSEMATRIX_H_ */
