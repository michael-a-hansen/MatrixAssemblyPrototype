#ifndef DIAGONALMATRIX_H_
#define DIAGONALMATRIX_H_

#include "AssemblerBase.h"
#include "MapUtil.h"

namespace ast {

/*
 * @class DiagonalMatrix
 *
 * This class represents a diagonal matrix. It maps from ordinals (diagonal
 * indices) to tags of fields representing the diagonal elements. This class has
 * all operations defined - Place, AddIn, SubIn, and RMult.
 *
 * A DiagonalMatrix is constructed without arguments, and elements should be
 * added with the parentheses operator, as M(i) = aTag; where i is the diagonal
 * index, and aTag is the field tag that goes in element (i,i).
 *
 * The is_completely_entered() function MUST BE CALLED before a DiagonalMatrix
 * can be used in assembly.
 *
 * An example construction and filling procedure is:
 *
 * \code
 * DiagonalMatrix<FieldT> D;
 * D(0,0) = aTag;
 * D(1,1) = bTag;
 * D.is_completely_entered();
 * \endcode
 */
template <typename FieldT>
class DiagonalMatrix : public AssemblerBase<FieldT> {
 protected:
  using OrdinalTagMapType = maputil::OrdinalTagMapType<matrix::OrdinalType>;

  bool isCompletelyEntered_ = false;
  OrdinalTagMapType diagMap_;

 public:
  /*
   * types from AssemblerBase
   */
  using FieldType = typename AssemblerBase<FieldT>::FieldType;
  using MatrixType = typename AssemblerBase<FieldT>::MatrixType;
  using VectorType = typename AssemblerBase<FieldT>::VectorType;
  using VectorPtrType = typename AssemblerBase<FieldT>::VectorPtrType;
  using OrdinalType = typename AssemblerBase<FieldT>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<FieldT>::TagFieldMapType;

  /*
   * @brief construct an empty DiagonalMatrix
   */
  DiagonalMatrix() {}

  /*
   * @brief Obtain a reference to the field tag associated with a certain
   * element, identified by an index along the diagonal. Use this to insert tags
   * into the diagonal, as D(i) = aTag, where i is the diagonal index and aTag
   * represents the field that goes into D(i,i). If a tag is already associated
   * with a diagonal index passed as an argument here, then it will be
   * overwritten by assignment.
   * @param diagIdx the diagonal index
   */
  Tag& operator()(const OrdinalType diagIdx) {
    if (isCompletelyEntered_) {
      throw std::runtime_error(
          "DiagonalMatrix is fully entered! Can't add to it!");
    }
    auto iter = diagMap_.find(diagIdx);
    if (iter == diagMap_.end()) {
      diagMap_[diagIdx] = Tag();
    }
    return diagMap_[diagIdx];
  }

  /*
   * @brief finalize the matrix and prevent further entries from being added,
   * and set required tags
   */
  void is_completely_entered() {
    isCompletelyEntered_ = true;
    AssemblerBase<FieldT>::set_required_tags(maputil::extract_values(diagMap_));
  }

  /*
   * @brief assemble this diagonal matrix by emplacement
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
    if (diagMap_.find(rowIdx) != diagMap_.end()) {
      *row[rowIdx] = fieldReqs.at(diagMap_.at(rowIdx));
    }
  }

  /*
   * @brief perform addition of this diagonal matrix
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
    if (diagMap_.find(rowIdx) != diagMap_.end()) {
      *row[rowIdx] += fieldReqs.at(diagMap_.at(rowIdx));
    }
  }

  /*
   * @brief perform subtraction of this diagonal matrix
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
    if (diagMap_.find(rowIdx) != diagMap_.end()) {
      *row[rowIdx] -= fieldReqs.at(diagMap_.at(rowIdx));
    }
  }

  /*
   * @brief perform right multiplication by this diagonal matrix
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
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      if (diagMap_.find(colIdx) != diagMap_.end()) {
        *row[colIdx] *= fieldReqs.at(diagMap_.at(colIdx));
      }
    }
  }
};
}

#endif /* DIAGONALMATRIX_H_ */
