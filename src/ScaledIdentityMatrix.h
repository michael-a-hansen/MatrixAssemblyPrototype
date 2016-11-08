/*
 * ScaledIdentityMatrix.h
 *
 *  Created on: Nov 7, 2016
 *      Author: mike
 */

#ifndef SCALEDIDENTITYMATRIX_H_
#define SCALEDIDENTITYMATRIX_H_

#include "AssemblerBase.h"

namespace ast {

enum class ScalingType { SCALEMULTIPLY, SCALEDIVIDE };

/*
 * @class ScaledIdentityMatrix
 *
 * This class represents a scaled identity matrix. It uses either a field
 * obtained through a tag or a double precision number as the scale. And it has
 * a ScalingType object that determines if the scale is a multiplication or
 * division. The scale type defaults to multiplication while the tag defaults to
 * an empty tag, which results in an unscaled identity matrix. This class has
 * all operations defined - Place, AddIn, SubIn, and RMult.
 */
template <typename FieldT>
class ScaledIdentityMatrix : public AssemblerBase<FieldT> {
 protected:
  const Tag scaleTag_;
  const ScalingType scalingType_;
  const bool noScaling_;

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
   * @brief construct a ScaledIdentityMatrix from a scale (tag or double) and a
   * scaling type (divide or multiply)
   * @param scale the double precision scaling factor or field tag used to
   * obtain the scaling field, defaults to an empty tag
   * @param scalingType the method of scaling, multiplication or division,
   * defaults to multiplication
   */
  ScaledIdentityMatrix(
      const Tag scaleTag = Tag(),
      const ScalingType scalingType = ScalingType::SCALEMULTIPLY)
      : AssemblerBase<FieldT>(scaleTag_),
        scaleTag_(scaleTag),
        scalingType_(scalingType),
        noScaling_(scaleTag.is_empty()) {}

  /*
   * @brief assemble this scaled identity matrix by emplacement
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
    if (noScaling_) {
      *row[rowIdx] = 1.0;
    } else {
      const FieldT& s = fieldReqs.at(scaleTag_);
      if (scalingType_ == ScalingType::SCALEMULTIPLY) {
        *row[rowIdx] = s;
      } else {
        *row[rowIdx] = 1.0 / s;
      }
    }
  }

  /*
   * @brief perform addition of this scaled identity matrix
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
    if (noScaling_) {
      *row[rowIdx] += 1.0;
    } else {
      const FieldT& s = fieldReqs.at(scaleTag_);
      if (scalingType_ == ScalingType::SCALEMULTIPLY) {
        *row[rowIdx] += s;
      } else {
        *row[rowIdx] += 1.0 / s;
      }
    }
  }

  /*
   * @brief perform subtraction of this scaled identity matrix
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
    if (noScaling_) {
      *row[rowIdx] -= 1.0;
    } else {
      const FieldT& s = fieldReqs.at(scaleTag_);
      if (scalingType_ == ScalingType::SCALEMULTIPLY) {
        *row[rowIdx] -= s;
      } else {
        *row[rowIdx] -= 1.0 / s;
      }
    }
  }

  /*
   * @brief perform right multiplication by this scaled identity matrix
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
    if (!noScaling_) {
      const FieldT& s = fieldReqs.at(scaleTag_);
      if (scalingType_ == ScalingType::SCALEMULTIPLY) {
        for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
          *row[colIdx] *= s;
        }
      } else {
        for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
          *row[colIdx] *= 1.0 / s;
        }
      }
    }
  }
};
}

#endif /* SCALEDIDENTITYMATRIX_H_ */
