#ifndef TRANSFORMER_H_
#define TRANSFORMER_H_

#include "AssemblerBase.h"
#include <vector>

namespace ast {

/*
 * @brief this method is used by the Transformer class below.
 *
 * This function performs right multiplication by the 4x4 matrix:
 * B = [1 0 0 0,
 *      2 1 0 0,
 *      0 2 1 0,
 *      0 0 2 1].
 *
 * The first row of A*B for a matrix A is then
 * [A11 + 2*A12, A12 + 2*A13, A13 + 2*A14, A14].
 */
template <typename FieldType>
void apply_to_row(std::vector<FieldType*>& row) {
  const std::size_t nrows = row.size();
  std::vector<FieldType> tmp(nrows);
  for (std::size_t i = 0; i < nrows - 1; ++i) {
    tmp[i] = *row[i] + 2 * *row[i + 1];
  }
  tmp[nrows - 1] = *row[nrows - 1];
  for (std::size_t i = 0; i < nrows; ++i) {
    *row[i] = tmp[i];
  }
}

/*
 * @class Transformer
 *
 * This type of assembler is only capable of being right multiplied on another
 * assembler. It cannot be emplaced, added, or subtracted, because the
 * assemble(...,Place), ... methods are not defined. If an emplacement is
 * attempted, the base class (AssemblerBase) assemble(...,Place) method will be
 * called and a runtime error will be thrown.
 *
 * This class is empty and is used as an example of calling an external method
 * to do the right multiplication. Here the apply_to_row function above is used.
 */
template <typename FieldT>
class Transformer : public AssemblerBase<FieldT> {
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
   * @brief construct a Transformer
   */
  Transformer() : AssemblerBase<FieldT>() {}

  using AssemblerBase<FieldT>::assemble;

  /*
   * @brief right multiply by this Transformer
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
    apply_to_row(row);
  }
};
}

#endif /* TRANSFORMER_H_ */
