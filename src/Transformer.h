#ifndef TRANSFORMER_H_
#define TRANSFORMER_H_

#include "AssemblerBase.h"
#include <vector>

namespace ast {

// A * B, B is a sparse matrix never formed
// B = [1 0 0 0,
//      2 1 0 0,
//      0 2 1 0,
//      0 0 2 1],
//
// AB = [A11 + 2*A12, A12 + 2*A13, A13 + 2*A14, A14,
//       A21 + 2*A22, A22 + 2*A23, A23 + 2*A24, A24,
//       A31 + 2*A32, A32 + 2*A33, A33 + 2*A34, A34,
//       A41 + 2*A42, A42 + 2*A43, A43 + 2*A44, A44],
//
// so each row i of A is turned into the following row in AB
//      [Ai1 + 2*Ai2, Ai2 + 2*Ai3, Ai3 + 2*Ai4, Ai4]
//
// we therefore do a scratch row and loop instead of a scratch matrix
//
// apply the sparse transformation to a single row
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

// a Transformer only does right multiplication, and it uses the function above
template <typename FieldT>
class Transformer : public AssemblerBase<FieldT> {
 protected:
 public:
  Transformer() : AssemblerBase<FieldT>() {}

  using FieldType = typename AssemblerBase<FieldT>::FieldType;
  using MatrixType = typename AssemblerBase<FieldT>::MatrixType;
  using VectorType = typename AssemblerBase<FieldT>::VectorType;
  using VectorPtrType = typename AssemblerBase<FieldT>::VectorPtrType;
  using OrdinalType = typename AssemblerBase<FieldT>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<FieldT>::TagFieldMapType;

  using AssemblerBase<FieldT>::assemble;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const RMult) const {
    assert(rowIdx < row.size());
    std::cout << "Transformer: RMult->\n";
    apply_to_row(row);
    std::cout << "Transformer: RMult.\n";
  }
};
}

#endif /* TRANSFORMER_H_ */
