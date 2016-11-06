#ifndef SPARSEACCUMULATOR_H_
#define SPARSEACCUMULATOR_H_

#include "AssemblerBase.h"
#include "MapUtil.h"

namespace ast {
using ElementMapType = maputil::OrdinalPairTagMapType<matrix::OrdinalType>;

template <typename FieldT>
class SparseAccumulator : public AssemblerBase<FieldT> {
 protected:
  const ElementMapType& elementMap_;

 public:
  SparseAccumulator(const ElementMapType& elementMap)
      : AssemblerBase<FieldT>(maputil::extract_values(elementMap)),
        elementMap_(elementMap) {}

  using FieldType = typename AssemblerBase<FieldT>::FieldType;
  using MatrixType = typename AssemblerBase<FieldT>::MatrixType;
  using VectorType = typename AssemblerBase<FieldT>::VectorType;
  using VectorPtrType = typename AssemblerBase<FieldT>::VectorPtrType;
  using OrdinalType = typename AssemblerBase<FieldT>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<FieldT>::TagFieldMapType;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const Place) const override {
    std::cout << "SparseAccumulator: Place->\n";
    assert(rowIdx < row.size());
    const OrdinalType nrows = row.size();
    assert(AssemblerBase<FieldT>::nTags_ <= nrows);

    for (OrdinalType colIdx = 0; colIdx < nrows; colIdx++) {
      auto pair = std::make_pair(rowIdx, colIdx);
      if (elementMap_.find(pair) != elementMap_.end()) {
        *row[colIdx] = fieldReqs.at(elementMap_.at(pair));
      }
    }
    std::cout << "SparseAccumulator: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const AddIn) const {
    std::cout << "SparseAccumulator: AddIn->\n";
    assert(rowIdx < row.size());
    const OrdinalType nrows = row.size();
    assert(AssemblerBase<FieldT>::nTags_ <= nrows);

    for (OrdinalType colIdx = 0; colIdx < nrows; colIdx++) {
      auto pair = std::make_pair(rowIdx, colIdx);
      if (elementMap_.find(pair) != elementMap_.end()) {
        *row[colIdx] += fieldReqs.at(elementMap_.at(pair));
      }
    }
    std::cout << "SparseAccumulator: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const SubIn) const {
    std::cout << "SparseAccumulator: SubIn->\n";
    assert(rowIdx < row.size());
    const OrdinalType nrows = row.size();
    assert(AssemblerBase<FieldT>::nTags_ <= nrows);

    for (OrdinalType colIdx = 0; colIdx < nrows; colIdx++) {
      auto pair = std::make_pair(rowIdx, colIdx);
      if (elementMap_.find(pair) != elementMap_.end()) {
        *row[colIdx] -= fieldReqs.at(elementMap_.at(pair));
      }
    }
    std::cout << "SparseAccumulator: SubIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const RMult) const {
    std::cout << "SparseAccumulator: RMult->\n";
    const OrdinalType nrows = row.size();
    assert(AssemblerBase<FieldT>::nTags_ <= nrows);

    for (OrdinalType colIdxToAssign = 0; colIdxToAssign < nrows;
         colIdxToAssign++) {
      FieldT temp = 0.0;
      for (OrdinalType colIdx = 0; colIdx < nrows; colIdx++) {
        auto pair = std::make_pair(colIdx, colIdxToAssign);
        if (elementMap_.find(pair) != elementMap_.end()) {
          temp += *row[colIdx] * fieldReqs.at(elementMap_.at(pair));
        }
      }
      *row[colIdxToAssign] = temp;
    }
    std::cout << "SparseAccumulator: RMult.\n";
  }
};
}

#endif /* SPARSEACCUMULATOR_H_ */
