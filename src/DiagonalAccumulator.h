#ifndef DIAGONALACCUMULATOR_H_
#define DIAGONALACCUMULATOR_H_

#include "AssemblerBase.h"
#include "MapUtil.h"

namespace ast {

using OrdinalTagMapType = maputil::OrdinalTagMapType<matrix::OrdinalType>;

template <typename FieldT>
class DiagonalAccumulator : public AssemblerBase<FieldT> {
 protected:
  const OrdinalTagMapType& diagMap_;

 public:
  DiagonalAccumulator(const OrdinalTagMapType& diagMap)
      : AssemblerBase<FieldT>(maputil::extract_values(diagMap)),
        diagMap_(diagMap) {}

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
    std::cout << "DiagAccumulator: Place->\n";
    assert(AssemblerBase<FieldT>::nTags_ <= row.size());

    if (diagMap_.find(rowIdx) != diagMap_.end()) {
      *row[rowIdx] = fieldReqs.at(diagMap_.at(rowIdx));
    }
    std::cout << "DiagAccumulator: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "DiagAccumulator: AddIn->\n";
    assert(AssemblerBase<FieldT>::nTags_ <= row.size());

    if (diagMap_.find(rowIdx) != diagMap_.end()) {
      *row[rowIdx] += fieldReqs.at(diagMap_.at(rowIdx));
    }
    std::cout << "DiagAccumulator: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "DiagAccumulator: SubIn->\n";
    assert(AssemblerBase<FieldT>::nTags_ <= row.size());

    if (diagMap_.find(rowIdx) != diagMap_.end()) {
      *row[rowIdx] -= fieldReqs.at(diagMap_.at(rowIdx));
    }
    std::cout << "DiagAccumulator: SubIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const RMult) const {
    assert(rowIdx < row.size());
    std::cout << "DiagAccumulator: RMult->\n";
    const OrdinalType nrows = row.size();
    assert(AssemblerBase<FieldT>::nTags_ <= row.size());

    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      if (diagMap_.find(colIdx) != diagMap_.end()) {
        *row[colIdx] *= fieldReqs.at(diagMap_.at(colIdx));
      }
    }
    std::cout << "DiagAccumulator: RMult.\n";
  }
};
}

#endif /* DIAGONALACCUMULATOR_H_ */
