#ifndef COMPOUNDS_H_
#define COMPOUNDS_H_

#include "AssemblerBase.h"
#include <cassert>

namespace ast {
template <typename FieldT>
class RowPlaceHolder : public AssemblerBase<FieldT> {
 protected:
  const std::vector<FieldT*>& tmp_;

 public:
  RowPlaceHolder(const std::vector<FieldT*>& tmp)
      : AssemblerBase<FieldT>(), tmp_(tmp) {}
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
    assert(rowIdx < row.size());
    std::cout << "TEMPORARY ROW: Place->\n";
    const OrdinalType nrows = row.size();
    for (OrdinalType i = 0; i < nrows; ++i) {
      *row[i] = *tmp_[i];
    }
    std::cout << "TEMPORARY ROW: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "TEMPORARY ROW: AddIn->\n";
    const OrdinalType nrows = row.size();
    for (OrdinalType i = 0; i < nrows; ++i) {
      *row[i] += *tmp_[i];
    }
    std::cout << "TEMPORARY ROW: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "TEMPORARY ROW: SubIn->\n";
    const OrdinalType nrows = row.size();
    for (OrdinalType i = 0; i < nrows; ++i) {
      *row[i] -= *tmp_[i];
    }
    std::cout << "TEMPORARY ROW: SubIn.\n";
  }
};

template <typename LeftAssemblerT, typename RightAssemblerT>
class Summation : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const RightAssemblerT right_;

 public:
  Summation(const LeftAssemblerT& left, const RightAssemblerT& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(std::is_same<typename LeftAssemblerT::FieldType,
                               typename RightAssemblerT::FieldType>::value,
                  "different field types in Summation - error!");
  }

  using FieldType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::FieldType;
  using MatrixType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::MatrixType;
  using VectorType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorType;
  using VectorPtrType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorPtrType;
  using OrdinalType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<
      typename LeftAssemblerT::FieldType>::TagFieldMapType;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const Place) const override {
    assert(rowIdx < row.size());
    std::cout << "Summation: Place->\n";
    left_.assemble(row, rowIdx, fieldReqs, Place());
    right_.assemble(row, rowIdx, fieldReqs, AddIn());
    std::cout << "Summation: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "Summation: AddIn->\n";
    left_.assemble(row, rowIdx, fieldReqs, AddIn());
    right_.assemble(row, rowIdx, fieldReqs, AddIn());
    std::cout << "Summation: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "Summation: SubIn->\n";
    left_.assemble(row, rowIdx, fieldReqs, SubIn());
    right_.assemble(row, rowIdx, fieldReqs, SubIn());
    std::cout << "Summation: SubIn.\n";
  }
};

template <typename LeftAssemblerT, typename RightAssemblerT>
Summation<LeftAssemblerT, RightAssemblerT> operator+(LeftAssemblerT left_,
                                                     RightAssemblerT right_) {
  return Summation<LeftAssemblerT, RightAssemblerT>(left_, right_);
}

template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
class ProdSum : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const Summation<RightLeftAssemblerT, RightRightAssemblerT> right_;

 public:
  ProdSum(const LeftAssemblerT& left,
          const Summation<RightLeftAssemblerT, RightRightAssemblerT>& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(
        std::is_same<
            typename LeftAssemblerT::FieldType,
            typename Summation<RightLeftAssemblerT,
                               RightRightAssemblerT>::FieldType>::value,
        "different field types in ProdSum - error!");
  }

  using FieldType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::FieldType;
  using MatrixType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::MatrixType;
  using VectorType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorType;
  using VectorPtrType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorPtrType;
  using OrdinalType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<
      typename LeftAssemblerT::FieldType>::TagFieldMapType;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                Place) const {
    assert(rowIdx < row.size());
    std::cout << "ProdSum: Place->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    std::vector<FieldType> lft(nrows, 0.0);
    std::vector<FieldType> res(nrows, 0.0);
    std::vector<FieldType *> tmpPtrs(nrows), resPtrs(nrows), lftPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
      resPtrs[colIdx] = &res[colIdx];
      lftPtrs[colIdx] = &lft[colIdx];
    }
    left_.assemble(lftPtrs, rowIdx, fieldReqs, Place());
    std::cout << "rowIdx = " << rowIdx << '\n';
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      std::cout << "lftPtrs[" << colIdx << "] = " << *lftPtrs[colIdx] << ", ";
    }
    std::cout << '\n';
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        std::cout << "tmpPtrs[" << colIdx << "] = " << *tmpPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
        std::cout << "resPtrs[" << colIdx << "] = " << *resPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, Place());
    std::cout << "ProdSum: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "ProdSum: AddIn->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    std::vector<FieldType> lft(nrows, 0.0);
    std::vector<FieldType> res(nrows, 0.0);
    std::vector<FieldType *> tmpPtrs(nrows), resPtrs(nrows), lftPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
      resPtrs[colIdx] = &res[colIdx];
      lftPtrs[colIdx] = &lft[colIdx];
    }
    left_.assemble(lftPtrs, rowIdx, fieldReqs, Place());
    std::cout << "rowIdx = " << rowIdx << '\n';
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      std::cout << "lftPtrs[" << colIdx << "] = " << *lftPtrs[colIdx] << ", ";
    }
    std::cout << '\n';
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        std::cout << "tmpPtrs[" << colIdx << "] = " << *tmpPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
        std::cout << "resPtrs[" << colIdx << "] = " << *resPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, AddIn());
    std::cout << "ProdSum: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "ProdSum: SubIn->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    std::vector<FieldType> lft(nrows, 0.0);
    std::vector<FieldType> res(nrows, 0.0);
    std::vector<FieldType *> tmpPtrs(nrows), resPtrs(nrows), lftPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
      resPtrs[colIdx] = &res[colIdx];
      lftPtrs[colIdx] = &lft[colIdx];
    }
    left_.assemble(lftPtrs, rowIdx, fieldReqs, Place());
    std::cout << "rowIdx = " << rowIdx << '\n';
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      std::cout << "lftPtrs[" << colIdx << "] = " << *lftPtrs[colIdx] << ", ";
    }
    std::cout << '\n';
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        std::cout << "tmpPtrs[" << colIdx << "] = " << *tmpPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
        std::cout << "resPtrs[" << colIdx << "] = " << *resPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, SubIn());
    std::cout << "ProdSum: SubIn.\n";
  }
};

template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
ProdSum<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT> operator*(
    LeftAssemblerT left_,
    Summation<RightLeftAssemblerT, RightRightAssemblerT> right_) {
  return ProdSum<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT>(
      left_, right_);
}

template <typename LeftAssemblerT, typename RightAssemblerT>
class Difference : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const RightAssemblerT right_;

 public:
  Difference(const LeftAssemblerT& left, const RightAssemblerT& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(std::is_same<typename LeftAssemblerT::FieldType,
                               typename RightAssemblerT::FieldType>::value,
                  "different field types in Summation - error!");
  }

  using FieldType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::FieldType;
  using MatrixType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::MatrixType;
  using VectorType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorType;
  using VectorPtrType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorPtrType;
  using OrdinalType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<
      typename LeftAssemblerT::FieldType>::TagFieldMapType;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                const Place) const override {
    assert(rowIdx < row.size());
    std::cout << "Difference: Place->\n";
    left_.assemble(row, rowIdx, fieldReqs, Place());
    right_.assemble(row, rowIdx, fieldReqs, SubIn());
    std::cout << "Difference: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "Difference: AddIn->\n";
    left_.assemble(row, rowIdx, fieldReqs, SubIn());
    right_.assemble(row, rowIdx, fieldReqs, SubIn());
    std::cout << "Difference: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "Difference: SubIn->\n";
    left_.assemble(row, rowIdx, fieldReqs, SubIn());
    right_.assemble(row, rowIdx, fieldReqs, AddIn());
    std::cout << "Difference: SubIn.\n";
  }
};

template <typename LeftAssemblerT, typename RightAssemblerT>
Difference<LeftAssemblerT, RightAssemblerT> operator-(LeftAssemblerT left_,
                                                      RightAssemblerT right_) {
  return Difference<LeftAssemblerT, RightAssemblerT>(left_, right_);
}

template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
class ProdDiff : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const Difference<RightLeftAssemblerT, RightRightAssemblerT> right_;

 public:
  ProdDiff(const LeftAssemblerT& left,
           const Difference<RightLeftAssemblerT, RightRightAssemblerT>& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(
        std::is_same<
            typename LeftAssemblerT::FieldType,
            typename Difference<RightLeftAssemblerT,
                                RightRightAssemblerT>::FieldType>::value,
        "different field types in ProdSum - error!");
  }

  using FieldType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::FieldType;
  using MatrixType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::MatrixType;
  using VectorType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorType;
  using VectorPtrType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorPtrType;
  using OrdinalType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<
      typename LeftAssemblerT::FieldType>::TagFieldMapType;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                Place) const {
    assert(rowIdx < row.size());
    std::cout << "ProdDiff: Place->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    std::vector<FieldType> lft(nrows, 0.0);
    std::vector<FieldType> res(nrows, 0.0);
    std::vector<FieldType *> tmpPtrs(nrows), resPtrs(nrows), lftPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
      resPtrs[colIdx] = &res[colIdx];
      lftPtrs[colIdx] = &lft[colIdx];
    }
    left_.assemble(lftPtrs, rowIdx, fieldReqs, Place());
    std::cout << "rowIdx = " << rowIdx << '\n';
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      std::cout << "lftPtrs[" << colIdx << "] = " << *lftPtrs[colIdx] << ", ";
    }
    std::cout << '\n';
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        std::cout << "tmpPtrs[" << colIdx << "] = " << *tmpPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
        std::cout << "resPtrs[" << colIdx << "] = " << *resPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, Place());
    std::cout << "ProdDiff: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "ProdDiff: AddIn->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    std::vector<FieldType> lft(nrows, 0.0);
    std::vector<FieldType> res(nrows, 0.0);
    std::vector<FieldType *> tmpPtrs(nrows), resPtrs(nrows), lftPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
      resPtrs[colIdx] = &res[colIdx];
      lftPtrs[colIdx] = &lft[colIdx];
    }
    left_.assemble(lftPtrs, rowIdx, fieldReqs, Place());
    std::cout << "rowIdx = " << rowIdx << '\n';
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      std::cout << "lftPtrs[" << colIdx << "] = " << *lftPtrs[colIdx] << ", ";
    }
    std::cout << '\n';
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        std::cout << "tmpPtrs[" << colIdx << "] = " << *tmpPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
        std::cout << "resPtrs[" << colIdx << "] = " << *resPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, AddIn());
    std::cout << "ProdDiff: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "ProdDiff: SubIn->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    std::vector<FieldType> lft(nrows, 0.0);
    std::vector<FieldType> res(nrows, 0.0);
    std::vector<FieldType *> tmpPtrs(nrows), resPtrs(nrows), lftPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
      resPtrs[colIdx] = &res[colIdx];
      lftPtrs[colIdx] = &lft[colIdx];
    }
    left_.assemble(lftPtrs, rowIdx, fieldReqs, Place());
    std::cout << "rowIdx = " << rowIdx << '\n';
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      std::cout << "lftPtrs[" << colIdx << "] = " << *lftPtrs[colIdx] << ", ";
    }
    std::cout << '\n';
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        std::cout << "tmpPtrs[" << colIdx << "] = " << *tmpPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
        std::cout << "resPtrs[" << colIdx << "] = " << *resPtrs[colIdx] << ", ";
      }
      std::cout << '\n';
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, SubIn());
    std::cout << "ProdDiff: SubIn.\n";
  }
};

template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
ProdDiff<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT> operator*(
    LeftAssemblerT left_,
    Difference<RightLeftAssemblerT, RightRightAssemblerT> right_) {
  return ProdDiff<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT>(
      left_, right_);
}

template <typename LeftAssemblerT, typename RightAssemblerT>
class Product : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const RightAssemblerT right_;

 public:
  Product(const LeftAssemblerT& left, const RightAssemblerT& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(std::is_same<typename LeftAssemblerT::FieldType,
                               typename RightAssemblerT::FieldType>::value,
                  "different field types in Product - error!");
  }

  using FieldType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::FieldType;
  using MatrixType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::MatrixType;
  using VectorType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorType;
  using VectorPtrType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::VectorPtrType;
  using OrdinalType =
      typename AssemblerBase<typename LeftAssemblerT::FieldType>::OrdinalType;
  using TagFieldMapType = typename AssemblerBase<
      typename LeftAssemblerT::FieldType>::TagFieldMapType;

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                Place) const override {
    assert(rowIdx < row.size());
    std::cout << "Product: Place->\n";
    left_.assemble(row, rowIdx, fieldReqs, Place());
    right_.assemble(row, rowIdx, fieldReqs, RMult());
    std::cout << "Product: Place.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                RMult) const {
    assert(rowIdx < row.size());
    std::cout << "Product: RMult->\n";
    left_.assemble(row, rowIdx, fieldReqs, RMult());
    right_.assemble(row, rowIdx, fieldReqs, RMult());
    std::cout << "Product: RMult.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    std::cout << "Product: AddIn->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    tmp[rowIdx] = 1.0;
    std::vector<FieldType*> tmpPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
    }
    left_.assemble(tmpPtrs, rowIdx, fieldReqs, Place());
    right_.assemble(tmpPtrs, rowIdx, fieldReqs, RMult());
    RowPlaceHolder<FieldType> tempHolder(std::move(tmpPtrs));
    tempHolder.assemble(row, rowIdx, fieldReqs, AddIn());
    std::cout << "Product: AddIn.\n";
  }

  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    std::cout << "Product: SubIn->\n";
    const OrdinalType nrows = row.size();
    std::vector<FieldType> tmp(nrows, 0.0);
    tmp[rowIdx] = 1.0;
    std::vector<FieldType*> tmpPtrs(nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      tmpPtrs[colIdx] = &tmp[colIdx];
    }
    left_.assemble(tmpPtrs, rowIdx, fieldReqs, Place());
    right_.assemble(tmpPtrs, rowIdx, fieldReqs, RMult());
    RowPlaceHolder<FieldType> tempHolder(std::move(tmpPtrs));
    tempHolder.assemble(row, rowIdx, fieldReqs, SubIn());
    std::cout << "Product: SubIn.\n";
  }
};

template <typename LeftAssemblerT, typename RightAssemblerT>
Product<LeftAssemblerT, RightAssemblerT> operator*(LeftAssemblerT left_,
                                                   RightAssemblerT right_) {
  return Product<LeftAssemblerT, RightAssemblerT>(left_, right_);
}
}

#endif /* COMPOUNDS_H_ */
