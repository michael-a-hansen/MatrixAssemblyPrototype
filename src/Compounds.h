#ifndef COMPOUNDS_H_
#define COMPOUNDS_H_

#include "AssemblerBase.h"
#include <cassert>

namespace ast {

/*
 * @class RowPlaceHolder
 *
 * This class serves as a place-holder of a dense row for compound operations
 * that must respect a certain order, such as a product of an assembler and a
 * summation. The summation must be evaluated row-wise first, and a
 * RowPlaceHolder is used to store the result. A RowPlaceHolder can do
 * emplacement, addition, and subtraction. It cannot be used to RMult, because
 * it can not know anything about columns. It stores a vector of field pointers.
 */
template <typename FieldT>
class RowPlaceHolder : public AssemblerBase<FieldT> {
 protected:
  const std::vector<FieldT*>& tmp_;

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
   * @brief construct a RowPlaceHolder from a vector of field pointers.
   */
  RowPlaceHolder(const std::vector<FieldT*>& tmp)
      : AssemblerBase<FieldT>(), tmp_(tmp) {}

  /*
   * @brief emplace this row
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
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      *row[colIdx] = *tmp_[colIdx];
    }
  }

  /*
   * @brief add this row
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
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      *row[colIdx] += *tmp_[colIdx];
    }
  }

  /*
   * @brief subtract this row
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
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
    for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      *row[colIdx] -= *tmp_[colIdx];
    }
  }
};

/*
 * @class Summation
 *
 * This compound assembler class is a summation of assembler types (possibly
 * compounds themselves). This class supports emplacement, addition, and
 * subtraction. Right-multiplication by a summation is enabled by the ProdSum
 * compound assembler below, which behaves in a special manner.
 *
 * Through the Summation and ProdSum compounds, a summation may be used in
 * emplacement, addition, subtraction, and right multiplication.
 */
template <typename LeftAssemblerT, typename RightAssemblerT>
class Summation : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const RightAssemblerT right_;

 public:
  /*
   * types from AssemblerBase templated on the left assembler type
   */
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

  /*
   * @brief construct a summation object from two assemblers. This should not be
   * called from user code. Summations should be made only with the '+'
   * operator overloaded below.
   * @param left the assembler object on the left of '+'
   * @param right the assembler object on the right of '+'
   */
  Summation(const LeftAssemblerT& left, const RightAssemblerT& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(std::is_same<typename LeftAssemblerT::FieldType,
                               typename RightAssemblerT::FieldType>::value,
                  "different field types in Summation - error!");
  }

  /*
   * @brief assemble this assembler
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
    left_.assemble(row, rowIdx, fieldReqs, Place());
    right_.assemble(row, rowIdx, fieldReqs, AddIn());
  }

  /*
   * @brief add this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an AddIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    left_.assemble(row, rowIdx, fieldReqs, AddIn());
    right_.assemble(row, rowIdx, fieldReqs, AddIn());
  }

  /*
   * @brief subtract this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a SubIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    left_.assemble(row, rowIdx, fieldReqs, SubIn());
    right_.assemble(row, rowIdx, fieldReqs, SubIn());
  }
};

/*
 * @brief operator '+' overloaded to produce a summation object
 * @param left the assembler object on the left of '+'
 * @param right the assembler object on the right of '+'
 */
template <typename LeftAssemblerT, typename RightAssemblerT>
Summation<LeftAssemblerT, RightAssemblerT> operator+(LeftAssemblerT left,
                                                     RightAssemblerT right) {
  return Summation<LeftAssemblerT, RightAssemblerT>(left, right);
}

/*
 * @class ProdSum
 *
 * This compound assembler class is a product of an assembler type and a
 * summation type.
 *
 * Through the Summation and ProdSum compounds, a summation may be used in
 * emplacement, addition, subtraction, and right multiplication.
 */
template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
class ProdSum : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const Summation<RightLeftAssemblerT, RightRightAssemblerT> right_;

 public:
  /*
   * types from AssemblerBase
   */
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

  /*
   * @brief construct a product-summation object from an assembler on the left
   * and a summation assembler on the right
   * operator overloaded below.
   * @param left the assembler object on the left of '*'
   * @param right the summation assembler on the right of '*'
   */
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

  /*
   * @brief assemble this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a Place() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                Place) const {
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
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
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
      }
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, Place());
  }

  /*
   * @brief add this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an AddIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
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
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
      }
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, AddIn());
  }

  /*
   * @brief subtract this row
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a SubIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
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
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
      }
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
      }
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, SubIn());
  }
};

/*
 * @brief operator '*' overloaded to produce a ProdSum object
 * @param left the assembler object on the left of '*'
 * @param right the summation assembler on the right of '*'
 */
template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
ProdSum<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT> operator*(
    LeftAssemblerT left,
    Summation<RightLeftAssemblerT, RightRightAssemblerT> right) {
  return ProdSum<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT>(
      left, right);
}

/*
 * @class Difference
 *
 * This compound assembler class is a difference of assembler types (possibly
 * compounds themselves). This class supports emplacement, addition, and
 * subtraction. Right-multiplication by a summation is enabled by the ProdDiff
 * compound assembler below, which behaves in a special manner.
 *
 * Through the Difference and ProdDiff compounds, a difference may be used in
 * emplacement, addition, subtraction, and right multiplication.
 */
template <typename LeftAssemblerT, typename RightAssemblerT>
class Difference : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const RightAssemblerT right_;

 public:
  /*
   * types from AssemblerBase
   */
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

  /*
   * @brief construct a difference object from two assemblers. This should not
   * be called from user code. Differences should be made only with the '-'
   * operator overloaded below.
   * @param left the assembler object on the left of '-'
   * @param right the assembler object on the right of '-'
   */
  Difference(const LeftAssemblerT& left, const RightAssemblerT& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(std::is_same<typename LeftAssemblerT::FieldType,
                               typename RightAssemblerT::FieldType>::value,
                  "different field types in Difference - error!");
  }

  /*
   * @brief assemble this assembler
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
    left_.assemble(row, rowIdx, fieldReqs, Place());
    right_.assemble(row, rowIdx, fieldReqs, SubIn());
  }

  /*
   * @brief add this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an AddIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
    left_.assemble(row, rowIdx, fieldReqs, SubIn());
    right_.assemble(row, rowIdx, fieldReqs, SubIn());
  }

  /*
   * @brief subtract this row
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a SubIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
    left_.assemble(row, rowIdx, fieldReqs, SubIn());
    right_.assemble(row, rowIdx, fieldReqs, AddIn());
  }
};

/*
 * @brief operator '-' overloaded to produce a difference object
 * @param left the assembler object on the left of '-'
 * @param right the assembler object on the right of '-'
 */
template <typename LeftAssemblerT, typename RightAssemblerT>
Difference<LeftAssemblerT, RightAssemblerT> operator-(LeftAssemblerT left,
                                                      RightAssemblerT right) {
  return Difference<LeftAssemblerT, RightAssemblerT>(left, right);
}

/*
 * @class ProdDiff
 *
 * This compound assembler class is a product of an assembler type and a
 * difference type.
 *
 * Through the Difference and ProdDiff compounds, a difference may be used in
 * emplacement, addition, subtraction, and right multiplication.
 */
template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
class ProdDiff : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const Difference<RightLeftAssemblerT, RightRightAssemblerT> right_;

 public:
  /*
   * types from AssemblerBase
   */
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

  /*
   * @brief construct a product-difference object from an assembler on the left
   * and a difference assembler on the right
   * operator overloaded below.
   * @param left the assembler object on the left of '*'
   * @param right the difference assembler on the right of '*'
   */
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
        "different field types in ProdDiff - error!");
  }

  /*
   * @brief assemble this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a Place() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                Place) const {
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
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
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
      }
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, Place());
  }

  /*
   * @brief add this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an AddIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
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
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
      }
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, AddIn());
  }

  /*
   * @brief subtract this row
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a SubIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    const OrdinalType nrows = row.size();
    assert(rowIdx < nrows);
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
    for (OrdinalType rowIdxLooped = 0; rowIdxLooped < nrows; ++rowIdxLooped) {
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *tmpPtrs[colIdx] = 0;
      }
      right_.assemble(tmpPtrs, rowIdxLooped, fieldReqs, Place());
      for (OrdinalType colIdx = 0; colIdx < nrows; ++colIdx) {
        *resPtrs[colIdx] += *lftPtrs[rowIdxLooped] * *tmpPtrs[colIdx];
      }
    }
    RowPlaceHolder<FieldType> tempHolder(resPtrs);
    tempHolder.assemble(row, rowIdx, fieldReqs, SubIn());
  }
};

/*
 * @brief operator '*' overloaded to produce a ProdDiff object
 * @param left the assembler object on the left of '*'
 * @param right the difference assembler on the right of '*'
 */
template <typename LeftAssemblerT,
          typename RightLeftAssemblerT,
          typename RightRightAssemblerT>
ProdDiff<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT> operator*(
    LeftAssemblerT left,
    Difference<RightLeftAssemblerT, RightRightAssemblerT> right) {
  return ProdDiff<LeftAssemblerT, RightLeftAssemblerT, RightRightAssemblerT>(
      left, right);
}

/*
 * @class Product
 *
 * This compound assembler class is a product of assembler types (possibly
 * compounds themselves). This class supports emplacement, addition, and
 * subtraction.
 */
template <typename LeftAssemblerT, typename RightAssemblerT>
class Product : public AssemblerBase<typename LeftAssemblerT::FieldType> {
 protected:
  const LeftAssemblerT left_;
  const RightAssemblerT right_;

 public:
  /*
   * types from AssemblerBase
   */
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

  /*
   * @brief construct a product object from two assemblers. This should not
   * be called from user code. Differences should be made only with the '*'
   * operator overloaded below.
   * @param left the assembler object on the left of '*'
   * @param right the assembler object on the right of '*'
   */
  Product(const LeftAssemblerT& left, const RightAssemblerT& right)
      : AssemblerBase<typename LeftAssemblerT::FieldType>(
            union_taglists(left.required_tags(), right.required_tags())),
        left_(left),
        right_(right) {
    static_assert(std::is_same<typename LeftAssemblerT::FieldType,
                               typename RightAssemblerT::FieldType>::value,
                  "different field types in Product - error!");
  }

  /*
   * @brief assemble this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a Place() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                Place) const override {
    assert(rowIdx < row.size());
    left_.assemble(row, rowIdx, fieldReqs, Place());
    right_.assemble(row, rowIdx, fieldReqs, RMult());
  }

  /*
   * @brief multiply by this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an RMult() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                RMult) const {
    assert(rowIdx < row.size());
    left_.assemble(row, rowIdx, fieldReqs, RMult());
    right_.assemble(row, rowIdx, fieldReqs, RMult());
  }

  /*
   * @brief add this assembler
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - an AddIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                AddIn) const {
    assert(rowIdx < row.size());
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
  }

  /*
   * @brief subtract this row
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation - a SubIn() object here
   */
  void assemble(VectorPtrType& row,
                const OrdinalType rowIdx,
                const TagFieldMapType& fieldReqs,
                SubIn) const {
    assert(rowIdx < row.size());
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
  }
};

/*
 * @brief operator '*' overloaded to produce a product object
 * @param left the assembler object on the left of '*'
 * @param right the assembler object on the right of '*'
 */
template <typename LeftAssemblerT, typename RightAssemblerT>
Product<LeftAssemblerT, RightAssemblerT> operator*(LeftAssemblerT left,
                                                   RightAssemblerT right) {
  return Product<LeftAssemblerT, RightAssemblerT>(left, right);
}
}

#endif /* COMPOUNDS_H_ */
