#include "Compounds.h"
#include "DiagonalAccumulator.h"
#include "SparseAccumulator.h"
#include "Transformer.h"
#include <iostream>
#include <map>

using ast::Tag;
using ast::TagComparator;
using ast::TagList;

using matrix::Matrix;
using matrix::OrdinalType;
using matrix::MatCtorType;

using ScalarType = double;
using ScalarMapType = maputil::TagFieldMapType<ScalarType>;
using MatrixType = Matrix<ScalarType>;
using VectorType = typename Matrix<ScalarType>::VectorType;
using VectorPtrType = typename Matrix<ScalarType>::VectorPtrType;

using OrdinalType = matrix::OrdinalType;
using OrdinalPairType = std::pair<OrdinalType, OrdinalType>;
using OrdinalVecType = std::vector<OrdinalType>;
using OrdinalTagMapType = maputil::OrdinalTagMapType<OrdinalType>;
using OrdinalPairTagMapType = maputil::OrdinalPairTagMapType<OrdinalType>;

using DiagAccumType = ast::DiagonalAccumulator<ScalarType>;
using SparseAccumType = ast::SparseAccumulator<ScalarType>;
using TransformType = ast::Transformer<ScalarType>;

#define FIELD(map, name, value) map[Tag(name)] = value;
#define DIAG(map, ordinal, name) map[ordinal] = Tag(name);
#define SPARSE(map, rowidx, colidx, name) \
  map[std::make_pair(rowidx, colidx)] = Tag(name);

int main() {
  // build the map of names (tags) to fields (values)
  ScalarMapType scalarMap;

  FIELD(scalarMap, "a1", 2)
  FIELD(scalarMap, "a2", 2)
  FIELD(scalarMap, "a3", 2)
  FIELD(scalarMap, "a4", 2)

  FIELD(scalarMap, "b1", 3)
  FIELD(scalarMap, "b2", 3)

  FIELD(scalarMap, "c1", 10)
  FIELD(scalarMap, "c2", 10)
  FIELD(scalarMap, "c3", 10)
  FIELD(scalarMap, "c4", 10)

  // build the accumulators
  OrdinalTagMapType dRadVDiagMap;
  OrdinalTagMapType dRbdVDiagMap;
  OrdinalPairTagMapType dRcdUDiagMap;

  DIAG(dRadVDiagMap, 0, "a1")
  DIAG(dRadVDiagMap, 1, "a2")
  DIAG(dRadVDiagMap, 2, "a3")
  DIAG(dRadVDiagMap, 3, "a4")

  DIAG(dRbdVDiagMap, 0, "b1")
  DIAG(dRbdVDiagMap, 2, "b2")

  SPARSE(dRcdUDiagMap, 0, 1, "c1")
  SPARSE(dRcdUDiagMap, 0, 2, "c2")
  SPARSE(dRcdUDiagMap, 0, 3, "c3")

  DiagAccumType dRadV(dRadVDiagMap);
  DiagAccumType dRbdV(dRbdVDiagMap);
  SparseAccumType dRcdU(dRcdUDiagMap);

  // build the transform
  TransformType dVdU;

  // build the Jacobian matrix
  const OrdinalType nrows = 4;
  MatrixType jac1(nrows, MatCtorType::ZEROS);
  MatrixType jac2(nrows, MatCtorType::ZEROS);

  // assemble the Jacobian
  //  auto jacobian1 = dRcdU + (dRadV + dRadV) * dVdU;
  //  auto jacobian2 = (dRadV + dRadV) * dVdU + dRcdU;
  //  auto jacobian1 = dRcdU + (dRadV + dRadV) * dVdU + (dRadV + dRadV) * dVdU;
  //  auto jacobian2 = (dRadV + dRadV) * dVdU + (dRadV + dRadV) * dVdU + dRcdU;

  auto jacobian1 = (dRadV + dRbdV) * dVdU + (dRadV + dRbdV) * dVdU;
  auto jacobian2 = (dRadV + dRbdV) * dVdU + dRcdU;

  for (OrdinalType rowIdx = 0; rowIdx < nrows; ++rowIdx) {
    VectorPtrType row1 = jac1.get_row_ptrs(rowIdx);
    jacobian1.assemble(row1, rowIdx, scalarMap, ast::Place());
    VectorPtrType row2 = jac2.get_row_ptrs(rowIdx);
    jacobian2.assemble(row2, rowIdx, scalarMap, ast::Place());
  }
  std::cout << '\n';
  jac1.print();
  std::cout << '\n';
  jac2.print();
  std::cout << '\n';

  return 0;
}
