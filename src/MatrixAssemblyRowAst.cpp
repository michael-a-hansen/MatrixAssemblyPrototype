#include "Compounds.h"
#include "Transformer.h"
#include "DiagonalMatrix.h"
#include "SparseMatrix.h"
#include "ScaledIdentityMatrix.h"
#include "TestHelper.h"
#include <iostream>
#include <map>

using ast::Tag;
using ast::TagComparator;
using ast::TagList;

using matrix::Matrix;
using matrix::OrdinalType;
using matrix::MatCtorType;

using FieldType = double;
using FieldMapType = maputil::TagFieldMapType<FieldType>;
using MatrixType = Matrix<FieldType>;
using VectorType = typename Matrix<FieldType>::VectorType;
using VectorPtrType = typename Matrix<FieldType>::VectorPtrType;

using OrdinalType = matrix::OrdinalType;
using OrdinalPairType = std::pair<OrdinalType, OrdinalType>;
using OrdinalVecType = std::vector<OrdinalType>;
using OrdinalTagMapType = maputil::OrdinalTagMapType<OrdinalType>;
using OrdinalPairTagMapType = maputil::OrdinalPairTagMapType<OrdinalType>;

using DiagMatrixType = ast::DiagonalMatrix<FieldType>;
using SparseMatrixType = ast::SparseMatrix<FieldType>;
using ScaledIdentityType = ast::ScaledIdentityMatrix<FieldType>;
using TransformType = ast::Transformer<FieldType>;

int main() {
  FieldMapType fieldMap;

  Tag aTag("a"), twoATag("2a"), threeATag("3a"), minusATag("-a");
  FieldType aValue = 7.0;
  FieldType twoAValue = aValue * 2.0;
  FieldType threeAValue = aValue * 3.0;
  FieldType minusAValue = -aValue;
  fieldMap[aTag] = aValue;
  fieldMap[twoATag] = twoAValue;
  fieldMap[threeATag] = threeAValue;
  fieldMap[minusATag] = minusAValue;

  FIELD(fieldMap, "-1", -1.0);
  FIELD(fieldMap, "0", 0.0);
  FIELD(fieldMap, "2", 2.0);
  FIELD(fieldMap, "3", 3.0);

  FIELD(fieldMap, "D1", 1.0)
  FIELD(fieldMap, "D2", 2.0)
  FIELD(fieldMap, "D3", 3.0)
  FIELD(fieldMap, "D4", 4.0)
  FIELD(fieldMap, "-D1", -1.0)
  FIELD(fieldMap, "-D2", -2.0)
  FIELD(fieldMap, "-D3", -3.0)
  FIELD(fieldMap, "-D4", -4.0)
  FIELD(fieldMap, "D1*2", 2.0)
  FIELD(fieldMap, "D2*2", 4.0)
  FIELD(fieldMap, "D3*2", 6.0)
  FIELD(fieldMap, "D4*2", 8.0)
  FIELD(fieldMap, "D1*3", 3.0)
  FIELD(fieldMap, "D2*3", 6.0)
  FIELD(fieldMap, "D3*3", 9.0)
  FIELD(fieldMap, "D4*3", 12.0)

  FIELD(fieldMap, "M11", 1.0)
  FIELD(fieldMap, "M12", 2.0)
  FIELD(fieldMap, "M13", 3.0)
  FIELD(fieldMap, "M14", 4.0)
  FIELD(fieldMap, "M21", 5.0)
  FIELD(fieldMap, "M22", 6.0)
  FIELD(fieldMap, "M23", 7.0)
  FIELD(fieldMap, "M24", 8.0)
  FIELD(fieldMap, "M31", 9.0)
  FIELD(fieldMap, "M32", 10.0)
  FIELD(fieldMap, "M33", 11.0)
  FIELD(fieldMap, "M34", 12.0)
  FIELD(fieldMap, "M41", 13.0)
  FIELD(fieldMap, "M42", 14.0)
  FIELD(fieldMap, "M43", 15.0)
  FIELD(fieldMap, "M44", 16.0)
  FIELD(fieldMap, "M11*2", 2.0 * 1.0)
  FIELD(fieldMap, "M12*2", 2.0 * 2.0)
  FIELD(fieldMap, "M13*2", 2.0 * 3.0)
  FIELD(fieldMap, "M14*2", 2.0 * 4.0)
  FIELD(fieldMap, "M21*2", 2.0 * 5.0)
  FIELD(fieldMap, "M22*2", 2.0 * 6.0)
  FIELD(fieldMap, "M23*2", 2.0 * 7.0)
  FIELD(fieldMap, "M24*2", 2.0 * 8.0)
  FIELD(fieldMap, "M31*2", 2.0 * 9.0)
  FIELD(fieldMap, "M32*2", 2.0 * 10.0)
  FIELD(fieldMap, "M33*2", 2.0 * 11.0)
  FIELD(fieldMap, "M34*2", 2.0 * 12.0)
  FIELD(fieldMap, "M41*2", 2.0 * 13.0)
  FIELD(fieldMap, "M42*2", 2.0 * 14.0)
  FIELD(fieldMap, "M43*2", 2.0 * 15.0)
  FIELD(fieldMap, "M44*2", 2.0 * 16.0)
  FIELD(fieldMap, "M11*3", 3.0 * 1.0)
  FIELD(fieldMap, "M12*3", 3.0 * 2.0)
  FIELD(fieldMap, "M13*3", 3.0 * 3.0)
  FIELD(fieldMap, "M14*3", 3.0 * 4.0)
  FIELD(fieldMap, "M21*3", 3.0 * 5.0)
  FIELD(fieldMap, "M22*3", 3.0 * 6.0)
  FIELD(fieldMap, "M23*3", 3.0 * 7.0)
  FIELD(fieldMap, "M24*3", 3.0 * 8.0)
  FIELD(fieldMap, "M31*3", 3.0 * 9.0)
  FIELD(fieldMap, "M32*3", 3.0 * 10.0)
  FIELD(fieldMap, "M33*3", 3.0 * 11.0)
  FIELD(fieldMap, "M34*3", 3.0 * 12.0)
  FIELD(fieldMap, "M41*3", 3.0 * 13.0)
  FIELD(fieldMap, "M42*3", 3.0 * 14.0)
  FIELD(fieldMap, "M43*3", 3.0 * 15.0)
  FIELD(fieldMap, "M44*3", 3.0 * 16.0)
  FIELD(fieldMap, "-M11", -1.0)
  FIELD(fieldMap, "-M12", -2.0)
  FIELD(fieldMap, "-M13", -3.0)
  FIELD(fieldMap, "-M14", -4.0)
  FIELD(fieldMap, "-M21", -5.0)
  FIELD(fieldMap, "-M22", -6.0)
  FIELD(fieldMap, "-M23", -7.0)
  FIELD(fieldMap, "-M24", -8.0)
  FIELD(fieldMap, "-M31", -9.0)
  FIELD(fieldMap, "-M32", -10.0)
  FIELD(fieldMap, "-M33", -11.0)
  FIELD(fieldMap, "-M34", -12.0)
  FIELD(fieldMap, "-M41", -13.0)
  FIELD(fieldMap, "-M42", -14.0)
  FIELD(fieldMap, "-M43", -15.0)
  FIELD(fieldMap, "-M44", -16.0)

  // build matrices
  ScaledIdentityType I;
  ScaledIdentityType twoI(Tag("2"));
  ScaledIdentityType threeI(Tag("3"));
  ScaledIdentityType mI(Tag("-1"));
  ScaledIdentityType Z(Tag("0"));

  ScaledIdentityType S(aTag);
  ScaledIdentityType twoS(twoATag);
  ScaledIdentityType threeS(threeATag);
  ScaledIdentityType minusS(minusATag);

  DiagMatrixType D;
  D(0) = Tag("D1");
  D(1) = Tag("D2");
  D(2) = Tag("D3");
  D(3) = Tag("D4");
  D.is_completely_entered();
  DiagMatrixType twoD;
  twoD(0) = Tag("D1*2");
  twoD(1) = Tag("D2*2");
  twoD(2) = Tag("D3*2");
  twoD(3) = Tag("D4*2");
  twoD.is_completely_entered();
  DiagMatrixType threeD;
  threeD(0) = Tag("D1*3");
  threeD(1) = Tag("D2*3");
  threeD(2) = Tag("D3*3");
  threeD(3) = Tag("D4*3");
  threeD.is_completely_entered();
  DiagMatrixType minusD;
  minusD(0) = Tag("-D1");
  minusD(1) = Tag("-D2");
  minusD(2) = Tag("-D3");
  minusD(3) = Tag("-D4");
  minusD.is_completely_entered();

  SparseMatrixType M;
  M(0, 0) = Tag("M11");
  M(0, 1) = Tag("M12");
  M(0, 2) = Tag("M13");
  M(0, 3) = Tag("M14");
  M(1, 0) = Tag("M21");
  M(1, 1) = Tag("M22");
  M(1, 2) = Tag("M23");
  M(1, 3) = Tag("M24");
  M(2, 0) = Tag("M31");
  M(2, 1) = Tag("M32");
  M(2, 2) = Tag("M33");
  M(2, 3) = Tag("M34");
  M(3, 0) = Tag("M41");
  M(3, 1) = Tag("M42");
  M(3, 2) = Tag("M43");
  M(3, 3) = Tag("M44");
  M.is_completely_entered();
  SparseMatrixType twoM;
  twoM(0, 0) = Tag("M11*2");
  twoM(0, 1) = Tag("M12*2");
  twoM(0, 2) = Tag("M13*2");
  twoM(0, 3) = Tag("M14*2");
  twoM(1, 0) = Tag("M21*2");
  twoM(1, 1) = Tag("M22*2");
  twoM(1, 2) = Tag("M23*2");
  twoM(1, 3) = Tag("M24*2");
  twoM(2, 0) = Tag("M31*2");
  twoM(2, 1) = Tag("M32*2");
  twoM(2, 2) = Tag("M33*2");
  twoM(2, 3) = Tag("M34*2");
  twoM(3, 0) = Tag("M41*2");
  twoM(3, 1) = Tag("M42*2");
  twoM(3, 2) = Tag("M43*2");
  twoM(3, 3) = Tag("M44*2");
  twoM.is_completely_entered();
  SparseMatrixType threeM;
  threeM(0, 0) = Tag("M11*3");
  threeM(0, 1) = Tag("M12*3");
  threeM(0, 2) = Tag("M13*3");
  threeM(0, 3) = Tag("M14*3");
  threeM(1, 0) = Tag("M21*3");
  threeM(1, 1) = Tag("M22*3");
  threeM(1, 2) = Tag("M23*3");
  threeM(1, 3) = Tag("M24*3");
  threeM(2, 0) = Tag("M31*3");
  threeM(2, 1) = Tag("M32*3");
  threeM(2, 2) = Tag("M33*3");
  threeM(2, 3) = Tag("M34*3");
  threeM(3, 0) = Tag("M41*3");
  threeM(3, 1) = Tag("M42*3");
  threeM(3, 2) = Tag("M43*3");
  threeM(3, 3) = Tag("M44*3");
  threeM.is_completely_entered();
  SparseMatrixType minusM;
  minusM(0, 0) = Tag("-M11");
  minusM(0, 1) = Tag("-M12");
  minusM(0, 2) = Tag("-M13");
  minusM(0, 3) = Tag("-M14");
  minusM(1, 0) = Tag("-M21");
  minusM(1, 1) = Tag("-M22");
  minusM(1, 2) = Tag("-M23");
  minusM(1, 3) = Tag("-M24");
  minusM(2, 0) = Tag("-M31");
  minusM(2, 1) = Tag("-M32");
  minusM(2, 2) = Tag("-M33");
  minusM(2, 3) = Tag("-M34");
  minusM(3, 0) = Tag("-M41");
  minusM(3, 1) = Tag("-M42");
  minusM(3, 2) = Tag("-M43");
  minusM(3, 3) = Tag("-M44");
  minusM.is_completely_entered();

  // run tests on compositions
  const OrdinalType nrows = 4;

  ast::TestHelper tester;
  ast::TestDisplay testDisplay = ast::TestDisplay::PASSFAIL;

  std::cout << "- IdentityMatrix\n";
  tester.test(I + I, twoI, fieldMap, nrows, "AddIn x 1", testDisplay);
  tester.test(I + I + I, threeI, fieldMap, nrows, "AddIn x 2", testDisplay);
  tester.test(I - I, Z, fieldMap, nrows, "SubIn x 1", testDisplay);
  tester.test(I - I - I, mI, fieldMap, nrows, "SubIn x 2", testDisplay);
  tester.test(I * I, I, fieldMap, nrows, "RMult", testDisplay);
  tester.test(I * I * I, I, fieldMap, nrows, "RMult x 2", testDisplay);
  std::cout << "--------------------\n\n";

  std::cout << "- ScaledIdentityMatrix\n";
  tester.test(I * S, S, fieldMap, nrows, "RMult", testDisplay);
  tester.test(I * S * S, S * S, fieldMap, nrows, "RMult x 2", testDisplay);
  tester.test(S + S, twoS, fieldMap, nrows, "AddIn", testDisplay);
  tester.test(S + S + S, threeS, fieldMap, nrows, "AddIn x 2", testDisplay);
  tester.test(S - S, Z, fieldMap, nrows, "SubIn", testDisplay);
  tester.test(S - S - S, minusS, fieldMap, nrows, "SubIn x 2", testDisplay);
  std::cout << "--------------------\n\n";

  std::cout << "- DiagonalMatrix\n";
  tester.test(I * D, D, fieldMap, nrows, "RMult", testDisplay);
  tester.test(I * D * D, D * D, fieldMap, nrows, "RMult x 2", testDisplay);
  tester.test(D + D, twoD, fieldMap, nrows, "AddIn", testDisplay);
  tester.test(D + D + D, threeD, fieldMap, nrows, "AddIn x 2", testDisplay);
  tester.test(D - D, Z, fieldMap, nrows, "SubIn", testDisplay);
  tester.test(D - D - D, minusD, fieldMap, nrows, "SubIn x 2", testDisplay);
  std::cout << "--------------------\n\n";

  std::cout << "- SparseMatrix\n";
  tester.test(I * M, M, fieldMap, nrows, "RMult", testDisplay);
  tester.test(I * M * M, M * M, fieldMap, nrows, "RMult x 2", testDisplay);
  tester.test(M + M, twoM, fieldMap, nrows, "AddIn", testDisplay);
  tester.test(M + M + M, threeM, fieldMap, nrows, "AddIn x 2", testDisplay);
  tester.test(M - M, Z, fieldMap, nrows, "SubIn", testDisplay);
  tester.test(M - M - M, minusM, fieldMap, nrows, "SubIn x 2", testDisplay);
  std::cout << "--------------------\n\n";

  std::cout << "- ProdSum, ProdDiff\n";
  tester.test(I * (I + I), twoI, fieldMap, nrows, "RMult by Sum", testDisplay);
  tester.test(I * (I + I + I), threeI, fieldMap, nrows, "RMult by Sum x2",
              testDisplay);
  tester.test(I * (I + I - I), I, fieldMap, nrows, "RMult by Sum+Diff",
              testDisplay);
  tester.test(I * (I - I), Z, fieldMap, nrows, "RMult by Diff", testDisplay);
  tester.test(I * (I - I - I), mI, fieldMap, nrows, "RMult by Diff x2",
              testDisplay);
  tester.test(I * (I - I + I), I, fieldMap, nrows, "RMult by Diff+Sum",
              testDisplay);
  std::cout << "--------------------\n\n";

  tester.report();

  return 0;
}
