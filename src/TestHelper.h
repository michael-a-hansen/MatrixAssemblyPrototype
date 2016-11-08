/*
 * TestHelper.h
 *
 *  Created on: Nov 8, 2016
 *      Author: mike
 */

#ifndef TESTHELPER_H_
#define TESTHELPER_H_

namespace ast {

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

#define FIELD(map, name, value) map[Tag(name)] = value;

/*
 * @class TestDisplay
 *
 * A strong enumeration used to determine what information about a test is
 * displayed to std::cout upon test completion.
 *
 * PASSFAIL: show if the test passes or fails
 * FAILONLY: only show failing tests
 * VERBOSE: PASSFAIL plus show the assembled matrices
 * NONE: show nothing
 */
enum class TestDisplay { PASSFAIL, FAILONLY, VERBOSE, NONE };

/*
 * @class TestHelper
 *
 * This class facilitates testing of assemblers.
 */
class TestHelper {
 protected:
  int numTest_;
  int numPass_;

 public:
  /*
   * @brief build a TestHelper (empty)
   */
  TestHelper() : numTest_(0), numPass_(0) {}

  /*
   * @brief get the number of tests this helper has run
   */
  int num_test() { return numTest_; }
  /*
   * @brief get the number of passing tests this helper has run
   */
  int num_pass() { return numPass_; }

  /*
   * @brief write number of passing / number in total
   */
  void report() {
    std::cout << "-- passing " << numPass_ << " of " << numTest_ << " tests\n";
  }

  /*
   * @brief test two assemblers for equality, returns a bool if equal
   * @param testAssember the assembler of interest
   * @param answerAssembler the reference assembler, test should match this
   * @param FieldMap the tag-field map used by assemblers for field access
   * @param nrows the size of the matrix
   * @param testName the name of the test mentioned in pass/fail statement
   * @param disp the type of information displayed
   */
  template <typename TestAssemblerFieldT, typename ResultAssemblerFieldT>
  bool test(const ast::AssemblerBase<TestAssemblerFieldT>& testAssembler,
            const ast::AssemblerBase<ResultAssemblerFieldT>& answerAssembler,
            const FieldMapType& FieldMap,
            const int nrows,
            const std::string testName,
            const TestDisplay disp = TestDisplay::PASSFAIL) {
    MatrixType test(nrows, MatCtorType::ZEROS);
    MatrixType answer(nrows, MatCtorType::ZEROS);

    for (OrdinalType rowIdx = 0; rowIdx < nrows; ++rowIdx) {
      VectorPtrType testrow = test.get_row_ptrs(rowIdx);
      testAssembler.assemble(testrow, rowIdx, FieldMap, ast::Place());
      VectorPtrType ansrow = answer.get_row_ptrs(rowIdx);
      answerAssembler.assemble(ansrow, rowIdx, FieldMap, ast::Place());
    }

    const bool pass = (test == answer);

    if (disp == TestDisplay::PASSFAIL) {
      if (pass) {
        std::cout << "pass - " << testName << '\n';
      } else {
        std::cout << "FAIL - " << testName << '\n';
      }
    } else if (disp == TestDisplay::FAILONLY) {
      if (!pass)
        std::cout << "FAIL - " << testName << '\n';
    } else if (disp == TestDisplay::VERBOSE) {
      if (pass) {
        std::cout << "pass - " << testName << '\n';
      } else {
        std::cout << "FAIL - " << testName << '\n';
      }
      if (!pass) {
        std::cout << "Answer:\n";
        answer.print();
        std::cout << "Calculated:\n";
        test.print();
      }
    }

    numTest_++;
    if (pass) {
      numPass_++;
    }

    return pass;
  }
};
}

#endif /* TESTHELPER_H_ */
