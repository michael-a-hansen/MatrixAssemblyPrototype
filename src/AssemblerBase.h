/*
 * Assembler.h
 *
 *  Created on: Oct 26, 2016
 *      Author: mike
 */

#ifndef AssemblerBase_H_
#define AssemblerBase_H_

#include "MapUtil.h"
#include "Matrix.h"
#include "Operations.h"
#include "Tag.h"
#include <cassert>

namespace ast {

/* @class: AssemblerBase
 *
 * This class is the base for all matrix assemblers - the compositional units in
 * matrix assembly. This base class stores a vector of required tags and the
 * number of those tags (and getters). It defines fundamental types of field,
 * vector, and matrix objects, as well as for ordinals and tag-field maps.
 *
 * This class has a virtual method, assemble, which is defined for an Operation
 * argument of type 'Place.' This method should never be called unless from a
 * derived class that overrides the method, so its default implementation in
 * this base class is to throw an error.
 *
 * All derived classes must be template classes. If they are to be used in an
 * emplacement operation, they must override the assemble(...,Place) method.
 * If they are never to be used in emplacement, for instance a state
 * transformation assembler that is only used for right multiplication, then do
 * not override the assemble(...,Place) method! This will provide compile-time
 * safety on compositions of assemblers.
 *
 * There are four operation types defined between assemblers:
 *
 * 1. Place - emplacing elements into a matrix:
 * \code assemble(..., Place) \endcode
 *
 * 2. AddIn - adding elements into a matrix:
 * \code assemble(..., AddIn) \endcode
 *
 * 3. SubIn - subtracting elements from a matrix:
 * \code assemble(..., SubIn) \endcode
 *
 * 4. RMult - performing right multiplication on a matrix:
 * \code assemble(..., RMult) \endcode
 *
 * These operations must be provided on a row-wise basis, meaning a reference to
 * a matrix row is the input and must be modified in-place to provide the row
 * after the operation. To assemble a matrix, the Place operation may be used on
 * a pointer/reference of AssemblerBase type. As mentioned above, if this method
 * is not overridden by a derived class, then the code will throw a runtime
 * error.
 *
 * If a derived assembler is to be used only for right multiplication, for
 * instance, then define only the assemble(...,RMult) method to get compile-time
 * checking that the assembler is only ever used in that mode.
 *
 */
template <typename FieldT>
class AssemblerBase {
 protected:
  const TagList requiredTags_;
  const std::size_t nTags_;

 public:
  AssemblerBase() : nTags_(0) {}
  AssemblerBase(const TagList requiredTags)
      : requiredTags_(requiredTags), nTags_(requiredTags.size()) {}

  virtual ~AssemblerBase() {}

  using FieldType = FieldT;
  using VectorType = typename matrix::Matrix<FieldType>::VectorType;
  using VectorPtrType = typename matrix::Matrix<FieldType>::VectorPtrType;
  using MatrixType = typename matrix::Matrix<FieldType>;
  using OrdinalType = typename matrix::OrdinalType;
  using TagFieldMapType = typename maputil::TagFieldMapType<FieldType>;

  /*
   * @brief get the list of required tags for this assembler object
   */
  TagList required_tags() const { return requiredTags_; }

  /*
   * @brief get the number of required tags for this assembler object
   */
  std::size_t n_tags() const { return nTags_; }

  /*
   * @brief assemble into a matrix row by emplacement
   * @param a vector of field pointers representing a matrix row
   * @param the index of the row being assembled
   * @param the field requests needed by the assembler or by all composed
   * assembler objects
   * @param the mode of operation, in this case it must be a Place() object
   * (literally just use \code Place() \endcode here)
   *
   * In this class' documentation it is discussed that AssemblerBase will throw
   * a runtime error if its assemble(...,Place) method is called without being
   * overridden.
   *
   */
  virtual void assemble(VectorPtrType& row,
                        const OrdinalType rowIdx,
                        const TagFieldMapType& fieldReqs,
                        const Place) const {
    throw std::runtime_error(
        "AssemblerBase::assemble(...,Place) was called. The assembler you "
        "called assemble(...,Place) does not override assemble(...,Place), "
        "which it must do. Throwing an error...");
  }
};
}

#endif /* AssemblerBase_H_ */
