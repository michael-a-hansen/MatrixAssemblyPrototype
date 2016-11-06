#ifndef MAPUTIL_H_
#define MAPUTIL_H_

#include "Tag.h"
#include <map>

namespace maputil {

template <typename KeyType,
          typename ValueType,
          typename ComparatorType = std::less<KeyType>>
using BaseMapType = std::map<KeyType, ValueType, ComparatorType>;

template <typename FieldType>
using TagFieldMapType = BaseMapType<ast::Tag, FieldType, ast::TagComparator>;

template <typename OrdinalType>
using OrdinalTagMapType = BaseMapType<OrdinalType, ast::Tag>;

template <typename OrdinalType>
using OrdinalPairTagMapType =
    BaseMapType<std::pair<OrdinalType, OrdinalType>, ast::Tag>;

template <typename KeyType, typename ValueType, typename ComparatorType>
std::vector<KeyType> extract_keys(
    const BaseMapType<KeyType, ValueType, ComparatorType>& input_map) {
  std::vector<KeyType> retval;
  for (auto const& element : input_map) {
    retval.push_back(element.first);
  }
  return retval;
}

template <typename KeyType, typename ValueType, typename ComparatorType>
std::vector<ValueType> extract_values(
    const BaseMapType<KeyType, ValueType, ComparatorType>& input_map) {
  std::vector<ValueType> retval;
  for (auto const& element : input_map) {
    retval.push_back(element.second);
  }
  return retval;
}
}

#endif /* MAPUTIL_H_ */
