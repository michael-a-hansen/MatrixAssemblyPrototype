#ifndef MAPUTIL_H_
#define MAPUTIL_H_

#include "Tag.h"
#include <map>
#include <set>

namespace maputil {

/*
 * @typedef BaseMapType
 *
 * This is the map type used throughout the entire code, defined here to
 * (possibly) facilitate changing this type later on.
 */
template <typename KeyType,
          typename ValueType,
          typename ComparatorType = std::less<KeyType>>
using BaseMapType = std::map<KeyType, ValueType, ComparatorType>;

/*
 * @typedef TagFieldMapType
 *
 * A map between Tags and Fields
 */
template <typename FieldType>
using TagFieldMapType = BaseMapType<ast::Tag, FieldType, ast::TagComparator>;

/*
 * @typedef OrdinalTagMapType
 *
 * A map between ordinals and Tags (e.g., for diagonal matrix)
 */
template <typename OrdinalType>
using OrdinalTagMapType = BaseMapType<OrdinalType, ast::Tag>;

/*
 * @typedef OrdinalPairTagMapType
 *
 * A map between ordinal pairs and Tags (e.g., for sparse matrix)
 */
template <typename OrdinalType>
using OrdinalPairTagMapType =
    BaseMapType<std::pair<OrdinalType, OrdinalType>, ast::Tag>;

/*
 * @brief extract the keys of a map into a vector
 */
template <typename KeyType, typename ValueType, typename ComparatorType>
std::vector<KeyType> extract_keys(
    const BaseMapType<KeyType, ValueType, ComparatorType>& input_map) {
  std::vector<KeyType> retval;
  for (auto const& element : input_map) {
    retval.push_back(element.first);
  }
  return retval;
}

/*
 * @brief extract the values of a map into a vector
 */
template <typename KeyType, typename ValueType, typename ComparatorType>
std::vector<ValueType> extract_values(
    const BaseMapType<KeyType, ValueType, ComparatorType>& input_map) {
  std::vector<ValueType> retval;
  for (auto const& element : input_map) {
    retval.push_back(element.second);
  }
  return retval;
}

/*
 * @brief convert a map<pair<ordinal1,ordinal2>, tag> to
 * map<ordinal1, map<ordinal2,tag>> so that the left ordinal may be iterated.
 *
 * 'rowify' means taking a map between (i,j) pairs and tags, and converting it
 * to a map from between row index and a column-tag map. This way we can
 * iterate over columns of active rows one at a time without having to consider
 * the entire map and do comparisons of pairs.
 *
 * This is useful for emplacement, addition, and subtraction of matrix rows.
 */
template <typename OrdinalType, typename FieldType>
BaseMapType<OrdinalType, BaseMapType<OrdinalType, FieldType>>
rowify_element_map(const BaseMapType<std::pair<OrdinalType, OrdinalType>,
                                     FieldType>& elementMap) {
  std::set<OrdinalType> activeRows;
  std::map<OrdinalType, OrdinalTagMapType<OrdinalType>> activeRowMaps;

  for (auto m : elementMap) {
    activeRows.insert(m.first.first);
  }
  for (auto row : activeRows) {
    OrdinalTagMapType<OrdinalType> columnMap;
    for (auto m : elementMap) {
      if (m.first.first == row) {
        columnMap[m.first.second] = m.second;
      }
    }
    activeRowMaps[row] = columnMap;
  }
  return activeRowMaps;
}

/*
 * @brief convert a map<pair<ordinal1,ordinal2>, tag> to
 * map<ordinal2, map<ordinal1,tag>> so that the right ordinal may be iterated.
 *
 * 'colify' means taking a map between (i,j) pairs and tags, and converting it
 * to a map from between column index and a row-tag map. This way we can
 * iterate over rows of active columns one at a time without having to consider
 * the entire map and do comparisons of pairs.
 *
 * This is useful for right-multiplication when we need fast access to columns.
 */
template <typename OrdinalType, typename FieldType>
BaseMapType<OrdinalType, BaseMapType<OrdinalType, FieldType>>
colify_element_map(const BaseMapType<std::pair<OrdinalType, OrdinalType>,
                                     FieldType>& elementMap) {
  std::set<OrdinalType> activeCols;
  std::map<OrdinalType, OrdinalTagMapType<OrdinalType>> activeColMaps;

  for (auto m : elementMap) {
    activeCols.insert(m.first.second);
  }
  for (auto col : activeCols) {
    OrdinalTagMapType<OrdinalType> rowMap;
    for (auto m : elementMap) {
      if (m.first.second == col) {
        rowMap[m.first.first] = m.second;
      }
    }
    activeColMaps[col] = rowMap;
  }
  return activeColMaps;
}
}

#endif /* MAPUTIL_H_ */
