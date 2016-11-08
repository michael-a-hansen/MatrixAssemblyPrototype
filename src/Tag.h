#ifndef TAG_H_
#define TAG_H_

#include <string>
#include <vector>

namespace ast {

/*
 * @class Tag
 *
 * A simple wrapper for a string meant to mimic the behavior of destination
 * code. A Tag is used to identify fields.
 */
class Tag {
 protected:
  std::string name_;
  bool isEmpty_;

 public:
  /*
   * @brief construct an empty Tag
   */
  Tag() : name_(""), isEmpty_(true) {}

  /*
   * @brief construct a non-empty Tag from a string name
   */
  Tag(const std::string name) : name_(name), isEmpty_(false) {}

  /*
   * @brief obtain the Tag name
   */
  std::string name() const { return name_; }

  /*
   * @brief check if the Tag is empty
   */
  bool is_empty() const { return isEmpty_; }
};

/*
 * @class TagComparator
 *
 * A comparator for sorting Tags.
 */
struct TagComparator {
  bool operator()(const Tag& lhs, const Tag& rhs) const {
    return lhs.name() < rhs.name();
  }
};

/*
 * @typedef TagList
 *
 * A std::vector of Tags.
 */
using TagList = std::vector<Tag>;

/*
 * @brief construct a TagList as the union of two TagLists
 * @param left the first TagList
 * @param right the second TagList
 *
 * The new TagList is constructed by copying 'left,' reserving room for 'right,'
 * and then adding each element of 'right.'
 */
TagList union_taglists(const TagList& left, const TagList& right) {
  TagList retval(left);
  retval.reserve(right.size());
  for (auto r : right) {
    retval.push_back(r);
  }
  return retval;
}

/*
 * @brief appends a Tag to the end of a TagList
 * @param list the existing TagList
 * @param newTag the new Tag to be added to the list
 */
TagList append_tag_to_list(const TagList& list, const Tag& newTag) {
  TagList retval(list);
  retval.reserve(1);
  retval.push_back(newTag);
  return retval;
}

/*
 * @brief make a TagList of a single Tag if that Tag is not empty
 */
TagList make_nonempty_single_tag_list(const Tag& tag) {
  TagList retval;
  if (tag.is_empty()) {
    retval.push_back(tag);
  }
  return retval;
}
}

#endif /* TAG_H_ */
