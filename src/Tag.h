#ifndef TAG_H_
#define TAG_H_

#include <string>
#include <vector>

namespace ast {

class Tag {
 protected:
  std::string name_;

 public:
  Tag() : name_("") {}
  Tag(const std::string name) : name_(name) {}
  std::string name() const { return name_; }
};

struct TagComparator {
  bool operator()(const Tag& lhs, const Tag& rhs) const {
    return lhs.name() < rhs.name();
  }
};

using TagList = std::vector<Tag>;

TagList union_taglists(TagList left, TagList right) {
  TagList retval;
  retval.reserve(left.size() + right.size());
  for (auto l : left) {
    retval.push_back(l);
  }
  for (auto r : right) {
    retval.push_back(r);
  }
  return retval;
}
}

#endif /* TAG_H_ */
