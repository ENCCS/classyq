#pragma once

#include "CRTP.hpp"

namespace classyq {
template <typename T> struct Printable : CRTP<T, Printable> {
  static constexpr bool is_printable = true;

  void print(std::ostream &os) const { os << this->underlying().get(); }
};

#if FLUENT_HOSTED == 1
template <typename T, typename Parameter, template <typename> class... Skills>
typename std::enable_if<NamedType<T, Parameter, Skills...>::is_printable, std::ostream &>::type operator<<(
    std::ostream &os,
    NamedType<T, Parameter, Skills...> const &object) {
  object.print(os);
  return os;
}
#endif
} // namespace classyq