#pragma once

namespace classyq {
template <typename T, template <typename> class CRTPType> struct CRTP {
  constexpr T &underlying() { return static_cast<T &>(*this); }
  constexpr T const &underlying() const { return static_cast<T const &>(*this); }
};
} // namespace classyq