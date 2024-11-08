#pragma once

#include <utility>
#include <type_traits>

namespace detail {
/** Compile-time loop-application of a function to a variadic argument pack.
 *
 * @tparam F function to apply.
 * @tparam Args variadic argument pack.
 *
 * @param[in] F function to apply.
 * @param[in] Args variadic argument pack.
 */
template<typename F, typename... Args>
constexpr auto
constexpr_for(F&& f, Args&&... args) -> void
{
  (f(std::forward<Args>(args)), ...);
}

/** Compile-time loop-application of a function.
 *
 * @tparam Start loop start index.
 * @tparam End loop end index.
 * @tparam Increment loop increment.
 * @tparam F function to apply.
 *
 * @param[in] F function to apply.
 */
template<auto Start, auto End, auto Increment, typename F>
constexpr auto
constexpr_for(F&& f) -> void
{
  if constexpr (Start < End) {
    f(std::integral_constant<decltype(Start), Start>());
    constexpr_for<Start + Increment, End, Increment>(f);
  }
}
} // namespace detail