#pragma once

#include <cassert>
#include <cstddef>

#include "tedder/field.hpp"

namespace tedder::detail
{
  // ---------------------------------------------------------------- local frame:

  // The local frame of ADR 0003 §3 for B = I/h. It maps an offset
  // r = g.offset(x, x0), already translated and wrapped by the domain, to the
  // local coordinates xi = r / h. It holds only h: x0 and the domain reach it
  // through the offset.
  template <Real T, std::size_t D>
  struct IsotropicFrame
  {
    T h; // bandwidth, in the units of the offset

    // One division per component, never a multiplication by 1 / h. The caller
    // checks admits_bandwidth(g, h) where the frame is built; this assert only
    // catches h <= 0 and NaN.
    constexpr Point<T, D> to_local(const Point<T, D> &offset) const noexcept
    {
      assert(h > T{0});
      Point<T, D> xi{};
      for (std::size_t d = 0; d < D; ++d)
        xi[d] = offset[d] / h;
      return xi;
    }
  };

  // q^2 = xi . xi, summed in axis order. It takes xi, not a physical offset:
  // distance_squared(g, x, x0) / h^2 is the rejected form, which fails once
  // h^2 underflows (ADR 0003 §4).
  template <Real T, std::size_t D>
  constexpr T squared_norm(const Point<T, D> &xi) noexcept
  {
    T sum = T{0};
    for (std::size_t d = 0; d < D; ++d)
      sum += xi[d] * xi[d];
    return sum;
  }

  // The candidate test: q^2 < 1. The boundary q^2 == 1 is out, and a NaN q^2
  // fails, because every comparison with NaN is false.
  template <Real T>
  constexpr bool is_candidate(const T q2) noexcept
  {
    return q2 < T{1};
  }
} // namespace tedder::detail
