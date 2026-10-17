#pragma once

#include <cassert>
#include <concepts>

#include "tedder/field.hpp"

namespace tedder
{
// ---------------------------------------------------------------- kernels:

// A kernel turns the local distance q into a weight w(q). The concept checks
// only the call. A model also promises (ADR 0003 §4):
// - precondition 0 <= q <= 1, which excludes NaN; w(1) == 0;
// - on [0, 1] the weight is finite and >= 0; on [0, 1) it is > 0 in exact
//   arithmetic, but may round to 0 near q = 1;
// - no normalising constant: only relative weights matter;
// - equal q give equal weights, and a call writes nothing.
// The concept cannot tell q from q * q. Kernels take q.
template <class K, class T>
concept Kernel = Real<T> && requires(const K &k, const T q) {
  { k(q) } noexcept -> std::same_as<T>;
};

// Wendland's C2 function (d = 3, k = 1), scaled so that w(0) == 1:
//   w(q) = (1 - q)^4 (4q + 1)
// Kept factored: the expanded polynomial cancels near q = 1 and can return
// negative weights.
struct WendlandC2
{
    template <Real T> static constexpr T operator()(const T q) noexcept
    {
      assert(q >= T{0} && q <= T{1});
      const T t = T{1} - q;
      const T t2 = t * t;
      return t2 * t2 * (T{4} * q + T{1});
    }
};

// Wendland's C4 function (d = 3, k = 2), scaled so that w(0) == 1:
//   w(q) = (1 - q)^6 (35q^2 + 18q + 3) / 3
struct WendlandC4
{
    template <Real T> static constexpr T operator()(const T q) noexcept
    {
      assert(q >= T{0} && q <= T{1});
      const T t = T{1} - q;
      const T t3 = t * t * t;
      constexpr T c = T{35} / T{3};
      return t3 * t3 * ((c * q + T{6}) * q + T{1});
    }
};
} // namespace tedder
