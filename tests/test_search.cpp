// Tests for the local frame, q^2 and the candidate test: T1 and T2 of the
// Stage 02 plan (docs/stages/02-neighbour-search/search-plan.md). The
// contract they check is ADR 0003 §3 and §4.
//
// What is tested, for one sample:
//   r  = g.offset(x, x0)      physical offset (computed by the domain; not
//                             tested here, so the tests pass r directly)
//   xi = frame.to_local(r)    local coordinates: xi = r / h, per axis
//   q2 = squared_norm(xi)     q^2 = xi_0^2 + xi_1^2 + ..., in axis order
//   is_candidate(q2)          true exactly when q^2 < 1
//
// Every expected value is written by hand, never computed with the code under
// test. Each test runs twice, for float and for double (TestType).

#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>

#include <tedder/detail/frame.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace tedder;
using tedder::detail::is_candidate;
using tedder::detail::IsotropicFrame;
using tedder::detail::squared_norm;

// These run at compile time. If the functions stopped being constexpr, or
// returned a wrong value, this file would not compile.
static_assert(IsotropicFrame<double, 2>{0.5}.to_local({1.0, -0.25}) == Point<double, 2>{2.0, -0.5});
static_assert(squared_norm(Point<float, 2>{3.0f, -4.0f}) == 25.0f);
static_assert(is_candidate(0.0) && !is_candidate(1.0));

// ---------------------------------------------------------------- T1: frame and q^2

namespace
{
  // One to_local example: bandwidth h, offset r, and the xi we expect.
  template <Real T>
  struct FrameCase
  {
    T h;
    Point<T, 3> r;
    Point<T, 3> xi;
  };

  // One squared_norm example: xi, and the q^2 we expect.
  template <Real T>
  struct NormCase
  {
    Point<T, 3> xi;
    T q2;
  };

  // 5/3 rounded once to T, computed in exact rational arithmetic.
  // A hex literal spells out the binary digits, so the value is exact, with
  // no decimal-to-binary conversion:
  // 5/3 = 1.101010..._2 = 0x1.aaaa...; the last hex digit is rounded.
  // In decimal: 1.6666666666666667 (double), 1.6666666f (float).
  template <Real T>
  constexpr T five_thirds = 0x1.aaaaaaaaaaaabp+0;
  template <>
  constexpr float five_thirds<float> = 0x1.aaaaaap+0f;

  // The axis-order case: xi = (5s, 7s, b), with b large and s small. A
  // "step" is the gap between b^2 and the next representable number.
  template <Real T>
  struct OrderCase
  {
    Point<T, 3> xi;
    T q2;    // b^2 + 37 steps: the sum in axis order
    T wrong; // b^2 + 36 steps: the sum when b^2 is added earlier
  };

  // double: b = 0.75, s = 2^-27, b^2 = 0.5625, step = 2^-53.
  template <Real T>
  constexpr OrderCase<T> order_case{{0x1.4p-25, 0x1.cp-25, 0.75}, 0x1.2000000000025p-1, 0x1.2000000000024p-1};
  // float: b = 0.5, s = 2^-13, b^2 = 0.25, step = 2^-25.
  template <>
  constexpr OrderCase<float> order_case<float>{{0x1.4p-11f, 0x1.cp-11f, 0.5f}, 0x1.00004ap-2f, 0x1.000048p-2f};
} // namespace

// Dividing by a power of two only changes the exponent, so r / h is exact and
// the expected xi can be written down by hand. For example, with h = 4:
// (2, -1, -12) / 4 = (0.5, -0.25, -3). Includes negative components, a zero,
// and a 1-D frame.
TEMPLATE_TEST_CASE("Frame: to_local divides by a power-of-two h exactly", "[search][frame]", float, double)
{
  using T = TestType;
  const FrameCase<T> cases[] = {
      {T{0.125}, {T{0.0625}, T{-0.5}, T{0.375}}, {T{0.5}, T{-4}, T{3}}},
      {T{1}, {T{-0.75}, T{0.25}, T{0}}, {T{-0.75}, T{0.25}, T{0}}},
      {T{4}, {T{2}, T{-1}, T{-12}}, {T{0.5}, T{-0.25}, T{-3}}},
  };
  for (const FrameCase<T> &c : cases)
  {
    CAPTURE(c.h);
    CHECK(IsotropicFrame<T, 3>{c.h}.to_local(c.r) == c.xi);
  }
  CHECK(IsotropicFrame<T, 1>{T{2}}.to_local({T{-3}}) == Point<T, 1>{T{-1.5}});
}

// ADR 0003 §4 fixes xi = r / h: one division, rounded once. A tempting
// shortcut computes 1 / h once and multiplies by it. That rounds twice (once
// for 1 / h, once for the product) and can change the last bit. For r = 5 and
// h = 3 it does, in both types, so this case catches the shortcut.
TEMPLATE_TEST_CASE("Frame: to_local divides, it does not multiply by 1 / h", "[search][frame]", float, double)
{
  using T = TestType;
  CHECK(IsotropicFrame<T, 1>{T{3}}.to_local({T{5}})[0] == five_thirds<T>);
  // Guard: the shortcut really gives a different value, so this case can
  // tell the two forms apart.
  CHECK(T{5} * (T{1} / T{3}) != five_thirds<T>);
}

// Squares and sums chosen so that every intermediate value is exactly
// representable. Then the result has no rounding error and must match
// exactly. For example, (3, -4, 12): 9 + 16 + 144 = 169.
TEMPLATE_TEST_CASE("Frame: squared_norm of hand-chosen xi is exact", "[search][frame]", float, double)
{
  using T = TestType;
  const NormCase<T> cases[] = {
      {{T{0}, T{0}, T{0}}, T{0}},
      {{T{0.5}, T{-0.25}, T{0}}, T{0.3125}},
      {{T{-1}, T{-1}, T{-1}}, T{3}},
      {{T{3}, T{-4}, T{12}}, T{169}},
  };
  for (const NormCase<T> &c : cases)
    CHECK(squared_norm(c.xi) == c.q2);
  CHECK(squared_norm(Point<T, 1>{T{-0.75}}) == T{0.5625});
}

// Rounding makes a floating-point sum depend on the order of its terms. The
// plan fixes axis order 0, 1, 2, ..., so that every caller gets the same q^2.
//
// xi = (5s, 7s, b). Measured in steps of b^2, the two small squares are
// 12.5 and 24.5 steps: each is half a step off a whole number.
// - In axis order, the two small squares are added first, exactly:
//   12.5 + 24.5 = 37 steps. Adding that to b^2 is exact too.
// - In either other order, b^2 comes earlier, and each small square is added
//   to it on its own. Each addition lands exactly halfway between two
//   representable numbers and rounds to the even one, losing half a step.
//   Twice: the sum ends at b^2 + 36 steps.
// Swapping axes 0 and 1 cannot change the sum, because a + b and b + a round
// the same. So this case tells apart every order that could change q^2.
TEMPLATE_TEST_CASE("Frame: squared_norm sums in axis order", "[search][frame]", float, double)
{
  using T = TestType;
  const OrderCase<T> &c = order_case<T>;
  CHECK(squared_norm(c.xi) == c.q2);
  // Guard: both orders that add b^2 earlier give b^2 + 36 steps.
  const T s0 = c.xi[0] * c.xi[0];
  const T s1 = c.xi[1] * c.xi[1];
  const T b2 = c.xi[2] * c.xi[2];
  CHECK((s0 + b2) + s1 == c.wrong);
  CHECK((s1 + b2) + s0 == c.wrong);
}

// A second order check, on a longer vector. xi = (1, a, a, a, a), with a^2
// tiny:
// - In axis order, 1 + a^2 rounds back to 1, four times, so q^2 = 1.
//   a^2 is at most half the gap between 1 and the next representable number,
//   and an exact tie rounds to the even neighbour, which is 1.
// - In reverse order, the four a^2 are added first, and 4 a^2 is big enough
//   to survive, so the sum exceeds 1. Other orders that add several a^2
//   before the 1 fail the same way.
// float: a = 2^-12, a^2 = 2^-24. double: a = 2^-27, a^2 = 2^-54.
// (0x1p-12 is the hex literal for 1 * 2^-12.)
TEMPLATE_TEST_CASE("Frame: squared_norm does not sum in reverse order", "[search][frame]", float, double)
{
  using T = TestType;
  const T a = std::same_as<T, float> ? T{0x1p-12f} : T{0x1p-27};
  const Point<T, 5> xi{T{1}, a, a, a, a};
  CHECK(squared_norm(xi) == T{1});
  // Guard: summing in reverse order really gives a different value, so this
  // case can tell the orders apart.
  T reversed = T{0};
  for (std::size_t d = xi.size(); d-- > 0;)
    reversed += xi[d] * xi[d];
  CHECK(reversed != T{1});
}

// ---------------------------------------------------------------- T2: candidate test

namespace
{
  // Bandwidths for the boundary tests. They must be powers of two: then
  // nextafter(h, 0) / h is exactly the largest value below 1. For other h the
  // quotient can land two steps below 1 (plan, T2).
  template <Real T>
  constexpr T power_of_two_h[] = {T{0.125}, T{1}, T{4}, T{1024}};
} // namespace

// The support is the open ball q < 1: a sample at distance exactly h is not a
// candidate (ADR 0003 §4). With a power-of-two h, xi = +1 or -1 exactly, so
// q^2 = 1 exactly. Checked along each axis, in both directions, for each h.
TEMPLATE_TEST_CASE("Candidate test: an offset of exactly h along an axis is out", "[search][candidate]", float, double)
{
  using T = TestType;
  for (const T h : power_of_two_h<T>)
    for (std::size_t axis = 0; axis < 3; ++axis)
      for (const T sign : {T{1}, T{-1}})
      {
        CAPTURE(h, axis, sign);
        Point<T, 3> r{};
        r[axis] = sign * h;
        const T q2 = squared_norm(IsotropicFrame<T, 3>{h}.to_local(r));
        CHECK(q2 == T{1});
        CHECK_FALSE(is_candidate(q2));
      }
}

// The other side of the boundary: the largest offset below h,
// nextafter(h, 0). Then xi is the largest value below 1, q^2 < 1, and the
// sample is a candidate. Together with the previous test, this puts the
// boundary exactly between these two neighbouring offsets.
TEMPLATE_TEST_CASE("Candidate test: the largest offset below h is in", "[search][candidate]", float, double)
{
  using T = TestType;
  const T below_one = std::nextafter(T{1}, T{0});
  for (const T h : power_of_two_h<T>)
    for (std::size_t axis = 0; axis < 3; ++axis)
      for (const T sign : {T{1}, T{-1}})
      {
        CAPTURE(h, axis, sign);
        Point<T, 3> r{};
        r[axis] = sign * std::nextafter(h, T{0});
        const Point<T, 3> xi = IsotropicFrame<T, 3>{h}.to_local(r);
        CHECK(xi[axis] == sign * below_one);
        const T q2 = squared_norm(xi);
        CHECK(q2 < T{1});
        CHECK(is_candidate(q2));
      }
}

// The sample at x0 itself (r = 0, so q^2 = 0) is a candidate. NaN and
// infinity are not: q^2 < 1 is false for both, because every comparison with
// NaN is false.
TEMPLATE_TEST_CASE("Candidate test: the centre, NaN and infinity", "[search][candidate]", float, double)
{
  using T = TestType;
  const T q2 = squared_norm(IsotropicFrame<T, 3>{T{1}}.to_local({T{0}, T{0}, T{0}}));
  CHECK(q2 == T{0});
  CHECK(is_candidate(q2));
  // NaN is passed straight to is_candidate. Producing it from coordinates
  // would need non-finite coordinates, which callers must not pass (A5).
  CHECK_FALSE(is_candidate(std::numeric_limits<T>::quiet_NaN()));
  CHECK_FALSE(is_candidate(std::numeric_limits<T>::infinity()));
}

// Extreme bandwidths: the two range examples of ADR 0003 §4. In both, an
// intermediate value overflows or underflows, and the decision must still be
// right. numeric_limits<T>::min() is the smallest positive normal number, and
// max() the largest finite one.
TEMPLATE_TEST_CASE("Candidate test: q^2 overflows far outside and xi underflows far inside", "[search][candidate]", float, double)
{
  using T = TestType;
  const T tiny = std::numeric_limits<T>::min();
  const T huge = std::numeric_limits<T>::max();

  // Far outside: h = min(), r = 1. xi_0 = 1 / min() is huge but finite; its
  // square overflows to infinity. Infinity is not < 1, so the sample is
  // rejected, which is the right answer.
  const T q2_outside = squared_norm(IsotropicFrame<T, 3>{tiny}.to_local({T{1}, T{0}, T{0}}));
  CHECK(q2_outside == std::numeric_limits<T>::infinity());
  CHECK_FALSE(is_candidate(q2_outside));

  // Far inside: h = max(), r = min(). xi_0 = min() / max() is too small to
  // represent and underflows to 0. Then q^2 = 0, so the sample is kept, which
  // is the right answer.
  const Point<T, 3> xi_inside = IsotropicFrame<T, 3>{huge}.to_local({tiny, T{0}, T{0}});
  CHECK(xi_inside[0] == T{0});
  CHECK(is_candidate(squared_norm(xi_inside)));
}
