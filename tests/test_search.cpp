// Tests for the local frame, q^2 and the candidate test: T1 and T2 of the
// Stage 02 plan (docs/stages/02-neighbour-search/search-plan.md). The
// contract they check is ADR 0003 §3 and §4.

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

// The frame, q^2 and the test are usable in constant expressions.
static_assert(IsotropicFrame<double, 2>{0.5}.to_local({1.0, -0.25}) == Point<double, 2>{2.0, -0.5});
static_assert(squared_norm(Point<float, 2>{3.0f, -4.0f}) == 25.0f);
static_assert(is_candidate(0.0) && !is_candidate(1.0));

// ---------------------------------------------------------------- T1: frame and q^2

namespace
{
  // Every expected value is written by hand. With h a power of two, r / h only
  // changes the exponent, so the quotients are exact.
  template <Real T>
  struct FrameCase
  {
    T h;
    Point<T, 3> r;
    Point<T, 3> xi;
  };

  template <Real T>
  struct NormCase
  {
    Point<T, 3> xi;
    T q2;
  };

  // 5 / 3 rounded once to T, from exact rational arithmetic.
  template <Real T>
  constexpr T five_thirds = 0x1.aaaaaaaaaaaabp+0;
  template <>
  constexpr float five_thirds<float> = 0x1.aaaaaap+0f;
} // namespace

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

// ADR 0003 §4 fixes one division per component. Multiplying by a rounded 1 / h
// rounds twice and can differ in the last bit, as it does for 5 / 3.
TEMPLATE_TEST_CASE("Frame: to_local divides, it does not multiply by 1 / h", "[search][frame]", float, double)
{
  using T = TestType;
  CHECK(IsotropicFrame<T, 1>{T{3}}.to_local({T{5}})[0] == five_thirds<T>);
  // The case tells the two forms apart.
  CHECK(T{5} * (T{1} / T{3}) != five_thirds<T>);
}

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

// Rounding makes the sum depend on its order. With xi = (1, a, a, a, a), each
// a^2 is at most half an ulp of 1, so adding it to 1 leaves 1 (a tie rounds to
// even). Summed first, the four a^2 exceed half an ulp.
// float: a^2 = 2^-24; double: a^2 = 2^-54.
TEMPLATE_TEST_CASE("Frame: squared_norm sums in axis order", "[search][frame]", float, double)
{
  using T = TestType;
  const T a = std::same_as<T, float> ? T{0x1p-12f} : T{0x1p-27};
  const Point<T, 5> xi{T{1}, a, a, a, a};
  CHECK(squared_norm(xi) == T{1});
  // The case tells the orders apart.
  T reversed = T{0};
  for (std::size_t d = xi.size(); d-- > 0;)
    reversed += xi[d] * xi[d];
  CHECK(reversed != T{1});
}

// ---------------------------------------------------------------- T2: candidate test

namespace
{
  // Powers of two, so that nextafter(h, 0) / h is exactly the largest value
  // below 1 (plan, T2).
  template <Real T>
  constexpr T power_of_two_h[] = {T{0.125}, T{1}, T{4}, T{1024}};
} // namespace

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

TEMPLATE_TEST_CASE("Candidate test: the centre, NaN and infinity", "[search][candidate]", float, double)
{
  using T = TestType;
  const T q2 = squared_norm(IsotropicFrame<T, 3>{T{1}}.to_local({T{0}, T{0}, T{0}}));
  CHECK(q2 == T{0});
  CHECK(is_candidate(q2));
  // NaN is passed directly: non-finite coordinates would break a precondition.
  CHECK_FALSE(is_candidate(std::numeric_limits<T>::quiet_NaN()));
  CHECK_FALSE(is_candidate(std::numeric_limits<T>::infinity()));
}

// The range examples of ADR 0003 §4.
TEMPLATE_TEST_CASE("Candidate test: q^2 overflows far outside and xi underflows far inside", "[search][candidate]", float, double)
{
  using T = TestType;
  const T tiny = std::numeric_limits<T>::min();
  const T huge = std::numeric_limits<T>::max();

  // Far outside: h = min(), r = 1. xi_0 = 1 / min() is finite, its square is not.
  const T q2_outside = squared_norm(IsotropicFrame<T, 3>{tiny}.to_local({T{1}, T{0}, T{0}}));
  CHECK(q2_outside == std::numeric_limits<T>::infinity());
  CHECK_FALSE(is_candidate(q2_outside));

  // Far inside: h = max(), r = min(). xi_0 underflows to 0.
  const Point<T, 3> xi_inside = IsotropicFrame<T, 3>{huge}.to_local({tiny, T{0}, T{0}});
  CHECK(xi_inside[0] == T{0});
  CHECK(is_candidate(squared_norm(xi_inside)));
}
