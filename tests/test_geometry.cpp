#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include <tedder/field.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace tedder;

// ---------------------------------------------------------------- generic laws
//
// These hold for any Domain, evaluated at coordinates chosen to stay well
// inside a single minimum image (a quarter of max_bandwidth() when finite,
// a fixed moderate scale otherwise), so a Periodic domain behaves like an
// Euclidean one here and no wrap is exercised. Wrap-crossing behaviour has
// its own dedicated checks further down.

namespace
{
template <Domain G> void check_generic_domain_laws(const G &g)
{
  using T = typename G::scalar_type;
  using P = typename G::point_type;
  constexpr std::size_t D = G::dimension;

  const T mb = g.max_bandwidth();
  const T scale = std::isfinite(mb) ? mb / T{4} : T{10};
  const T eps = std::numeric_limits<T>::epsilon();
  const T tol = T{8} * eps * (scale > T{1} ? scale : T{1});

  P p{}, q{};
  for (std::size_t d = 0; d < D; ++d)
    {
      p[d] = static_cast<T>(d) * scale * T{0.25};
      q[d] = p[d] + (d % 2 == 0 ? scale * T{0.5} : -scale * T{0.5});
    }

  // offset(p,p) is zero in every component; distance_squared(p,p) is zero.
  const auto zero = g.offset(p, p);
  for (std::size_t d = 0; d < D; ++d)
    CHECK(zero[d] == T{0});
  CHECK(g.distance_squared(p, p) == T{0});

  // distance_squared is never negative and is symmetric in its arguments.
  const T d2pq = g.distance_squared(p, q);
  const T d2qp = g.distance_squared(q, p);
  CHECK(d2pq >= T{0});
  CHECK(d2pq == d2qp);

  // offset(q,p) is the negation of offset(p,q); p,q sit far from any tie.
  const auto opq = g.offset(p, q);
  const auto oqp = g.offset(q, p);
  for (std::size_t d = 0; d < D; ++d)
    CHECK(oqp[d] == -opq[d]);

  // distance_squared matches an independently summed squared norm of offset.
  T norm2{};
  for (std::size_t d = 0; d < D; ++d)
    norm2 += opq[d] * opq[d];
  CHECK(d2pq == norm2);

  // Translation invariance for a translation well inside the safe scale,
  // with positive, negative and mixed-sign components.
  P t{};
  for (std::size_t d = 0; d < D; ++d)
    t[d] = (d % 2 == 0 ? T{1} : T{-1}) * scale * T{0.125} * static_cast<T>(d + 1);

  P pt{}, qt{};
  for (std::size_t d = 0; d < D; ++d)
    {
      pt[d] = p[d] + t[d];
      qt[d] = q[d] + t[d];
    }

  const auto shifted = g.offset(pt, qt);
  for (std::size_t d = 0; d < D; ++d)
    CHECK(std::abs(shifted[d] - opq[d]) <= tol);
}

// A known non-zero distance, computed by hand, independent of the domain's
// own reduction. For Periodic domains the caller must pick a period large
// enough that this pair never wraps.
template <Domain G>
void check_known_distance(const G &g, const typename G::point_type &p, const typename G::point_type &q, typename G::scalar_type expected_distance)
{
  using T = typename G::scalar_type;
  const T expected_sq = expected_distance * expected_distance;
  CHECK(g.distance_squared(p, q) == expected_sq);
  CHECK(distance(g, p, q) == expected_distance);
}

template <class T, std::size_t D> struct KnownPair
{
  Point<T, D> p{};
  Point<T, D> q{};
  T distance{};
};

// 1D: 0 to 5. 2D: the 3-4-5 triangle. 3D: 2-3-6-7, another integer
// Pythagorean quadruple. All coordinates and the result are small exact
// integers, representable identically in float and double.
template <Real T, std::size_t D> KnownPair<T, D> make_known_pair()
{
  KnownPair<T, D> kp;
  if constexpr (D == 1)
    {
      kp.p = Point<T, 1>{T{0}};
      kp.q = Point<T, 1>{T{5}};
      kp.distance = T{5};
    }
  else if constexpr (D == 2)
    {
      kp.p = Point<T, 2>{T{0}, T{0}};
      kp.q = Point<T, 2>{T{3}, T{4}};
      kp.distance = T{5};
    }
  else
    {
      static_assert(D == 3);
      kp.p = Point<T, 3>{T{0}, T{0}, T{0}};
      kp.q = Point<T, 3>{T{2}, T{3}, T{6}};
      kp.distance = T{7};
    }
  return kp;
}

template <Real T, std::size_t D> void run_euclidean_checks()
{
  Euclidean<T, D> g;
  check_generic_domain_laws(g);
  const auto kp = make_known_pair<T, D>();
  check_known_distance(g, kp.p, kp.q, kp.distance);
}

template <Real T, std::size_t D> void run_periodic_checks_no_wrap()
{
  Periodic<T, D> g;
  for (std::size_t d = 0; d < D; ++d)
    g.length[d] = T{100} + static_cast<T>(d); // generous: nothing here wraps
  REQUIRE(g.is_valid());

  check_generic_domain_laws(g);
  const auto kp = make_known_pair<T, D>();
  check_known_distance(g, kp.p, kp.q, kp.distance);
}
} // namespace

TEST_CASE("Domain: generic laws hold for Euclidean<double, D>", "[geometry][euclidean]")
{
  run_euclidean_checks<double, 1>();
  run_euclidean_checks<double, 2>();
  run_euclidean_checks<double, 3>();
}

TEST_CASE("Domain: generic laws hold for Euclidean<float, D>", "[geometry][euclidean]")
{
  run_euclidean_checks<float, 1>();
  run_euclidean_checks<float, 2>();
  run_euclidean_checks<float, 3>();
}

TEST_CASE("Domain: generic laws hold for Periodic<double, D> away from any wrap", "[geometry][periodic]")
{
  run_periodic_checks_no_wrap<double, 1>();
  run_periodic_checks_no_wrap<double, 2>();
  run_periodic_checks_no_wrap<double, 3>();
}

TEST_CASE("Domain: generic laws hold for Periodic<float, D> away from any wrap", "[geometry][periodic]")
{
  run_periodic_checks_no_wrap<float, 1>();
  run_periodic_checks_no_wrap<float, 2>();
  run_periodic_checks_no_wrap<float, 3>();
}

// ---------------------------------------------------------------- Periodic: wrapping

namespace
{
// Shifting a point by a whole number of periods on every axis must leave
// both the offset and the distance unchanged, away from half-period ties.
template <Real T, std::size_t D> void check_periodic_period_shift_invariance()
{
  Periodic<T, D> g;
  Point<T, D> p{}, q{};
  for (std::size_t d = 0; d < D; ++d)
    {
      g.length[d] = T{6} + static_cast<T>(d);
      p[d] = T{1} + static_cast<T>(d) * T{0.5}; // comfortably below half a period
      q[d] = T{0};
    }
  REQUIRE(g.is_valid());

  const auto o0 = g.offset(p, q);
  const T d20 = g.distance_squared(p, q);

  Point<T, D> p_shifted{};
  for (std::size_t d = 0; d < D; ++d)
    p_shifted[d] = p[d] + static_cast<T>(d + 1) * g.length[d]; // shift by (d+1) whole periods

  const auto o1 = g.offset(p_shifted, q);
  const T d21 = g.distance_squared(p_shifted, q);

  const T eps = std::numeric_limits<T>::epsilon();
  for (std::size_t d = 0; d < D; ++d)
    CHECK(std::abs(o1[d] - o0[d]) <= T{8} * eps * g.length[d]);
  CHECK(std::abs(d21 - d20) <= T{8} * eps * std::max(d20, T{1}));
}

// An open axis (length[d] <= 0) is exact, unwrapped subtraction, no matter
// how far apart the coordinates are. With D > 1 this also exercises a mixed
// box: axis 0 stays open while the remaining axes wrap.
template <Real T, std::size_t D> void check_periodic_open_axis_does_not_wrap()
{
  Periodic<T, D> g;
  for (std::size_t d = 0; d < D; ++d)
    g.length[d] = (d == 0) ? T{-1} : T{4};
  REQUIRE(g.is_valid());

  Point<T, D> p{}, q{};
  p[0] = T{123.5};

  const auto r = g.offset(p, q);
  CHECK(r[0] == p[0] - q[0]);
}

template <Real T> void check_periodic_wrapped_bound_1d()
{
  Periodic<T, 1> g;
  g.length[0] = T{5};
  REQUIRE(g.is_valid());

  const T eps = std::numeric_limits<T>::epsilon();
  const T half = g.length[0] / T{2};
  const T tol = T{8} * eps * g.length[0];

  const T multipliers[] = {T{-3.25}, T{-1.0}, T{-0.5}, T{0.0}, T{0.5}, T{0.9375}, T{1.0}, T{2.75}};
  for (T m : multipliers)
    {
      const Point<T, 1> p{m * g.length[0]};
      const Point<T, 1> q{};
      const auto r = g.offset(p, q);
      CHECK(std::abs(r[0]) <= half + tol);
    }
}

// Faces and corners of a mixed 2D box: axis 0 wraps, axis 1 is open.
template <Real T> void check_periodic_wrapped_bound_2d()
{
  Periodic<T, 2> g;
  g.length = {T{4}, T{-1}};
  REQUIRE(g.is_valid());

  const T eps = std::numeric_limits<T>::epsilon();
  const T half0 = g.length[0] / T{2};
  const T tol0 = T{8} * eps * g.length[0];

  const T multipliers[] = {T{-2.75}, T{-1.0}, T{-0.5}, T{0.0}, T{0.5}, T{0.9375}, T{1.0}, T{3.25}};
  const T open_axis_values[] = {T{0}, T{57.25}};

  for (T m : multipliers)
    for (T open_val : open_axis_values)
      {
        const Point<T, 2> p{m * g.length[0], open_val};
        const Point<T, 2> q{};
        const auto r = g.offset(p, q);

        CHECK(std::abs(r[0]) <= half0 + tol0); // wrapping axis: minimum image
        CHECK(r[1] == p[1] - q[1]);            // open axis: exact, unwrapped
      }
}

// Faces and corners of a mixed 3D box: axes 0,1 wrap, axis 2 is open.
template <Real T> void check_periodic_wrapped_bound_3d()
{
  Periodic<T, 3> g;
  g.length = {T{4}, T{6}, T{-1}};
  REQUIRE(g.is_valid());

  const T eps = std::numeric_limits<T>::epsilon();
  const T half0 = g.length[0] / T{2};
  const T half1 = g.length[1] / T{2};
  const T tol0 = T{8} * eps * g.length[0];
  const T tol1 = T{8} * eps * g.length[1];

  const T multipliers[] = {T{-2.75}, T{-0.5}, T{0.0}, T{0.5}, T{3.25}};
  const T open_axis_values[] = {T{0}, T{-41.5}};

  for (T m0 : multipliers)
    for (T m1 : multipliers)
      for (T open_val : open_axis_values)
        {
          const Point<T, 3> p{m0 * g.length[0], m1 * g.length[1], open_val};
          const Point<T, 3> q{};
          const auto r = g.offset(p, q);

          CHECK(std::abs(r[0]) <= half0 + tol0);
          CHECK(std::abs(r[1]) <= half1 + tol1);
          CHECK(r[2] == p[2] - q[2]);
        }
}

// admits_bandwidth is strict at max_bandwidth(), for a domain where that
// limit is a normal, finite, representable value.
template <Real T, std::size_t D> void check_admits_bandwidth_boundary()
{
  Periodic<T, D> g;
  for (std::size_t d = 0; d < D; ++d)
    g.length[d] = T{2} + static_cast<T>(d); // axis 0 is always the smallest: limit == 1
  REQUIRE(g.is_valid());

  const T limit = g.max_bandwidth();
  REQUIRE(std::isfinite(limit));
  REQUIRE(limit > T{0});

  CHECK(admits_bandwidth(g, std::nextafter(limit, T{0})));
  CHECK_FALSE(admits_bandwidth(g, limit));
  CHECK_FALSE(admits_bandwidth(g, std::nextafter(limit, std::numeric_limits<T>::infinity())));

  CHECK_FALSE(admits_bandwidth(g, T{0}));
  CHECK_FALSE(admits_bandwidth(g, T{-1}));
  CHECK_FALSE(admits_bandwidth(g, std::numeric_limits<T>::quiet_NaN()));
  CHECK_FALSE(admits_bandwidth(g, std::numeric_limits<T>::infinity()));
  CHECK_FALSE(admits_bandwidth(g, -std::numeric_limits<T>::infinity()));
}
} // namespace

TEST_CASE("Periodic: shifting by whole periods leaves distance and offset unchanged", "[geometry][periodic]")
{
  check_periodic_period_shift_invariance<double, 1>();
  check_periodic_period_shift_invariance<double, 2>();
  check_periodic_period_shift_invariance<double, 3>();
  check_periodic_period_shift_invariance<float, 1>();
  check_periodic_period_shift_invariance<float, 2>();
  check_periodic_period_shift_invariance<float, 3>();
}

TEST_CASE("Periodic: an open axis never wraps", "[geometry][periodic]")
{
  check_periodic_open_axis_does_not_wrap<double, 1>();
  check_periodic_open_axis_does_not_wrap<double, 2>();
  check_periodic_open_axis_does_not_wrap<double, 3>();
  check_periodic_open_axis_does_not_wrap<float, 1>();
  check_periodic_open_axis_does_not_wrap<float, 2>();
  check_periodic_open_axis_does_not_wrap<float, 3>();
}

TEST_CASE("Periodic: every wrapped offset component stays within half a period", "[geometry][periodic]")
{
  check_periodic_wrapped_bound_1d<double>();
  check_periodic_wrapped_bound_2d<double>();
  check_periodic_wrapped_bound_3d<double>();
  check_periodic_wrapped_bound_1d<float>();
  check_periodic_wrapped_bound_2d<float>();
  check_periodic_wrapped_bound_3d<float>();
}

TEST_CASE("Periodic: admits_bandwidth is strict at the boundary", "[geometry][periodic][bandwidth]")
{
  check_admits_bandwidth_boundary<double, 1>();
  check_admits_bandwidth_boundary<double, 2>();
  check_admits_bandwidth_boundary<double, 3>();
  check_admits_bandwidth_boundary<float, 1>();
  check_admits_bandwidth_boundary<float, 2>();
  check_admits_bandwidth_boundary<float, 3>();
}

// max_bandwidth() must be the smallest positive period divided by two, not
// merely the first one encountered while scanning axes. Every test above
// this point uses periods in increasing order, so it cannot distinguish
// "smallest positive" from "first positive" or "first axis". Each case below
// places the shortest positive period in a different position, including one
// with an open axis ahead of it.
namespace
{
template <Real T> void check_max_bandwidth_smallest_period_2d()
{
  {
    Periodic<T, 2> g;
    g.length = {T{2}, T{6}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
  {
    Periodic<T, 2> g;
    g.length = {T{6}, T{2}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
  {
    // Open axis precedes the shortest period: distinguishes "first positive"
    // from "first axis" as well.
    Periodic<T, 2> g;
    g.length = {T{-1}, T{2}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
}

template <Real T> void check_max_bandwidth_smallest_period_3d()
{
  {
    Periodic<T, 3> g;
    g.length = {T{2}, T{6}, T{4}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
  {
    Periodic<T, 3> g;
    g.length = {T{6}, T{2}, T{4}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
  {
    Periodic<T, 3> g;
    g.length = {T{6}, T{4}, T{2}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
  {
    // Open axis precedes the shortest period.
    Periodic<T, 3> g;
    g.length = {T{0}, T{6}, T{2}};
    REQUIRE(g.is_valid());
    CHECK(g.max_bandwidth() == T{1});
  }
}

// The consequence, not only the value: a bandwidth admissible under a
// "first positive period" mutation must actually be rejected.
template <Real T> void check_admits_bandwidth_follows_smallest_period()
{
  Periodic<T, 3> g;
  g.length = {T{6}, T{2}, T{4}};
  REQUIRE(g.is_valid());
  CHECK_FALSE(admits_bandwidth(g, T{2}));
}
} // namespace

TEST_CASE("Periodic: max_bandwidth is governed by the smallest positive period, not the first",
          "[geometry][periodic][bandwidth]")
{
  check_max_bandwidth_smallest_period_2d<double>();
  check_max_bandwidth_smallest_period_2d<float>();
  check_max_bandwidth_smallest_period_3d<double>();
  check_max_bandwidth_smallest_period_3d<float>();
  check_admits_bandwidth_follows_smallest_period<double>();
  check_admits_bandwidth_follows_smallest_period<float>();
}

// max_bandwidth() must reflect the current length on every call, not a value
// cached from an earlier state.
TEST_CASE("Periodic: max_bandwidth follows a changed period, not a cached value", "[geometry][periodic][bandwidth]")
{
  Periodic<double, 2> g;
  g.length = {4.0, 4.0};
  REQUIRE(g.is_valid());
  CHECK(g.max_bandwidth() == 2.0);

  g.length[0] = 2.0;
  CHECK(g.max_bandwidth() == 1.0);

  g.length[0] = 10.0; // axis 1 (still 4) now governs again
  CHECK(g.max_bandwidth() == 2.0);
}

// ---------------------------------------------------------------- Periodic: validity (B2)

namespace
{
template <Real T, std::size_t D> void check_is_valid_positive_cases()
{
  const Periodic<T, D> default_g;
  CHECK(default_g.is_valid());

  Periodic<T, D> all_open;
  for (std::size_t d = 0; d < D; ++d)
    all_open.length[d] = T{-1}; // finite negative: still open, still valid
  CHECK(all_open.is_valid());

  Periodic<T, D> all_wrap;
  for (std::size_t d = 0; d < D; ++d)
    all_wrap.length[d] = T{2} + static_cast<T>(d);
  CHECK(all_wrap.is_valid());

  Periodic<T, D> mixed;
  for (std::size_t d = 0; d < D; ++d)
    mixed.length[d] = (d % 2 == 0) ? T{3} + static_cast<T>(d) : T{0};
  CHECK(mixed.is_valid());
}

// One axis gets a non-finite value; every other axis keeps a valid, finite,
// wrapping length, so a single bad position invalidates the whole domain.
template <Real T, std::size_t D> void check_invalid_axis_position(std::size_t bad_axis, T bad_value)
{
  Periodic<T, D> g;
  for (std::size_t d = 0; d < D; ++d)
    g.length[d] = T{2} + static_cast<T>(d);
  g.length[bad_axis] = bad_value;

  CHECK_FALSE(g.is_valid());
  CHECK(g.max_bandwidth() == T{0});

  // Every tested bandwidth is rejected, including a small positive value
  // that would be admissible if this domain were valid.
  CHECK_FALSE(admits_bandwidth(g, T{0.5}));
  CHECK_FALSE(admits_bandwidth(g, T{1}));
  CHECK_FALSE(admits_bandwidth(g, T{0}));
  CHECK_FALSE(admits_bandwidth(g, T{-1}));
  CHECK_FALSE(admits_bandwidth(g, std::numeric_limits<T>::quiet_NaN()));
  CHECK_FALSE(admits_bandwidth(g, std::numeric_limits<T>::infinity()));
}

template <Real T, std::size_t D> void check_invalid_axis_all_positions()
{
  for (std::size_t d = 0; d < D; ++d)
    {
      check_invalid_axis_position<T, D>(d, std::numeric_limits<T>::quiet_NaN());
      check_invalid_axis_position<T, D>(d, std::numeric_limits<T>::infinity());
      check_invalid_axis_position<T, D>(d, -std::numeric_limits<T>::infinity());
    }
}

// Validity must be recomputed from current `length`, never cached.
template <Real T> void check_validity_mutation_not_cached()
{
  Periodic<T, 2> g;
  CHECK(g.is_valid()); // default: valid

  g.length[0] = T{4};
  CHECK(g.is_valid()); // still valid: finite positive

  g.length[1] = std::numeric_limits<T>::infinity();
  CHECK_FALSE(g.is_valid()); // now invalid

  g.length[1] = T{-2}; // finite negative: open, valid again
  CHECK(g.is_valid());

  g.length[0] = std::numeric_limits<T>::quiet_NaN();
  CHECK_FALSE(g.is_valid()); // invalid again

  g.length[0] = T{0};
  CHECK(g.is_valid()); // back to valid
}
} // namespace

TEST_CASE("Periodic: is_valid is true for default, open, wrapping and mixed configurations", "[geometry][periodic][validity]")
{
  check_is_valid_positive_cases<double, 1>();
  check_is_valid_positive_cases<double, 2>();
  check_is_valid_positive_cases<double, 3>();
  check_is_valid_positive_cases<float, 1>();
  check_is_valid_positive_cases<float, 2>();
  check_is_valid_positive_cases<float, 3>();
}

TEST_CASE("Periodic: a non-finite value on any single axis invalidates the whole domain", "[geometry][periodic][validity]")
{
  check_invalid_axis_all_positions<double, 1>();
  check_invalid_axis_all_positions<double, 2>();
  check_invalid_axis_all_positions<double, 3>();
  check_invalid_axis_all_positions<float, 1>();
  check_invalid_axis_all_positions<float, 2>();
  check_invalid_axis_all_positions<float, 3>();
}

TEST_CASE("Periodic: validity follows current length, not a cached value", "[geometry][periodic][validity]")
{
  check_validity_mutation_not_cached<double>();
  check_validity_mutation_not_cached<float>();
}

// ---------------------------------------------------------------- Periodic: reduction (B3)

namespace
{
// offset(p,q) on a single axis, checked against an independently derived
// expected value, plus its exact reversal (antisymmetry holds bit-exactly
// here: round() is odd and a-b/b-a negate exactly, so the same branch is
// always taken for both argument orders).
template <Real T> void check_offset(T L, T p, T q, T expected)
{
  Periodic<T, 1> g;
  g.length[0] = L;
  REQUIRE(g.is_valid());

  const Point<T, 1> pp{p}, qq{q};
  CHECK(g.offset(pp, qq)[0] == expected);
  CHECK(g.offset(qq, pp)[0] == -expected);
}
} // namespace

TEST_CASE("Periodic: F3 regressions from field-plan.md", "[geometry][periodic][reduction]")
{
  // Authoritative against the current header: if any of these four fails,
  // stop and report it rather than adjusting the expected value.
  check_offset<double>(3.0, 1e16, 0.0, 1.0);
  check_offset<float>(3.0f, 1e8f, 0.0f, 1.0f);
  check_offset<float>(std::nextafterf(1.0f, 2.0f), 4000.25f, 0.0f, 0.249523162841796875f);
  check_offset<double>(1e308, 1e308, -1e308, 0.0);
}

TEST_CASE("Periodic: reduction with negative coordinates far outside the box", "[geometry][periodic][reduction]")
{
  // -1000000 is an exact multiple of the period; both residues are small
  // exact integers derived independently of the header's own algorithm.
  check_offset<double>(4.0, -1000000.0, 3.0, 1.0);
  check_offset<float>(4.0f, -1000000.0f, 3.0f, 1.0f);

  check_offset<double>(6.0, -53.0, -1.0, 2.0);
  check_offset<float>(6.0f, -53.0f, -1.0f, 2.0f);
}

TEST_CASE("Periodic: reduction for several non-unit periods", "[geometry][periodic][reduction]")
{
  check_offset<double>(2.5, 10.0, 1.0, -1.0);
  check_offset<float>(2.5f, 10.0f, 1.0f, -1.0f);

  check_offset<double>(7.0, 22.0, 2.0, -1.0);
  check_offset<float>(7.0f, 22.0f, 2.0f, -1.0f);
}

TEST_CASE("Periodic: reduction stays bounded where raw subtraction would overflow", "[geometry][periodic][reduction]")
{
  // double is already covered by the authoritative 1e308 case above. For
  // float, build the overflowing pair from the period itself (L*2 and its
  // negation), so the relationship is exact regardless of literal rounding.
  const float L = 1e38f;
  const float p = L * 2.0f; // doubling is exact: no rounding introduced
  const float q = -p;
  check_offset<float>(L, p, q, 0.0f);
}

TEST_CASE("Periodic: representation ties follow the documented convention, not a global sign", "[geometry][periodic][reduction]")
{
  // At an exact half-period separation both images are equally short, so the
  // sign is a convention. Antisymmetry, representation-invariance and
  // translation-invariance cannot all hold at a tie. The current reduction
  // keeps antisymmetry and translation-invariance and gives up
  // representation-invariance: for L = 2, offset(0,1) is +1 and offset(2,1)
  // is -1.
  Periodic<double, 1> g;
  g.length[0] = 2.0;
  REQUIRE(g.is_valid());

  CHECK(g.offset({2.0}, {1.0})[0] == -1.0);
}

// The reduction has two boundaries: |t| = L, where the fast path hands over
// to the std::remainder fallback, and |t| = L/2, where the minimum image
// flips. L = 4 is a power of two, so every value here is exact, and the
// neighbouring representable values straddle each boundary exactly.
namespace
{
template <Real T> void check_reduction_boundary_neighbours()
{
  const T L = T{4};
  const T below = L - std::nextafter(L, T{0});                                   // ulp just below L
  const T above = std::nextafter(L, std::numeric_limits<T>::infinity()) - L;     // ulp just above L

  Periodic<T, 2> g;
  g.length = {L, T{-1}};
  REQUIRE(g.is_valid());

  const Point<T, 2> q{T{0}, T{0}};

  // Just below L: fast path, offset is -below.
  {
    const Point<T, 2> p{std::nextafter(L, T{0}), T{0}};
    CHECK(g.offset(p, q)[0] == -below);
  }
  // Exactly L: fast path, offset is 0.
  {
    const Point<T, 2> p{L, T{0}};
    CHECK(g.offset(p, q)[0] == T{0});
  }
  // Just above L: first value to take the slow (std::remainder) path.
  {
    const Point<T, 2> p{std::nextafter(L, std::numeric_limits<T>::infinity()), T{0}};
    CHECK(g.offset(p, q)[0] == above);
  }
  // Just below L/2: fast path, unchanged.
  {
    const T half_below = std::nextafter(L / 2, T{0});
    const Point<T, 2> p{half_below, T{0}};
    CHECK(g.offset(p, q)[0] == half_below);
  }
  // Just above L/2: fast path, minimum image flips.
  {
    const T half_above = std::nextafter(L / 2, L);
    const Point<T, 2> p{half_above, T{0}};
    CHECK(g.offset(p, q)[0] == half_above - L);
  }
}
} // namespace

TEST_CASE("Periodic: reduction is exact at the values adjacent to both boundaries", "[geometry][periodic][reduction]")
{
  check_reduction_boundary_neighbours<double>();
  check_reduction_boundary_neighbours<float>();
}
