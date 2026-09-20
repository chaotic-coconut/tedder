#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>

#include <tedder/field.hpp>

#include <catch2/catch_test_macros.hpp>

using Catch::Matchers::WithinAbs;

using namespace tedder;

// Every model must satisfy the concept, at every type and dimension we use.
static_assert(Domain<Euclidean<double, 1>>);
static_assert(Domain<Euclidean<double, 2>>);
static_assert(Domain<Euclidean<double, 3>>);
static_assert(Domain<Euclidean<float, 2>>);
static_assert(Domain<Periodic<double, 1>>);
static_assert(Domain<Periodic<double, 2>>);
static_assert(Domain<Periodic<double, 3>>);
static_assert(Domain<Periodic<float, 3>>);

// And things that are not domains must not, or the concept checks nothing.
namespace
{
struct Empty
{
};
} // namespace

static_assert(!Domain<Empty>);

TEST_CASE("Euclidean: known values", "[domain][euclidean]")
{
  Euclidean<double, 2> g;
  const Point<double, 2> p{0.0, 0.0}, q{3.0, 4.0};

  CHECK(g.distance_squared(p, q) == 25.0);
  CHECK(distance(g, p, q) == 5.0);
  CHECK(g.distance_squared(p, p) == 0.0);
}

TEST_CASE("Euclidean: offset is p - q", "[domain][euclidean]")
{
  Euclidean<double, 2> g;
  const auto r = g.offset({1.0, 2.0}, {4.0, 6.0});
  CHECK(r[0] == -3.0);
  CHECK(r[1] == -4.0);
}

TEST_CASE("Euclidean: no bandwidth limit", "[domain][euclidean]")
{
  CHECK(Euclidean<double, 2>{}.max_bandwidth() == std::numeric_limits<double>::infinity());
  CHECK(Euclidean<float, 3>{}.max_bandwidth() == std::numeric_limits<float>::infinity());
}

TEST_CASE("admits_bandwidth rejects degenerate values", "[domain][bandwidth]")
{
  Euclidean<double, 2> g;
  CHECK(admits_bandwidth(g, 1e9));
  CHECK_FALSE(admits_bandwidth(g, 0.0));
  CHECK_FALSE(admits_bandwidth(g, -1.0));
  CHECK_FALSE(admits_bandwidth(g, std::numeric_limits<double>::quiet_NaN()));
  CHECK_FALSE(admits_bandwidth(g, std::numeric_limits<double>::infinity()));
}

// ---------------------------------------------------------------- Periodic

TEST_CASE("Periodic: default construction does not wrap", "[domain][periodic]")
{
  // All-zero lengths must behave exactly as Euclidean, not produce NaN.
  Periodic<double, 3> z{};
  Euclidean<double, 3> e;
  const Point<double, 3> p{0.1, 0.2, 0.3}, q{0.9, 0.8, 0.7};

  CHECK(z.offset(p, q) == e.offset(p, q));
  CHECK(z.distance_squared(p, q) == e.distance_squared(p, q));
  CHECK(z.max_bandwidth() == std::numeric_limits<double>::infinity());
}

TEST_CASE("Periodic: offset is p - q when nothing needs wrapping", "[domain][periodic]")
{
  Periodic<double, 2> g{{10.0, 10.0}};
  const auto r = g.offset({1.0, 2.0}, {4.0, 6.0});

  CHECK(r[0] == -3.0);
  CHECK(r[1] == -4.0);
}

TEST_CASE("Periodic: wraps across the boundary", "[domain][periodic]")
{
  Periodic<double, 1> g{{1.0}};
  const Point<double, 1> a{0.1}, b{0.9};

  // Straight line is 0.8, through the boundary it is 0.2.
  CHECK_THAT(distance(g, a, b), WithinAbs(0.2, 1e-15));
  CHECK_THAT(g.offset(a, b)[0], WithinAbs(0.2, 1e-15));
}

TEST_CASE("Periodic: only the wrapping axes wrap", "[domain][periodic]")
{
  Periodic<double, 2> g{{1.0, 0.0}}; // x wraps, y does not
  const auto r = g.offset({0.1, 0.1}, {0.9, 0.9});

  CHECK_THAT(r[0], WithinAbs(0.2, 1e-15));  // wrapped
  CHECK_THAT(r[1], WithinAbs(-0.8, 1e-15)); // not wrapped
}

TEST_CASE("Periodic: half a box apart is still antisymmetric", "[domain][periodic]")
{
  // Both images are equidistant, so which one wins is a tie broken by std::round.
  // Pinned here so a change of rounding is caught. Note this separation sits on
  // the limit that admits_bandwidth rejects.
  Periodic<double, 1> g{{2.0}};
  const Point<double, 1> x{0.0}, y{1.0};

  CHECK(g.offset(x, y)[0] == 1.0);
  CHECK(g.offset(y, x)[0] == -1.0);
  CHECK(g.distance_squared(x, y) == 1.0);
}

TEST_CASE("Periodic: points far outside the box still wrap", "[domain][periodic]")
{
  Periodic<double, 1> g{{1.0}};

  CHECK_THAT(g.offset({5.1}, {0.9})[0], WithinAbs(0.2, 1e-14));
  CHECK_THAT(g.offset({-3.1}, {0.9})[0], WithinAbs(0.0, 1e-14));
}

TEST_CASE("Periodic: negative lengths mean no wrapping", "[domain][periodic]")
{
  Periodic<double, 1> g{{-1.0}};
  const Point<double, 1> a{0.1}, b{0.9};

  CHECK(g.offset(a, b)[0] == -0.8);
  CHECK(g.max_bandwidth() == std::numeric_limits<double>::infinity());
}

TEST_CASE("Periodic: max_bandwidth is the smallest half period", "[domain][periodic]")
{
  CHECK(Periodic<double, 3>{{2.0, 4.0, 6.0}}.max_bandwidth() == 1.0);
  CHECK(Periodic<double, 3>{{0.0, 4.0, 6.0}}.max_bandwidth() == 2.0); // ignores the open axis
  CHECK(Periodic<double, 3>{}.max_bandwidth() == std::numeric_limits<double>::infinity());
}

TEST_CASE("Periodic: bandwidth guard is strict at the limit", "[domain][periodic][bandwidth]")
{
  Periodic<double, 3> g{{2.0, 4.0, 0.0}};
  REQUIRE(g.max_bandwidth() == 1.0);

  CHECK(admits_bandwidth(g, 0.999));
  CHECK_FALSE(admits_bandwidth(g, 1.0)); // equal to the limit is not admissible
  CHECK_FALSE(admits_bandwidth(g, 2.0));
}

TEST_CASE("Periodic: works at compile time", "[domain][periodic]")
{
  constexpr Periodic<double, 2> g{{2.0, 2.0}};
  static_assert(g.max_bandwidth() == 1.0);
  static_assert(admits_bandwidth(g, 0.5));
  // static_assert(g.distance_squared({0.25, 0.0}, {1.75, 0.0}) == 0.25); // not on clang
  SUCCEED();
}
