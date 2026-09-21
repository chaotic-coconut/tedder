// Systematic coverage of the Domain concept's requirements (F1). Each
// negative fixture changes exactly one aspect of a minimal valid domain and
// must make Domain<Fixture> evaluate to a clean `false` via static_assert,
// never a hard template error. Fixtures live at namespace scope: a local
// class cannot have a `static constexpr` member. distance_squared is not part of
// the concept (it is a free function over offset), so no fixture declares it
// unless the test is about exactly that.
//
// Every fixture below was probed on this project's GCC/Clang toolchains and
// confirmed to compile cleanly (no hard errors) before being kept here.

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

#include <tedder/field.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace tedder;

namespace
{
// ---------------------------------------------------------------- positive controls

struct ValidCustom
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;

    point_type offset(const point_type &p, const point_type &q) const noexcept
    {
      point_type r{};
      for (std::size_t i = 0; i < 2; ++i)
        r[i] = p[i] - q[i];
      return r;
    }
    double max_bandwidth() const noexcept { return std::numeric_limits<double>::infinity(); }
};

// A valid domain whose member distance_squared lies. The free functions must
// derive everything from offset and ignore it.
struct ValidWithBogusMember
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;

    point_type offset(const point_type &p, const point_type &q) const noexcept { return {p[0] - q[0], p[1] - q[1]}; }
    double distance_squared(const point_type &, const point_type &) const noexcept { return 999.0; }
    double max_bandwidth() const noexcept { return std::numeric_limits<double>::infinity(); }
};

static_assert(noexcept(distance(std::declval<const ValidCustom &>(), std::declval<const Point<double, 2> &>(),
                                std::declval<const Point<double, 2> &>())));
static_assert(noexcept(admits_bandwidth(std::declval<const ValidCustom &>(), std::declval<double>())));

// Point<T,D> is an alias for std::array<T,D>: spelling it out names the same
// type, so this must also satisfy Domain. Not a loophole.
struct ValidArraySpelling
{
    using scalar_type = double;
    using point_type = std::array<double, 2>;
    static constexpr std::size_t dimension = 2;

    point_type offset(const point_type &p, const point_type &q) const noexcept
    {
      point_type r{};
      for (std::size_t i = 0; i < 2; ++i)
        r[i] = p[i] - q[i];
      return r;
    }
    double max_bandwidth() const noexcept { return std::numeric_limits<double>::infinity(); }
};

// ---------------------------------------------------------------- scalar_type defects

struct StringScalar
{
    using scalar_type = std::string;
    using point_type = Point<double, 2>; // irrelevant: Real<scalar_type> fails first
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() const noexcept { return 0.0; }
};

struct IntScalar
{
    using scalar_type = int;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() const noexcept { return 0.0; }
};

// ---------------------------------------------------------------- dimension defects

// A genuinely non-constant dimension: `static int`, not `static constexpr`.
struct RuntimeDimension
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static int dimension; // defined out of line below, deliberately not constexpr

    point_type offset(const point_type &p, const point_type &q) const noexcept
    {
      point_type r{};
      for (std::size_t i = 0; i < 2; ++i)
        r[i] = p[i] - q[i];
      return r;
    }
    double max_bandwidth() const noexcept { return 0.0; }
};
int RuntimeDimension::dimension = 2;

// dimension == 0, otherwise fully self-consistent (point_type really is
// Point<double,0>) so the zero dimension is the sole defect.
struct ZeroDimension
{
    using scalar_type = double;
    using point_type = Point<double, 0>;
    static constexpr std::size_t dimension = 0;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() const noexcept { return 0.0; }
};

// A negative dimension stored as `int`, per the table's literal wording.
struct NegativeDimension
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr int dimension = -1;

    point_type offset(const point_type &p, const point_type &q) const noexcept
    {
      point_type r{};
      for (std::size_t i = 0; i < 2; ++i)
        r[i] = p[i] - q[i];
      return r;
    }
    double max_bandwidth() const noexcept { return 0.0; }
};

// ---------------------------------------------------------------- point_type defects

// dimension says 2, point_type is Point<T,3>.
struct PointTypeTooWide
{
    using scalar_type = double;
    using point_type = Point<double, 3>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() const noexcept { return 0.0; }
};

// scalar_type is double, point_type uses float.
struct PointTypeWrongScalar
{
    using scalar_type = double;
    using point_type = Point<float, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() const noexcept { return 0.0; }
};

// A foreign point type unrelated to std::array, otherwise fully consistent.
struct ForeignPoint
{
    double x, y;
};

struct ForeignPointDomain
{
    using scalar_type = double;
    using point_type = ForeignPoint;
    static constexpr std::size_t dimension = 2;

    point_type offset(const point_type &p, const point_type &q) const noexcept { return point_type{p.x - q.x, p.y - q.y}; }
    double max_bandwidth() const noexcept { return 0.0; }
};

// ---------------------------------------------------------------- missing methods

struct MissingOffset
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    double max_bandwidth() const noexcept { return 0.0; }
};

struct MissingMaxBandwidth
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
};

// ---------------------------------------------------------------- wrong return type

struct OffsetWrongReturn
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    int offset(const point_type &, const point_type &) const noexcept { return 0; }
    double max_bandwidth() const noexcept { return 0.0; }
};

struct MaxBandwidthWrongReturn
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    int max_bandwidth() const noexcept { return 0; }
};

// ---------------------------------------------------------------- not const

struct OffsetNotConst
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) noexcept { return {}; } // missing const
    double max_bandwidth() const noexcept { return 0.0; }
};

struct MaxBandwidthNotConst
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() noexcept { return 0.0; } // missing const
};

// ---------------------------------------------------------------- not noexcept

struct OffsetNotNoexcept
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const { return {}; } // missing noexcept
    double max_bandwidth() const noexcept { return 0.0; }
};

struct MaxBandwidthNotNoexcept
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;
    point_type offset(const point_type &, const point_type &) const noexcept { return {}; }
    double max_bandwidth() const { return 0.0; } // missing noexcept
};
// A dimension that is a class type with a constexpr conversion to size_t.
// Each has an operator> that disagrees with the converted value: these
// fixtures only pass if the concept compares the converted value, never
// operator> itself (which need not even be constexpr).
struct DimTwo
{
    constexpr operator std::size_t() const noexcept { return 2; }
    bool operator>(int) const noexcept { return false; } // lies, and not constexpr
};
struct DimZero
{
    constexpr operator std::size_t() const noexcept { return 0; }
    bool operator>(int) const noexcept { return true; } // lies the other way
};
struct ClassDimTwo : Euclidean<double, 2>
{
    static constexpr DimTwo dimension{};
};
struct ClassDimZero : Euclidean<double, 2>
{
    static constexpr DimZero dimension{};
};
} // namespace

static_assert(Domain<ClassDimTwo>);   // converted value 2 is positive
static_assert(!Domain<ClassDimZero>); // converted value 0, despite operator>

// ---------------------------------------------------------------- verdicts

static_assert(Domain<ValidCustom>);
static_assert(Domain<ValidArraySpelling>);
static_assert(Domain<Euclidean<double, 2>>);
static_assert(Domain<Periodic<float, 3>>);

static_assert(!Domain<StringScalar>);
static_assert(!Domain<IntScalar>);

static_assert(!Domain<RuntimeDimension>);
static_assert(!Domain<ZeroDimension>);
static_assert(!Domain<NegativeDimension>);

static_assert(!Domain<PointTypeTooWide>);
static_assert(!Domain<PointTypeWrongScalar>);
static_assert(!Domain<ForeignPointDomain>);

static_assert(!Domain<MissingOffset>);
static_assert(!Domain<MissingMaxBandwidth>);

static_assert(!Domain<OffsetWrongReturn>);
static_assert(!Domain<MaxBandwidthWrongReturn>);

static_assert(!Domain<OffsetNotConst>);
static_assert(!Domain<MaxBandwidthNotConst>);

static_assert(!Domain<OffsetNotNoexcept>);
static_assert(!Domain<MaxBandwidthNotNoexcept>);

static_assert(!Domain<int>);

namespace other
{
// Own namespace, so an unqualified distance_squared call in the library would
// find the overload below by argument-dependent lookup.
struct OwnDomain
{
    using scalar_type = double;
    using point_type = Point<double, 2>;
    static constexpr std::size_t dimension = 2;

    point_type offset(const point_type &p, const point_type &q) const noexcept { return {p[0] - q[0], p[1] - q[1]}; }
    double max_bandwidth() const noexcept { return std::numeric_limits<double>::infinity(); }
};

double distance_squared(const OwnDomain &, const OwnDomain::point_type &, const OwnDomain::point_type &) noexcept { return 999.0; }
} // namespace other

static_assert(Domain<ValidWithBogusMember>);
static_assert(Domain<other::OwnDomain>);

// A trivial runtime presence so this translation unit contributes a
// TEST_CASE to the suite (the checks above are all compile-time).
TEST_CASE("Domain: concept fixtures compiled and evaluated as expected", "[domain][contract]") { SUCCEED(); }

TEST_CASE("Domain: a custom Domain works through the free helper functions", "[domain][contract]")
{
  const ValidCustom g;
  const Point<double, 2> p{3.0, 4.0}, q{0.0, 0.0};

  CHECK(distance(g, p, q) == 5.0);
  CHECK(admits_bandwidth(g, 1.0));
  CHECK_FALSE(admits_bandwidth(g, -1.0));
}

TEST_CASE("Domain: a member distance_squared is ignored by the free functions", "[domain][contract]")
{
  const ValidWithBogusMember g;
  const Point<double, 2> p{3.0, 4.0}, q{0.0, 0.0};

  CHECK(tedder::distance_squared(g, p, q) == 25.0);
  CHECK(tedder::distance(g, p, q) == 5.0);
}

TEST_CASE("Domain: distance ignores an ADL-visible distance_squared overload", "[domain][contract]")
{
  const other::OwnDomain g;
  const Point<double, 2> p{6.0, 8.0}, q{0.0, 0.0};

  CHECK(tedder::distance_squared(g, p, q) == 100.0);
  CHECK(tedder::distance(g, p, q) == 10.0);
}
