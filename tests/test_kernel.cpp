// Tests for the Kernel concept and the Wendland kernels: T1 and T2 of the
// Stage 01 plan (docs/stages/01-kernels/kernel-plan.md). The contract they
// check is ADR 0003 §4.

#include <array>
#include <cmath>
#include <concepts>
#include <limits>

#include <tedder/kernel.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace tedder;

// ---------------------------------------------------------------- T1: concept

namespace
{
  // Every fixture declares its return type. With a deduced return type, checking
  // the concept would instantiate the body, and an ill-formed body would be a
  // hard error instead of a clean false.
  struct NoCall
  {
  };

  struct NotNoexcept
  {
    double operator()(double q) const { return q; }
  };

  // A float kernel that computes and returns double.
  struct ReturnsDouble
  {
    double operator()(float q) const noexcept { return 1.0 - q; }
  };

  struct NonConstCall
  {
    double operator()(double q) noexcept { return q; }
  };

  // Could modify the caller's q. The concept calls with a const T, so it is
  // rejected.
  struct TakesReference
  {
    double operator()(double &q) const noexcept { return q; }
  };

  struct TakesConstReference
  {
    double operator()(const double &q) const noexcept { return q; }
  };

  using NoexceptLambda = decltype([](double q) noexcept { return q; });
  using PlainLambda = decltype([](double q) { return q; });
} // namespace

static_assert(Kernel<WendlandC2, float> && Kernel<WendlandC2, double>);
static_assert(Kernel<WendlandC4, float> && Kernel<WendlandC4, double>);
static_assert(Kernel<NoexceptLambda, double>);
static_assert(Kernel<TakesConstReference, double>);

static_assert(!Kernel<WendlandC2, int>);
static_assert(!Kernel<NoCall, double>);
static_assert(!Kernel<NotNoexcept, double>);
static_assert(!Kernel<ReturnsDouble, float>);
static_assert(!Kernel<NonConstCall, double>);
static_assert(!Kernel<TakesReference, double>);
static_assert(!Kernel<PlainLambda, double>);
// A double lambda called with a float returns double, not float.
static_assert(!Kernel<NoexceptLambda, float>);

// The kernels are usable in constant expressions, and their ends are exact.
static_assert(WendlandC2{}(0.0f) == 1.0f && WendlandC2{}(1.0f) == 0.0f);
static_assert(WendlandC2{}(0.0) == 1.0 && WendlandC2{}(1.0) == 0.0);
static_assert(WendlandC4{}(0.0f) == 1.0f && WendlandC4{}(1.0f) == 0.0f);
static_assert(WendlandC4{}(0.0) == 1.0 && WendlandC4{}(1.0) == 0.0);

// ---------------------------------------------------------------- T2: values

namespace
{
  // Calls a kernel the way the fitter will: through a const reference, under
  // the concept.
  template <class K, Real T>
    requires Kernel<K, T>
  T weight(const K &k, T q) noexcept
  {
    return k(q);
  }

  struct Reference
  {
    double q;        // converted to T before the call
    double w_float;  // exact weight at static_cast<float>(q), rounded to double
    double w_double; // exact weight at q, rounded to double
  };

  // Computed offline with Python's fractions.Fraction from the expanded
  // polynomials, not from the factored forms in kernel.hpp:
  //   C2: 1 - 10q^2 + 20q^3 - 15q^4 + 4q^5
  //   C4: 1 - (28/3)q^2 + 70q^4 - (448/3)q^5 + 140q^6 - 64q^7 + (35/3)q^8
  // The two columns differ where the float and double values of q differ.
  // Every weight is above min() / epsilon() for float, so relative accuracy
  // applies (plan, "Relative accuracy").
  constexpr std::array<Reference, 13> c2_reference{{
      {0.1, 0.9185399978274107, 0.91854},
      {0.125, 0.8792724609375, 0.8792724609375},
      {0.25, 0.6328125, 0.6328125},
      {0.3, 0.5282199754667284, 0.52822},
      {0.3333333333333333, 0.4609053301713105, 0.4609053497942387},
      {0.45, 0.25621751785010133, 0.2562175},
      {0.5, 0.1875, 0.1875},
      {0.6, 0.0870399816894544, 0.08704000000000002},
      {0.75, 0.015625, 0.015625},
      {0.875, 0.0010986328125, 0.0010986328125},
      {0.9, 0.0004600004291535902, 0.0004599999999999996},
      {0.99, 4.9599811172754563e-08, 4.960000000000017e-08},
      {0.9990234375, 4.543920795185841e-12, 4.543920795185841e-12},
  }};

  constexpr std::array<Reference, 13> c4_reference{{
      {0.1, 0.9123070475362837, 0.91230705},
      {0.125, 0.8672034541765848, 0.8672034541765848},
      {0.25, 0.5747222900390625, 0.5747222900390625},
      {0.3, 0.4529486219502929, 0.45294865},
      {0.3333333333333333, 0.3771782539630093, 0.3771782756693594},
      {0.45, 0.1678139001680182, 0.16781388378906248},
      {0.5, 0.10807291666666667, 0.10807291666666667},
      {0.6, 0.03604478906250125, 0.03604480000000001},
      {0.75, 0.0029449462890625, 0.0029449462890625},
      {0.875, 5.7915846506754555e-05, 5.7915846506754555e-05},
      {0.9, 1.585002202988931e-05, 1.584999999999998e-05},
      {0.99, 1.837439513803595e-11, 1.8374500000000097e-11},
      {0.9990234375, 1.6165915793104553e-17, 1.6165915793104553e-17},
  }};

  // Relative tolerance n * epsilon of T, from the plan's first-order bounds
  // with u = epsilon / 2. For q < 1/2, t = 1 - q carries one rounding.
  // - C2: t^4 has 4 times that plus 3, (4q + 1) and the final product one each:
  //   9u. For q >= 1/2, t is exact: 5u.
  // - C4: t^6 has 6 times that plus 5, the quadratic 5 (the rounded 35/3, two
  //   multiplications and two additions of positive terms), the product one:
  //   17u, or 11u for q >= 1/2.
  // The comparison is done in double. For double the reference adds half an
  // epsilon of rounding; for float it adds 2^-53, which is negligible.
  // n = 4.5 + 0.5 = 5 for C2 and 8.5 + 0.5 = 9 for C4.
  template <class K, Real T>
  void check_reference(const std::array<Reference, 13> &table, double n)
  {
    const double eps = static_cast<double>(std::numeric_limits<T>::epsilon());
    for (const Reference &r : table)
    {
      const T q = static_cast<T>(r.q);
      const double expected = std::same_as<T, float> ? r.w_float : r.w_double;
      const double got = static_cast<double>(weight(K{}, q));
      CAPTURE(r.q, got, expected);
      CHECK(std::abs(got - expected) <= n * eps * expected);
    }
  }

  template <class K, Real T>
  void check_ends()
  {
    CHECK(weight(K{}, T{0}) == T{1});
    CHECK(weight(K{}, -T{0}) == T{1}); // -0 satisfies 0 <= q
    CHECK(weight(K{}, T{1}) == T{0});
  }

  // Finite and non-negative on [0, 1]; positive on the grid points of [0, 1),
  // whose exact weights are all normal numbers.
  template <class K, Real T>
  void check_grid()
  {
    constexpr int steps = 1000;
    for (int i = 0; i <= steps; ++i)
    {
      const T q = static_cast<T>(i) / static_cast<T>(steps);
      const T w = weight(K{}, q);
      CAPTURE(q, w);
      CHECK(std::isfinite(w));
      CHECK(w >= T{0});
      if (i < steps)
        CHECK(w > T{0});
    }
  }

  // The 2,000 largest values below 1. An expanded polynomial cancels here and
  // returns negative weights; the factored forms must not. The weight may be
  // subnormal or 0 (C4 in float), so only finite and >= 0 are required.
  template <class K, Real T>
  void check_near_one()
  {
    T q = T{1};
    for (int i = 0; i < 2000; ++i)
    {
      q = std::nextafter(q, T{0});
      const T w = weight(K{}, q);
      CAPTURE(q, w);
      CHECK(std::isfinite(w));
      CHECK(w >= T{0});
    }
  }
} // namespace

TEMPLATE_TEST_CASE("Kernel: w(0) == 1 and w(1) == 0 exactly", "[kernel]", float, double)
{
  check_ends<WendlandC2, TestType>();
  check_ends<WendlandC4, TestType>();
}

TEMPLATE_TEST_CASE("Kernel: weights match exact reference values", "[kernel]", float, double)
{
  check_reference<WendlandC2, TestType>(c2_reference, 5.0);
  check_reference<WendlandC4, TestType>(c4_reference, 9.0);
}

TEMPLATE_TEST_CASE("Kernel: weights are finite, non-negative and positive inside", "[kernel]", float, double)
{
  check_grid<WendlandC2, TestType>();
  check_grid<WendlandC4, TestType>();
}

TEMPLATE_TEST_CASE("Kernel: weights stay non-negative just below q = 1", "[kernel]", float, double)
{
  check_near_one<WendlandC2, TestType>();
  check_near_one<WendlandC4, TestType>();
}
