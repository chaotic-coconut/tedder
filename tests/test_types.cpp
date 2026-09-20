#include <array>
#include <concepts>
#include <cstddef>
#include <utility>

#include <tedder/field.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace tedder;

// ---------------------------------------------------------------- layout invariants
//
// Point<T,D> is std::array<T,D>: no padding, no vtable, no hidden state.
// Mirrors the type/dimension combinations tests/test_domain.cpp exercises.

static_assert(sizeof(Point<double, 1>) == 1 * sizeof(double));
static_assert(sizeof(Point<double, 2>) == 2 * sizeof(double));
static_assert(sizeof(Point<double, 3>) == 3 * sizeof(double));
static_assert(sizeof(Point<float, 2>) == 2 * sizeof(float));
static_assert(sizeof(Point<float, 3>) == 3 * sizeof(float));

static_assert(is_flat<double, 1>);
static_assert(is_flat<double, 2>);
static_assert(is_flat<double, 3>);
static_assert(is_flat<float, 2>);
static_assert(is_flat<float, 3>);

// ---------------------------------------------------------------- Matrix

TEST_CASE("Matrix: rectangular storage is row-major", "[types][matrix]")
{
  Matrix<double, 2, 3> m;

  // Six distinct entries, unrelated to their (row, col) indices, so the
  // comparison below cannot pass by coincidentally matching the storage
  // formula it is meant to check.
  m[0, 0] = 2.5;
  m[0, 1] = -1.0;
  m[0, 2] = 100.25;
  m[1, 0] = 7.75;
  m[1, 1] = -3.5;
  m[1, 2] = 42.0;

  const double expected[6] = {2.5, -1.0, 100.25, 7.75, -3.5, 42.0};
  for (std::size_t i = 0; i < 6; ++i)
    CHECK(m.data()[i] == expected[i]);

  const Matrix<double, 2, 3> &cm = m;
  CHECK(cm[0, 0] == 2.5);
  CHECK(cm[0, 2] == 100.25);
  CHECK(cm[1, 1] == -3.5);
  CHECK(cm[1, 2] == 42.0);
}

TEST_CASE("Matrix: rectangular storage is row-major (transposed shape)", "[types][matrix]")
{
  Matrix<double, 3, 2> m;

  m[0, 0] = 5.5;
  m[0, 1] = -2.25;
  m[1, 0] = 13.0;
  m[1, 1] = 0.125;
  m[2, 0] = -8.75;
  m[2, 1] = 64.0;

  const double expected[6] = {5.5, -2.25, 13.0, 0.125, -8.75, 64.0};
  for (std::size_t i = 0; i < 6; ++i)
    CHECK(m.data()[i] == expected[i]);

  const Matrix<double, 3, 2> &cm = m;
  CHECK(cm[0, 1] == -2.25);
  CHECK(cm[2, 0] == -8.75);
}

TEST_CASE("Matrix: writing one cell leaves every other cell unchanged", "[types][matrix]")
{
  Matrix<double, 2, 3> m;
  m[0, 0] = 1.0;
  m[0, 1] = 2.0;
  m[0, 2] = 3.0;
  m[1, 0] = 4.0;
  m[1, 1] = 5.0;
  m[1, 2] = 6.0;

  m[1, 1] = 99.0; // overwrite a single cell

  CHECK(m[0, 0] == 1.0);
  CHECK(m[0, 1] == 2.0);
  CHECK(m[0, 2] == 3.0);
  CHECK(m[1, 0] == 4.0);
  CHECK(m[1, 1] == 99.0);
  CHECK(m[1, 2] == 6.0);
}

TEST_CASE("Matrix: float storage is row-major", "[types][matrix]")
{
  Matrix<float, 2, 2> m;
  m[0, 0] = 1.5f;
  m[0, 1] = -2.25f;
  m[1, 0] = 4.0f;
  m[1, 1] = 0.0625f;

  const float expected[4] = {1.5f, -2.25f, 4.0f, 0.0625f};
  for (std::size_t i = 0; i < 4; ++i)
    CHECK(m.data()[i] == expected[i]);
}

TEST_CASE("Matrix: works at compile time", "[types][matrix]")
{
  constexpr Matrix<double, 2, 2> m{{1.0, 2.0, 3.0, 4.0}};
  static_assert(m[0, 1] == 2.0);
  static_assert(m[1, 0] == 3.0);
  SUCCEED();
}

static_assert(sizeof(Matrix<double, 2, 3>) == 2 * 3 * sizeof(double));
static_assert(sizeof(Matrix<double, 3, 2>) == 3 * 2 * sizeof(double));
static_assert(sizeof(Matrix<float, 2, 2>) == 2 * 2 * sizeof(float));
static_assert(sizeof(Matrix<double, 2, 2>) == 2 * 2 * sizeof(double));

static_assert(Matrix<double, 2, 3>::rows == 2 && Matrix<double, 2, 3>::cols == 3);
static_assert(Matrix<double, 3, 2>::rows == 3 && Matrix<double, 3, 2>::cols == 2);
static_assert(Matrix<float, 2, 2>::rows == 2 && Matrix<float, 2, 2>::cols == 2);
static_assert(Matrix<double, 2, 2>::rows == 2 && Matrix<double, 2, 2>::cols == 2);

// ---------------------------------------------------------------- LocalFit

TEST_CASE("LocalFit: rectangular Jacobian keeps D and C distinct", "[types][localfit]")
{
  // D=2, C=3: jacobian is Matrix<T,C,D>, so a swapped declaration
  // (Matrix<T,D,C>) would give 2 rows and 3 columns instead.
  using Fit = LocalFit<double, 2, 3>;
  using Jacobian = decltype(std::declval<Fit>().jacobian);

  static_assert(Jacobian::rows == 3);
  static_assert(Jacobian::cols == 2);

  Fit fit;
  fit.jacobian[2, 1] = 9.5;
  CHECK(fit.jacobian[2, 1] == 9.5);
  CHECK(std::as_const(fit).jacobian[2, 1] == 9.5);
}

namespace
{
// jacobian is Matrix<T,C,D>: C rows, D columns. Filling every entry
// independently is the regression check for a transposed declaration that
// read or wrote past the end of a smaller dimension (a heap overread when
// C > D).
template <std::size_t D, std::size_t C> void check_localfit_jacobian_shape_and_fill()
{
  using Fit = LocalFit<double, D, C>;
  using Jacobian = decltype(std::declval<Fit>().jacobian);

  static_assert(Jacobian::rows == C);
  static_assert(Jacobian::cols == D);

  Fit fit;
  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t d = 0; d < D; ++d)
      fit.jacobian[c, d] = static_cast<double>(c * 1000 + d) + 1.0;

  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t d = 0; d < D; ++d)
      CHECK(fit.jacobian[c, d] == static_cast<double>(c * 1000 + d) + 1.0);
}

template <std::size_t D, std::size_t C> void check_localfit_defaults()
{
  LocalFit<double, D, C> fit;

  for (std::size_t c = 0; c < C; ++c)
    CHECK(fit.value[c] == 0.0);

  for (std::size_t c = 0; c < C; ++c)
    for (std::size_t d = 0; d < D; ++d)
      CHECK(fit.jacobian[c, d] == 0.0);

  CHECK(fit.neighbors == 0);
  CHECK(fit.bandwidth == 0.0);
}
} // namespace

TEST_CASE("LocalFit: every jacobian entry is independently addressable", "[types][localfit]")
{
  check_localfit_jacobian_shape_and_fill<3, 1>(); // C < D
  check_localfit_jacobian_shape_and_fill<2, 3>(); // C > D, the historical overread case
  check_localfit_jacobian_shape_and_fill<3, 3>(); // C == D
}

TEST_CASE("LocalFit: value-initialised defaults are all zero", "[types][localfit]")
{
  check_localfit_defaults<3, 1>();
  check_localfit_defaults<2, 3>();
  check_localfit_defaults<3, 3>();
}

// ---------------------------------------------------------------- SampleView

TEST_CASE("SampleView: subview selects the right pair and stays aligned", "[types][sampleview]")
{
  const std::array<Point<double, 2>, 4> points{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
  const std::array<Value<double, 1>, 4> values{{{10.0}, {11.0}, {12.0}, {13.0}}};

  const SampleView<double, 2, 1> view(points, values);
  REQUIRE(view.size() == 4);

  const auto sub = view.subview(1, 2);

  REQUIRE(sub.size() == 2);
  CHECK(sub.point(0) == points[1]);
  CHECK(sub.point(1) == points[2]);
  CHECK(sub.value(0) == values[1]);
  CHECK(sub.value(1) == values[2]);
}

TEST_CASE("SampleView: size and emptiness", "[types][sampleview]")
{
  SECTION("empty")
  {
    const std::array<Point<double, 2>, 0> points{};
    const std::array<Value<double, 1>, 0> values{};
    const SampleView<double, 2, 1> view(points, values);
    CHECK(view.size() == 0);
    CHECK(view.empty());
  }
  SECTION("one sample")
  {
    const std::array<Point<double, 2>, 1> points{{{1.5, -2.5}}};
    const std::array<Value<double, 1>, 1> values{{{7.0}}};
    const SampleView<double, 2, 1> view(points, values);
    CHECK(view.size() == 1);
    CHECK_FALSE(view.empty());
    CHECK(view.point(0) == points[0]);
    CHECK(view.value(0) == values[0]);
  }
  SECTION("four samples")
  {
    const std::array<Point<double, 2>, 4> points{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
    const std::array<Value<double, 1>, 4> values{{{10.0}, {11.0}, {12.0}, {13.0}}};
    const SampleView<double, 2, 1> view(points, values);
    CHECK(view.size() == 4);
    CHECK_FALSE(view.empty());
  }
}

TEST_CASE("SampleView: points() and values() stay index-aligned", "[types][sampleview]")
{
  const std::array<Point<double, 2>, 4> points{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
  const std::array<Value<double, 1>, 4> values{{{10.0}, {11.0}, {12.0}, {13.0}}};
  const SampleView<double, 2, 1> view(points, values);

  const auto pts = view.points();
  const auto vals = view.values();
  REQUIRE(pts.size() == 4);
  REQUIRE(vals.size() == 4);
  for (std::size_t i = 0; i < 4; ++i)
    {
      CHECK(pts[i] == points[i]);
      CHECK(vals[i] == values[i]);
      CHECK(pts[i] == view.point(i));
      CHECK(vals[i] == view.value(i));
    }
}

TEST_CASE("SampleView: subview variants stay aligned", "[types][sampleview]")
{
  const std::array<Point<double, 2>, 4> points{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
  const std::array<Value<double, 1>, 4> values{{{10.0}, {11.0}, {12.0}, {13.0}}};
  const SampleView<double, 2, 1> view(points, values);

  auto check = [&](std::size_t first, std::size_t count)
  {
    const auto sub = view.subview(first, count);
    REQUIRE(sub.size() == count);
    for (std::size_t i = 0; i < count; ++i)
      {
        CHECK(sub.point(i) == points[first + i]);
        CHECK(sub.value(i) == values[first + i]);
      }
  };

  check(0, 4); // full
  check(0, 2); // prefix
  check(2, 2); // suffix
  check(1, 2); // middle window
  check(0, 0); // empty at the start
  check(4, 0); // empty at the end

  // Nested subview of a subview.
  const auto mid = view.subview(1, 3);   // samples 1,2,3
  const auto nested = mid.subview(1, 2); // samples 2,3 of the original
  REQUIRE(nested.size() == 2);
  CHECK(nested.point(0) == points[2]);
  CHECK(nested.point(1) == points[3]);
  CHECK(nested.value(0) == values[2]);
  CHECK(nested.value(1) == values[3]);
}

namespace
{
// Every generic function in the library will hold a sample set through a
// const reference; subview() must compile and work through one.
SampleView<double, 2, 1> take_subview_by_const_ref(const SampleView<double, 2, 1> &view, std::size_t first, std::size_t count)
{
  return view.subview(first, count);
}
} // namespace

TEST_CASE("SampleView: subview works through a const reference", "[types][sampleview]")
{
  const std::array<Point<double, 2>, 4> points{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
  const std::array<Value<double, 1>, 4> values{{{10.0}, {11.0}, {12.0}, {13.0}}};
  const SampleView<double, 2, 1> view(points, values);

  const auto sub = take_subview_by_const_ref(view, 1, 2);
  REQUIRE(sub.size() == 2);
  CHECK(sub.point(0) == points[1]);
  CHECK(sub.point(1) == points[2]);
  CHECK(sub.value(0) == values[1]);
  CHECK(sub.value(1) == values[2]);
}

TEST_CASE("SampleView: copies see the same elements", "[types][sampleview]")
{
  const std::array<Point<double, 2>, 4> points{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
  const std::array<Value<double, 1>, 4> values{{{10.0}, {11.0}, {12.0}, {13.0}}};
  const SampleView<double, 2, 1> view(points, values);

  SampleView<double, 2, 1> copy = view;
  REQUIRE(copy.size() == view.size());
  for (std::size_t i = 0; i < view.size(); ++i)
    {
      CHECK(copy.point(i) == view.point(i));
      CHECK(copy.value(i) == view.value(i));
    }
}

static_assert(SampleView<double, 2, 1>::dimension == 2);
static_assert(SampleView<double, 2, 1>::components == 1);
static_assert(std::same_as<SampleView<double, 2, 1>::scalar_type, double>);

namespace
{
// Namespace-scope storage so a constexpr fixture never references a local.
static constexpr std::array<Point<double, 2>, 4> g_cpoints{{{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}}};
static constexpr std::array<Value<double, 1>, 4> g_cvalues{{{10.0}, {11.0}, {12.0}, {13.0}}};

constexpr bool check_sampleview_constexpr()
{
  SampleView<double, 2, 1> view(g_cpoints, g_cvalues);
  if (view.size() != 4)
    return false;

  const auto sub = view.subview(1, 2);
  if (sub.size() != 2)
    return false;
  if (sub.point(0) != g_cpoints[1] || sub.point(1) != g_cpoints[2])
    return false;
  if (sub.value(0) != g_cvalues[1] || sub.value(1) != g_cvalues[2])
    return false;

  return true;
}
} // namespace

TEST_CASE("SampleView: works at compile time", "[types][sampleview]")
{
  static_assert(check_sampleview_constexpr());
  SUCCEED();
}
