#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace tedder
{
template <class T>
concept Real = std::floating_point<T>;

// ---------------------------------------------------------------- domain:

template <Real T, std::size_t D> using Point = std::array<T, D>;

template <Real T, std::size_t D> using PointCloud = std::vector<Point<T, D>>;

template <Real T, std::size_t D> using PointView = std::span<const Point<T, D>>;

// ---------------------------------------------------------------- codomain:

// C=1 scalar, C=D vector, C=2 complex-as-R2, C=D*D tensor
template <Real T, std::size_t C> using Value = std::array<T, C>;

template <Real T, std::size_t C> using ValueField = std::vector<Value<T, C>>;

template <Real T, std::size_t C> using ValueView = std::span<const Value<T, C>>;

// ---------------------------------------------------------------- geometry:

// Anything that can measure distance between points.
// Euclidean and Periodic both satisfy it.
// offset(p, q) is p - q. On a periodic axis, wrapped to the nearest copy.
template <class G>
concept Domain = requires(const G &g, const typename G::point_type &p, const typename G::point_type &q) {
  typename G::scalar_type;
  typename G::point_type;
  { G::dimension } -> std::convertible_to<std::size_t>;
  // offset(p, q) is p - q, wrapped to the nearest image on a periodic domain
  { g.offset(p, q) } -> std::same_as<typename G::point_type>;
  { g.distance_squared(p, q) } -> std::same_as<typename G::scalar_type>;
  { g.max_bandwidth() } -> std::same_as<typename G::scalar_type>;
};

template <Real T, std::size_t D> struct Euclidean
{
    using scalar_type = T;
    using point_type = Point<T, D>;
    static constexpr std::size_t dimension = D;

    constexpr point_type offset(const point_type &p, const point_type &q) const noexcept
    {
      point_type r{};
      for (std::size_t i = 0; i < D; ++i)
        r[i] = p[i] - q[i];
      return r;
    }

    constexpr scalar_type distance_squared(const point_type &p, const point_type &q) const noexcept
    {
      const point_type r = offset(p, q);
      scalar_type s{};
      for (std::size_t i = 0; i < D; ++i)
        s += r[i] * r[i];
      return s;
    }

    constexpr scalar_type max_bandwidth() const noexcept { return std::numeric_limits<scalar_type>::infinity(); }
};

// Box that wraps around on some axes.
// length[d] <= 0 means axis d does not wrap.
// Default (all zeros) means nothing wraps, same as Euclidean.
template <Real T, std::size_t D> struct Periodic
{
    using scalar_type = T;
    using point_type = Point<T, D>;
    static constexpr std::size_t dimension = D;

    // length[d] <= 0 means axis d does not wrap.
    // All zeros (the default) is therefore plain Euclidean, not NaN.
    point_type length{};

    constexpr point_type offset(const point_type &p, const point_type &q) const noexcept
    {
      point_type r{};
      for (std::size_t d = 0; d < D; ++d)
        {
          scalar_type t = p[d] - q[d];
          if (length[d] > scalar_type{})
            t -= length[d] * std::round(t / length[d]);
          r[d] = t;
        }
      return r;
    }

    constexpr scalar_type distance_squared(const point_type &p, const point_type &q) const noexcept
    {
      const point_type r = offset(p, q);
      scalar_type s{};
      for (std::size_t i = 0; i < D; ++i)
        s += r[i] * r[i];
      return s;
    }

    // Largest bandwidth that still works, as a strict limit.
    // h must be smaller than this, not equal to it.
    // Returns infinity when no axis wraps.
    constexpr scalar_type max_bandwidth() const noexcept
    {
      scalar_type m = std::numeric_limits<scalar_type>::infinity();
      for (std::size_t d = 0; d < D; ++d)
        if (length[d] > scalar_type{})
          m = std::min(m, length[d] / 2);
      return m;
    }
};

template <Domain G>
constexpr typename G::scalar_type distance(const G &g, const typename G::point_type &p, const typename G::point_type &q) noexcept
{
  return std::sqrt(g.distance_squared(p, q));
}

// True if h is a usable bandwidth for this domain.
// Also rejects zero, negative and NaN.
template <Domain G> constexpr bool admits_bandwidth(const G &g, typename G::scalar_type h) noexcept
{
  return h > typename G::scalar_type{} && h < g.max_bandwidth();
}

// ---------------------------------------------------------------- bundles

// Points and their values, side by side. Owns nothing.
// Both spans must outlive this object and must be the same length.
template <Real T, std::size_t D, std::size_t C> class SampleView
{
  public:
    using scalar_type = T;
    static constexpr std::size_t dimension = D, components = C;

    // precondition: pts.size() == vals.size()
    constexpr SampleView(PointView<T, D> pts, ValueView<T, C> vals) noexcept : points_(pts), values_(vals)
    {
      assert(pts.size() == vals.size() && "SampleView: length mismatch");
    }

    constexpr std::size_t size() const noexcept { return points_.size(); }
    constexpr bool empty() const noexcept { return points_.empty(); }

    constexpr PointView<T, D> points() const noexcept { return points_; }
    constexpr ValueView<T, C> values() const noexcept { return values_; }

    constexpr const Point<T, D> &point(std::size_t i) const noexcept { return points_[i]; }
    constexpr const Value<T, C> &value(std::size_t i) const noexcept { return values_[i]; }

    constexpr SampleView subview(std::size_t first, std::size_t count) const noexcept
    {
      return SampleView(points_.subspan(first, count), values_.subspan(first, count));
    }

  private:
    PointView<T, D> points_;
    ValueView<T, C> values_;
};

// Small fixed-size matrix. Stored flat, row by row.
template <Real T, std::size_t R, std::size_t C> struct Matrix
{
    static_assert(R > 0 && C > 0);
    std::array<T, R * C> m{}; // flat, row-major
    static constexpr std::size_t rows = R, cols = C;

    constexpr T &operator[](std::size_t r, std::size_t c) noexcept { return m[r * cols + c]; }
    constexpr const T &operator[](std::size_t r, std::size_t c) const noexcept { return m[r * cols + c]; }
    constexpr T *data() noexcept { return m.data(); }
    constexpr const T *data() const noexcept { return m.data(); }
};

// Result of one local fit. Determined entirely by D and C.
// jacobian[c, d] is the derivative of component c along axis d.
template <Real T, std::size_t D, std::size_t C> struct LocalFit
{
    Value<T, C> value{};        // f(x)
    Matrix<T, C, D> jacobian{}; // J[c][d] = d f_c / d x_d
    std::size_t neighbors = 0;
    T bandwidth{};
};

// ---------------------------------------------------------------- invariants

// True if array<T,N> has no padding, so it can be written straight to a file.
template <Real T, std::size_t N> inline constexpr bool is_flat = (sizeof(std::array<T, N>) == N * sizeof(T));

} // namespace tedder
