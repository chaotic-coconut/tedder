# tedder

[![CI](https://github.com/chaotic-coconut/tedder/actions/workflows/ci.yml/badge.svg)](https://github.com/chaotic-coconut/tedder/actions/workflows/ci.yml)

`tedder` is a header-only C++23 library for reconstructing scalar and vector fields from scattered samples with local polynomial regression. The aim is to estimate both the field value and its spatial derivatives.

> [!WARNING]
> This project is pre-alpha. Reconstruction is not implemented yet, and the API is still taking shape.

## Status

Pre-alpha. Reconstruction is not implemented.

The foundation types are in place and covered by tests. Stage 00 fixed
three defects in them: the `Domain` concept accepted types that were not
domains, non-finite periods produced NaN distances, and periodic reduction
lost the minimum image for some large finite coordinates.

Present:

* fixed-size points, values, and matrices
* non-owning views of points and samples
* Euclidean and per-axis periodic domains
* distance and bandwidth calculations
* `LocalFit`, the result type for a reconstruction

Next: compactly supported kernels, then brute-force neighbour search.

## Requirements

Building tedder requires a C++23 compiler and CMake 3.24 or newer. The commands below also assume Ninja.

CI covers GCC 13 and Clang 18. Catch2 is downloaded when the tests are configured.

`float` and `double` are tested. Other floating-point types are not
rejected but are not supported.

## Build and test

```sh
cmake --preset dev-clang
cmake --build --preset dev-clang
ctest --preset dev-clang --output-on-failure
```

## Use with CMake

```cmake
add_subdirectory(path/to/tedder)
target_link_libraries(your_target PRIVATE tedder::tedder)
```

The current API already supports mixed periodic and open geometry:

```cpp
#include <tedder/field.hpp>

using Point = tedder::Point<double, 3>;

// x and y have period 2; z remains open.
const tedder::Periodic<double, 3> domain{{2.0, 2.0, 0.0}};

const Point p{0.1, 0.0, 1.0};
const Point q{1.9, 0.0, 4.0};

const double d = tedder::distance(domain, p, q); // approximately 3.0067

const double h_max = domain.max_bandwidth();     // 1.0
const bool valid = tedder::admits_bandwidth(domain, 0.5);
```

Here the separation along `x` is `0.2`, because the shorter path crosses the periodic boundary. A bandwidth of exactly `1.0` is not admitted.

## Name

A tedder spreads cut hay over its surroundings. A reconstruction kernel does something similar with the influence of each sample.

## License

[MIT](LICENSE) © 2026 Viacheslav Kruglov.
