# tedder

[![CI](https://github.com/chaotic-coconut/tedder/actions/workflows/ci.yml/badge.svg)](https://github.com/chaotic-coconut/tedder/actions/workflows/ci.yml)

`tedder` is an experimental, header-only C++23 library for reconstructing scalar and vector fields from scattered samples using local polynomial regression.

> [!WARNING]
> `tedder` is pre-alpha. Reconstruction is not implemented yet, and the API may change.

## Status

Currently implemented:

* fixed-size points, values, matrices, and non-owning sample views
* Euclidean and per-axis periodic domains
* distance and bandwidth calculations
* the `LocalFit` result type

Planned:

* compactly supported kernels
* brute-force and kd-tree neighbour search
* weighted polynomial fits of degree 0–2
* Jacobians, divergence, vorticity, and strain
* NumPy `.npy` input and output

## Requirements

* C++23 compiler
* CMake 3.24 or newer
* Ninja for the commands below

CI currently tests GCC 13 and Clang 18. Catch2 is downloaded automatically when building the tests.

## Build and test

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

## Use with CMake

```cmake
add_subdirectory(path/to/tedder)
target_link_libraries(your_target PRIVATE tedder::tedder)
```

Example using the current geometry API:

```cpp
#include <tedder/field.hpp>

using Point = tedder::Point<double, 3>;

// Periodic in x and y with period 2, open in z.
const tedder::Periodic<double, 3> domain{{2.0, 2.0, 0.0}};

const Point p{0.1, 0.0, 1.0};
const Point q{1.9, 0.0, 4.0};

// x wraps, so 1.8 apart the long way is 0.2 the short way. z does not wrap.
const double d = tedder::distance(domain, p, q);   // 3.0067, not 3.4986

// Half the shortest period: the largest bandwidth the domain allows.
const double h_max = domain.max_bandwidth();       // 1.0
const bool ok = tedder::admits_bandwidth(domain, 0.5);   // true, and 1.0 is false
```

## Name

A tedder spreads cut hay over a neighbourhood. The name reflects how a kernel spreads each sample’s influence.

## License

[MIT](LICENSE) © 2026 Viacheslav Kruglov.
