#!/bin/sh
# Configure, build and test CMake presets; stop at the first failure.
# Usage: tools/build-test [preset...]   (default: every preset below)
set -eu

cd "$(dirname "$0")/.."

[ "$#" -gt 0 ] || set -- dev-clang dev-gcc asan-ubsan

for preset in "$@"; do
  printf '\n=== %s ===\n' "$preset"
  cmake --preset "$preset"
  cmake --build --preset "$preset"
  ctest --preset "$preset" --output-on-failure --no-tests=error
done

printf '\n=== passed: %s ===\n' "$*"
