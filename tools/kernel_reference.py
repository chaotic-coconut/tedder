#!/usr/bin/env python3
"""Print the reference weights used in tests/test_kernel.cpp.

The Wendland kernels are evaluated in exact rational arithmetic, from their
expanded polynomials rather than the factored forms in kernel.hpp. For each
q, the weight is computed at the float value of q and at the double value of
q, and each exact weight is rounded once to the nearest double.

The output is the two C++ tables, ready to compare with the test file:

    python3 tools/kernel_reference.py

Standard library only; not part of the build.
"""

import struct
from fractions import Fraction

# Coefficients of q^0, q^1, ... of the expanded polynomials.
C2 = [1, 0, -10, 20, -15, 4]
C4 = [1, 0, Fraction(-28, 3), 0, 70, Fraction(-448, 3), 140, -64, Fraction(35, 3)]


def c2_factored(q):
    return (1 - q) ** 4 * (4 * q + 1)


def c4_factored(q):
    return (1 - q) ** 6 * (35 * q * q + 18 * q + 3) / 3


# The q values of the tables. Non-dyadic values below 1/2 make t = 1 - q
# round, in float and in double, which the test tolerance accounts for.
QS = [0.1, 0.125, 0.25, 0.3, 1 / 3, 0.45, 0.5, 0.6, 0.75, 0.875, 0.9, 0.99, 0.9990234375]

# Relative accuracy is only required where the exact weight is at least
# min() / epsilon() of float (plan, "Relative accuracy").
FLOAT_THRESHOLD = Fraction(2) ** -126 / Fraction(2) ** -23


def to_float(x):
    """static_cast<float>(x): round a double to the nearest float."""
    return struct.unpack("f", struct.pack("f", x))[0]


def weight(coeffs, q):
    """Exact weight. Fraction(q) converts the binary value of q exactly."""
    q = Fraction(q)
    return sum(Fraction(c) * q**i for i, c in enumerate(coeffs))


def table(name, coeffs, factored):
    lines = [f"  constexpr std::array<Reference, {len(QS)}> {name}{{{{"]
    for x in QS:
        w_float = weight(coeffs, to_float(x))
        w_double = weight(coeffs, x)
        # The expanded polynomial is the same function as the factored form.
        assert w_double == factored(Fraction(x))
        assert w_float == factored(Fraction(to_float(x)))
        assert min(w_float, w_double) >= FLOAT_THRESHOLD
        # float() rounds the exact value once; repr() prints the shortest
        # literal that reads back as the same double.
        lines.append(f"      {{{x!r}, {float(w_float)!r}, {float(w_double)!r}}},")
    lines.append("  }};")
    return "\n".join(lines)


if __name__ == "__main__":
    print(table("c2_reference", C2, c2_factored))
    print()
    print(table("c4_reference", C4, c4_factored))
