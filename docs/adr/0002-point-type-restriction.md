# ADR 0002 — require point_type to be Point<scalar_type, dimension>

Status: accepted, with corrections. Supersedes the `point_type` decision in
[ADR 0001](0001-domain-periodic-contracts.md); the rest of ADR 0001 stands.
Two statements below are corrected in the Corrections section at the end.

## Context

ADR 0001 rejected restricting `point_type`, on the grounds that it removed
the concept's only freedom for no present benefit.

Nothing downstream can use a foreign point type. Every consumer of a
`Domain` indexes into points and assumes `std::array` layout. The freedom
is nominal, and leaving it open means the first real consumer either adds
the restriction anyway or works around its absence.

## Decision

`Domain` requires `point_type` to be exactly
`Point<scalar_type, dimension>`, checked after the `Real<scalar_type>` and
`dimension > 0` requirements because `Point<scalar_type, dimension>` is
only well-formed once both hold.

## Consequences

A model spelling `point_type` as `std::array<T, D>` still satisfies the
concept: `Point<T,D>` is an alias for that type. A C array, a pointer, or
a foreign struct is rejected. Covered by
`tests/test_domain_contract.cpp`.

## Corrections

**Ordering rationale.** The Decision says `point_type` is "checked after the
`Real<scalar_type>` and `dimension > 0` requirements because
`Point<scalar_type, dimension>` is only well-formed once both hold." Only
`Real` is needed for that. `Point<T,D>` is an alias template whose scalar
parameter is constrained by `Real`, but `Point<double, 0>` is a valid
`std::array`. Positivity is a separate requirement of `Domain`, not a
condition for forming the point type.

**Necessity.** The Context says "Every consumer of a `Domain` indexes into
points and assumes `std::array` layout." Not every current consumer does:
`distance` only calls `distance_squared`, and `admits_bandwidth` never
touches a point. The restriction is a chosen representation contract,
adopted because the planned consumers such as neighbour search and fitting
will index points. It is not forced by the code that exists today. The
decision stands; only its justification was overstated.
