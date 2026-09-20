# ADR 0002 — require point_type to be Point<scalar_type, dimension>

Status: accepted. Supersedes the `point_type` decision in
[ADR 0001](0001-domain-periodic-contracts.md); the rest of ADR 0001 stands.

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
