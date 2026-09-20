# Developing tedder with Codex and Claude Code

Status: proposed working method. This document defines how AI should assist the
project; it does not replace the technical decisions in `PROJECT.md` or accepted
ADRs.

## Recommendation

Keep the setup deliberately small:

1. Use one short `AGENTS.md` as the shared, automatically loaded instruction file.
2. Make the human the default author of learning-heavy C++ code. Use AI as a
   tutor and reviewer unless a task is explicitly delegated.
3. Keep design and roadmap documents in the repository, but load them only when
   relevant.
4. Add one local `tedder-docs` skill backed by a bounded `ripgrep` script.
5. Let only one participant write code for a task. For core library work that
   participant is normally you; an AI can independently review the resulting diff.
6. Use subagents only for bounded, read-heavy work such as documentation lookup,
   code mapping, test-log analysis, and independent reviews.
7. Do not build a documentation MCP server or semantic-search system yet.

This is enough for the current repository. More infrastructure would cost more
context and maintenance than it saves.

## Learning contract

The purpose of tedder is not merely to produce a library. It is to learn modern
C++ while building one. Optimize the workflow for understanding and deliberate
practice, not implementation speed.

Use three modes of AI assistance:

| Mode | Use it for | AI may write production code? |
| --- | --- | --- |
| Tutor | Concepts, generic design, ranges, constraints, ownership, API tradeoffs, and numerical reasoning | No. Explain alternatives, ask questions, and use only small isolated examples. |
| Reviewer | Code you wrote, mathematical assumptions, diagnostics, portability, and missing tests | No. Report findings and possible directions; do not silently replace the implementation. |
| Delegate | Mechanical work with an accepted contract: tests, CI/CMake plumbing, formatting, documentation, repetitive adapters, and `.npy` headers or I/O scaffolding | Yes, but only when explicitly requested. |

The boundary is based on learning value, not file type. You should normally write:

- public concepts and APIs, including the future kernel concept;
- project-specific generic C++ and template machinery;
- neighbour-search and fitting abstractions;
- the local-polynomial implementation and its numerical policies.

It is reasonable to depend on established implementations of generic algorithms
that are not the subject of the exercise. In particular, tedder need not
reimplement SVD or a kd-tree. Eigen is the likely linear-algebra dependency and
nanoflann the accepted neighbour-search dependency; integrating them still
requires understanding and testing their contracts.

For an undecided design such as `Kernel`, do not ask an agent to invent the final
concept. Ask it to identify the smallest operations the current algorithm needs,
show two or three candidate interfaces, and explain what each would permit or
exclude. Then write a concrete kernel and generalize only after a second model or
call site creates real pressure.

## Modern development environment

Keep local development reproducible without hiding the compiler from you.

### Repository tooling

Add these in two small, reviewable changes:

1. **Compiler workflow:** a checked-in `CMakePresets.json` with `dev-clang`,
   `dev-gcc`, and `asan-ubsan` presets; `CMAKE_EXPORT_COMPILE_COMMANDS=ON`; and
   `CMakeUserPresets.json` ignored for machine-local overrides.
2. **Advisory analysis:** a conservative `.clang-tidy`, a formatting check, and
   warnings owned by tedder's test targets rather than exported to consumers.

Suggested daily loop:

```sh
cmake --preset dev-clang
cmake --build --preset dev-clang
ctest --preset dev-clang
```

Use the Clang debug compilation database for clangd. Keep GCC as an independent
compiler check rather than trying to make the editor switch databases constantly.
Use `ccache` when available, but do not make it a requirement. Run ASan and UBSan
together in a separate preset. Add TSan only when tedder contains concurrent code.

Keep `clang-tidy` advisory: do not enable automatic fixes over learning-heavy code.
Start with `bugprone-*`, `performance-*`, `portability-*`, and a small explicitly
chosen subset of `modernize-*`. Avoid enormous style profiles and do not use
`-Werror` in the normal developer preset.

### Editor and debugging

For Neovim, the essential integration is clangd reading the generated
`compile_commands.json`. Formatting on demand or save is useful; automatic
`clang-tidy --fix` is not. Keep GDB available and add LLDB only if it provides a
specific benefit for a Clang diagnostic or sanitizer failure.

### CI

Retain GCC and Clang Debug/Release builds. Add one Clang ASan+UBSan job, then a
format/tidy job only after its local commands are stable. CI should verify the
same presets where practical; it should not introduce an unrelated second build
system.

### AI tooling

Start with repository instructions, the small documentation skill, and ordinary
local tools. No project-specific MCP server is needed. The GitHub integration is
useful for remote CI and review state; source editing and tests should remain
local. Add a security plugin later, when `.npy` parsing or other untrusted input
enters the library.

## Repository knowledge layout

```text
AGENTS.md                         always-loaded project rules, kept short
README.md                         public description of current capabilities
PROJECT.md                        stable technical direction and mathematical conventions
ROADMAP.md                        stage status and links, not detailed design
docs/
  INDEX.md                        compact routing table for humans and agents
  adr/
    0001-periodic-validity.md     accepted decisions and their consequences
    0002-solver-policy.md
  stages/
    00-foundation.md              active or completed stage contract
    01-kernels.md
    02-neighbour-search.md
    ...
  archive/                        superseded material, excluded from normal search
tools/
  doc-search                      bounded local documentation search
.ai/skills/
  tedder-docs/SKILL.md            canonical shared skill
.agents/skills/tedder-docs        symlink to the canonical skill for Codex
.claude/skills/tedder-docs        symlink to the canonical skill for Claude Code
```

Do not store the same rule in `AGENTS.md`, `PROJECT.md`, and a stage file.

- `AGENTS.md` says **how to work**.
- `PROJECT.md` says **what tedder is intended to become**.
- ADRs say **what was decided and why**.
- A stage file says **what is being implemented now and how completion is
  judged**.
- `ROADMAP.md` says **where the project currently is**.

The existing `PLAN.md` and `FIELD_PLAN.md` belong to stage 00. After that work is
finished, preserve their lasting decisions as ADRs and mark the stage complete.
Do not keep generic `PLAN.md` files accumulating at the repository root.

## Shared agent instructions

Use `AGENTS.md` as the single shared project instruction file. Current Codex reads
it directly. Current Claude Code also reads `AGENTS.md` when there is no competing
project `CLAUDE.md`; this requires Claude Code 2.1.277 or later. If a Claude setup
cannot read it directly, create a minimal `CLAUDE.md` containing `@AGENTS.md`
rather than maintaining a second copy.

Keep `AGENTS.md` below roughly 100 lines. It should contain only durable rules,
for example:

```markdown
# tedder agent instructions

- This is an experimental header-only C++23 numerical library.
- `PROJECT.md` describes design direction; accepted ADRs override proposals.
- Before non-trivial work, read `docs/INDEX.md`, then run
  `tools/doc-search '<terms>'`. Do not read all project documentation.
- Work on one roadmap stage at a time. Read its stage file before editing.
- Do not turn undecided policies into APIs. Record a decision before implementing it.
- Preserve public numerical conventions, including offset sign and Jacobian layout.
- The project is a modern-C++ learning exercise. For core library work, default
  to tutor or reviewer mode and do not write production code unless explicitly asked.
- Fully delegated work is normally mechanical: tests for an accepted contract,
  build/CI plumbing, documentation, and file-format boilerplate.
- One participant owns code changes for a task. Review agents are read-only.
- Keep unrelated user changes intact.
- Run the smallest relevant tests while iterating and the full configured suite
  before handoff. Report commands actually run.
- Never claim floating-point identities are exact without a proof for the stated
  representation and operation sequence.
```

Do not import `PROJECT.md`, `ROADMAP.md`, or all ADRs from `AGENTS.md`. Imports
would place them in every session and defeat the retrieval scheme.

Official references: [Codex `AGENTS.md`](https://learn.chatgpt.com/docs/agent-configuration/agents-md),
[Claude Code memory and `AGENTS.md`](https://code.claude.com/docs/en/memory).

## Cheap documentation search

### First level: a small routing index

`docs/INDEX.md` should be a compact table, not a summary of every document:

```markdown
| Topic | Authoritative file | Status | Search terms |
| --- | --- | --- | --- |
| Overall design | `PROJECT.md` | active | local polynomial, degree, solver |
| Current work | `docs/stages/00-foundation.md` | active | Domain, Periodic, tests |
| Periodic validity | `docs/adr/0001-periodic-validity.md` | accepted | NaN, infinity, bandwidth |
| Solver policy | `docs/adr/0002-solver-policy.md` | proposed | LDLT, QR, SVD, Eigen |
```

The index lets an agent select likely files without opening every document.

### Second level: bounded lexical search

Start with `ripgrep`. The documentation set is small, technical terms are precise,
and lexical search is deterministic and nearly free. A suitable `tools/doc-search`
script is:

```bash
#!/usr/bin/env bash
set -euo pipefail

if (($# != 1)); then
  echo "usage: tools/doc-search 'regular expression'" >&2
  exit 2
fi

roots=(docs)
for file in AGENTS.md PROJECT.md ROADMAP.md; do
  [[ -f "$file" ]] && roots+=("$file")
done

set +e
rg -n -i \
  --glob '*.md' \
  --glob '!docs/archive/**' \
  --max-columns 240 \
  --max-columns-preview \
  -- "$1" "${roots[@]}" | awk 'NR <= 40 { print }'
status=${PIPESTATUS[0]}
set -e

((status == 0 || status == 1)) || exit "$status"
```

The agent workflow is:

1. Read `docs/INDEX.md`.
2. Run one focused query, for example
   `tools/doc-search 'periodic|minimum image|bandwidth'`.
3. Inspect headings in the selected file with `rg -n '^#{1,3} ' FILE`.
4. Read only the relevant section.
5. Cite the file and heading when using a decision.

The 40-line output bound matters more than sophisticated ranking: it prevents a
broad query from injecting a large amount of text into the model context.

### Third level: local full-text search, only if needed

Add SQLite FTS5 only when the repository has roughly 100 substantial Markdown
files or repeated searches fail because terminology has drifted. Index Markdown
by `(path, heading, body)` and return at most eight ranked headings with short
snippets. Rebuild the index locally; do not commit it.

Do not add embeddings, a vector database, or a documentation MCP server before
that threshold. They add dependencies, nondeterminism, tool descriptions, and
another stale index without solving a current problem.

## One shared documentation skill

Create a small `tedder-docs` skill whose full body is loaded only when project
documentation is needed:

```markdown
---
name: tedder-docs
description: Find tedder design decisions, stage contracts, and numerical conventions. Use before planning or changing public behaviour; do not use for ordinary source-code lookup.
---

1. Read `docs/INDEX.md`.
2. Run `tools/doc-search` with focused terms from the task.
3. Open at most three likely documents, and only the relevant sections.
4. Prefer accepted ADRs over stage proposals and historical notes.
5. Return the relevant file, heading, status, and concise conclusion.
6. If sources conflict, stop and report the conflict instead of choosing silently.
```

Use a canonical skill directory under `.ai/skills/` and symlink it into the two
tool-specific project skill directories. Both tools support symlinked skill
folders. If portability to systems without symlink support becomes important,
keep two tiny wrappers that both call the same `tools/doc-search` script.

Skills are appropriate here because only their short descriptions occupy startup
context; the full instructions load on invocation. Keep the skill small because
invoked skill content remains in the conversation. See
[Codex skills](https://learn.chatgpt.com/docs/build-skills) and
[Claude Code skills](https://code.claude.com/docs/en/skills).

Do not create more project skills until a workflow has been repeated at least
three times. A future `tedder-stage` skill may be useful, but first learn which
steps are genuinely stable.

## How to divide work between you, Codex, and Claude Code

Use the two products as tutors and independent checks, not as two simultaneous
authors. You are normally the author of core C++.

| Work | Recommended owner | Other system |
| --- | --- | --- |
| Learn a language feature or design a concept | You, with one AI tutoring | Other system only if a genuinely different view is useful |
| Clarify mathematics and API contract | Main conversation with you | Challenge assumptions read-only |
| Search external dependency documentation | Cheap research subagent | Verify only disputed claims |
| Implement core library code | You | One system reviews read-only |
| Implement explicitly delegated mechanical work | One Codex or Claude Code session | No edits |
| Review your implementation | One AI system | The second is unnecessary unless risk is high or findings conflict |
| Resolve review findings | You, unless the fix is explicitly delegated | Re-review only changed areas |
| Approve policy choices | You | Agents present evidence and alternatives |

Do not rotate core implementation merely to use both products. A reasonable
sequence is: discuss a narrow problem with one system, write the code yourself,
then ask the other system to review the diff without editing it. Delegate a full
implementation only when the work has low learning value and a settled contract.

Avoid asking both systems to “improve the library.” Give each task a stage ID,
files in scope, acceptance tests, and explicit exclusions.

## Subagents and custom agents

Use the main conversation for teaching and design because those phases benefit
from dialogue. Use subagents for work that is independent and produces noisy
intermediate output:

- locate relevant documentation and return at most five references;
- map affected code without editing;
- run a test matrix and summarize failures;
- review numerical correctness, API compatibility, or test gaps independently.

Do not spawn write-capable agents against learning-heavy code such as
`field.hpp`, kernels, or fitting machinery unless you explicitly delegate it.
Parallel writers save little time in a small library and create coordination
errors. Subagents also consume additional tokens, so one bounded research agent
is usually enough.

Initially, use the built-in explorer/reviewer behaviour with explicit prompts.
Create a permanent `numerical-reviewer` only after the same review prompt has been
useful several times. It should be read-only and focus on:

- floating-point range, cancellation, and tie conventions;
- rank and conditioning assumptions;
- coordinate and bandwidth scaling of derivatives;
- scalar/vector dimensional constraints;
- tests whose oracle merely repeats the implementation.

Codex project agents live under `.codex/agents/`; Claude Code project agents live
under `.claude/agents/`. Their configuration formats differ, so do not pretend one
agent definition is portable. Keep the short behavioural checklist above in an
ADR or reference file and make both definitions point reviewers to it.

Official references: [Codex subagents](https://learn.chatgpt.com/docs/agent-configuration/subagents),
[Claude Code subagents](https://code.claude.com/docs/en/sub-agents).

## MCP servers and plugins

Use fewer integrations than either tool permits.

| Integration | Recommendation |
| --- | --- |
| GitHub plugin | Keep it. It is already installed and is useful for remote PR, issue, and CI state. For local code and diffs, prefer `git`, `gh`, and repository tools. |
| Local documentation MCP | Do not build it. `docs/INDEX.md` plus bounded `rg` is cheaper and easier to audit. |
| General documentation MCP | Optional later for version-specific Eigen, nanoflann, Catch2, or CMake research. Give it only to a read-only research agent and verify important claims against primary documentation. |
| Codex Security plugin | Defer. Reconsider before a release or when parsers, file I/O, untrusted inputs, or a larger dependency surface appear. |
| Linear, Notion, Google Drive | Do not use for the source of truth. This is a one-person open-source repository; Git and GitHub Issues are enough. |
| Broad workflow plugins | Avoid initially. They can inject a methodology that conflicts with the project’s deliberately small stages. |

Add an MCP server only when it accesses an external source of truth that local
commands cannot access and the workflow recurs. Restrict its enabled tools and
approvals. MCP is an integration boundary, not a replacement for repository
documentation. Both Codex and Claude Code can discover MCP tools on demand, but
every server still adds configuration and another failure mode. See
[Codex MCP](https://learn.chatgpt.com/docs/extend/mcp) and
[Claude Code MCP](https://code.claude.com/docs/en/mcp).

## Development stages

`ROADMAP.md` should contain only this table plus links and current status. Each
active stage gets one detailed file under `docs/stages/`.

| Stage | Deliverable | Completion gate | AI use |
| --- | --- | --- | --- |
| 00 — Foundation repair | Implement accepted `field.hpp` and testing plans; add the documentation/search skeleton | GCC/Clang Debug and Release pass; sanitizer job passes; accepted regressions pass | You decide and write core API changes; AI may implement accepted tests and tooling, then review |
| 01 — Kernels | Compactly supported kernels, starting with the agreed Wendland subset | Values at zero/support boundary, non-negativity, symmetry, and promised smoothness tested | AI explains formulas and concept alternatives; you design and implement; AI may generate contract tests afterward |
| 02 — Reference neighbour search | Search concept/API plus brute-force implementation with fixed bandwidth | Exact comparison against a simple hand-built oracle; empty, boundary, mixed-domain, and periodic-seam cases | You design and implement; a research agent checks terminology and edge cases; reviewer checks units and inclusivity |
| 03 — Degree 0 and 1 fitting | Fixed-bandwidth local fits and explicit `FitError`; linear algebra supplied by the chosen dependency | Constant and affine reproduction within tolerance; scalar and vector fields; insufficient/rank-deficient neighbourhood failures; physical Jacobian scaling by `1/h` | Collaborative derivation; you implement project logic; AI builds independent oracles/tests and reviews conditioning claims |
| 04 — Degree 2 fitting | Quadratic basis and curvature-aware fit while retaining the existing primary result contract | Quadratic reproduction within tolerance, including mixed terms; Jacobian scaling; degenerate layouts fail explicitly | AI checks an independent basis/scaling derivation; you implement; AI may extend the test matrix |
| 05 — nanoflann backend | kd-tree implementation behind the same neighbour-search contract | Results agree with brute force on supported domains, tolerances, duplicates, and boundary cases | Research agent summarizes nanoflann's API; you implement the adapter; AI writes backend-equivalence tests if delegated |
| 06 — Derived quantities | Divergence, dimension-appropriate vorticity, and strain from value/Jacobian data | Analytic matrices/fields cover scalar/vector and dimensional constraints | You implement the small formulas; AI independently reviews them and may write repetitive dimension tests |
| 07 — Optional I/O and policies | NumPy I/O and only explicitly accepted advanced policies | Separate ADR and tests for each feature; no accidental default-policy changes | `.npy` boilerplate may be fully delegated after format and API decisions are accepted |

Stage 03 should choose and record the linear-algebra dependency before exposing a
solver-backed public API. Eigen is the current likely choice, not an irreversible
decision. Use pivoted QR as the robust comparison or fallback; an `LDLT` normal-
equations path must earn its place through conditioning tests and benchmarks.

Stage 05 has an unresolved design point: a conventional kd-tree does not
automatically reproduce a per-axis periodic metric. Decide and document supported
periodic behaviour before promising backend equivalence. Do not hide this in the
adapter implementation.

## Stage document template

```markdown
# Stage NN — Name

Status: proposed | accepted | in progress | blocked | complete
Owner: human author, or explicitly delegated AI implementer

## Goal
One paragraph.

## Inputs
Accepted ADRs and existing APIs this stage relies on.

## Contract
Observable behaviour, numerical conventions, and error behaviour.

## Acceptance tests
Concrete tests and supported toolchains.

## Out of scope
Ideas that must not enter this stage accidentally.

## Open decisions
Questions requiring evidence and explicit user approval.

## Handoff
Files changed, commands actually run, known limitations, and review status.
```

An AI may draft this file, but you accept the contract and open decisions before
implementation begins.

## Prompt patterns

Planning:

```text
Work on stage 03 only. Read AGENTS.md, docs/INDEX.md, and the stage file; use
tedder-docs for anything else. Do not edit code yet. Identify unresolved API or
numerical decisions, propose at most two alternatives, and state which tests
would distinguish them.
```

Learning implementation:

```text
Act as a tutor for stage 03. Do not write the production implementation. Explain
the next design decision, give at most three alternatives with consequences, and
ask me to choose or sketch an interface. Use tiny isolated examples only when
they clarify a C++ feature. After I write code, review it without editing.
```

Delegated implementation:

```text
Implement only the accepted mechanical task described in the stage file. This
task is explicitly delegated. Do not redesign the public API or numerical policy.
Keep unrelated files unchanged, run the relevant tests, and report assumptions
that the contract did not settle.
```

Independent review:

```text
Review this branch against main and the accepted stage contract. Do not edit.
Prioritize concrete correctness defects, numerical assumptions, API regressions,
and missing tests. Reproduce findings where practical. Return findings first,
with file references; omit style comments unless they hide a defect.
```

Documentation lookup:

```text
Use tedder-docs to find the authoritative decision about periodic tie handling.
Return only the file, heading, status, and conclusion. Report conflicts.
```

## Maintenance rules

- Update `docs/INDEX.md` in the same change that creates, supersedes, or renames a
  design document.
- Never silently edit an accepted ADR. Add a superseding ADR.
- At stage completion, keep the contract and handoff; discard raw AI transcripts
  and repetitive review chatter.
- Keep generated benchmark output outside automatically searched documentation.
- Review `AGENTS.md` after repeated agent mistakes, not after every task.
- Review skills after repeated workflows, not because a capability exists.
- Reconsider FTS5, additional agents, and MCP servers only when measured friction
  appears.

The result should feel boring: small files, explicit decisions, bounded searches,
one writer, one independent reviewer. That is a better foundation for numerical
C++ than a large autonomous-agent stack.
