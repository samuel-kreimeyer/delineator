# delineator-cli Specification

**Status:** Draft

## Summary
- **Problem:** <!-- What pain point or gap does this solve? -->
- **Outcome:** <!-- What is the user able to do when this is done? -->
- **Context:** <!-- Where does this live in the broader system? What depends on it or does it depend on? -->

## Goals
- Goal 1: <!-- Primary capability to deliver -->
- Goal 2:
- Goal 3:

## Non-goals
- Non-goal 1: <!-- Explicitly out of scope for this version -->
- Non-goal 2:

## Users & Use Cases
- **Primary user:** <!-- Who is the main person using this? -->
- **Secondary user:** <!-- Other relevant personas -->
- **Core workflow:** <!-- The main thing the user does with this -->
- **Key operational scenario:** <!-- Walk through a realistic end-to-end example -->

## User Stories
- US-1: As a [role], I want [capability] so that [benefit].
- US-2:
- US-3:

## Acceptance Criteria
- Given [context] and [user story]
- When [action or event]
- Then [observable outcome]

## Architecture
- **Core components:** <!-- List the major logical pieces -->
- **Boundaries and responsibilities:** <!-- What does this own vs. what do other packages/apps own? -->
- **Constraints:** <!-- Language, runtime, deployment, dependencies, monorepo rules -->
- **Deployment / runtime assumptions:** <!-- Where does this run? What can it assume about its environment? -->

## Data Model
<!-- Describe the main entities, their fields, and relationships. Use a table or bullet list.
     Note whether data is in-memory only, persisted to disk, or stored in a database. -->

| Entity | Fields | Notes |
| --- | --- | --- |
| `Foo` | `id`, `name`, ... | ... |

## Interfaces
- **User-facing interfaces:** <!-- CLI flags, web UI pages, API endpoints -->
- **Programmatic interfaces:** <!-- Public functions / classes other packages call -->
- **Inputs:** <!-- File formats, environment variables, config schema -->
- **Outputs:** <!-- Files produced, return values, side effects -->
- **Configuration:** <!-- How is the tool configured? Config file, flags, env vars? -->

## Error Handling
- <!-- How are invalid inputs reported? -->
- <!-- How are runtime failures surfaced? Exit codes, exceptions, log messages? -->
- <!-- What partial-failure behavior is acceptable? -->

## Testing Plan
- <!-- Unit tests: which modules / functions need coverage? -->
- <!-- Integration tests: which cross-component paths need an end-to-end test? -->
- <!-- Golden-file or snapshot tests? -->
- <!-- Coverage target -->

## Risks & Open Questions
- <!-- Technical unknowns that could affect the design -->
- <!-- Dependencies on external data, tools, or teams -->
- <!-- Decisions that need to be made before implementation can begin -->
