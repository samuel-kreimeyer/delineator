# delineator Specification

**Status:** Draft

## Summary
- **Problem:** Delineation of watersheds is a frequent task in hydologic and hydraulic engineering. Determining tributary areas to inlets may need to be frequently recalculated.
- **Outcome:** Given a terrain model and study points, drainage areas are calculated automatically using industry standard methods.
- **Context:** Delineator can be used independently, but will eventually be chained with other programs into engineering workflows. Outputs usable by multiple products is necessary.

## Goals
- Goal 1: Given a terrain model and study points, compute the upstream area of each study point.
- Goal 2: Propose locations for inlets based on plausible capture amounts
- Goal 3: Export results to shape or DXF

## Non-goals
- Non-goal 1: Full hydraulic analysis
- Non-goal 2: High configurability.

## Users & Use Cases
- **Primary user:** Hydraulic engineer designing drainage for transportation
- **Secondary user:** Application developers targeting engineers
- **Core workflow:** The program accepts a raster-based terrain model and inlet points.
- **Key operational scenario:** The roadway engineer provides the drainage engineer with a terrain model of the finished highway. The drainage engineer creates a schematic layout of the drainage system. The drainage engineer inputs a geoTIFF that covers the right of way for 1 mile of roadway (5000x200ft) and 50 pour points. The system computes correct, non-overlapping tributary areas for each inlet.

## User Stories
- US-1: As a hydraulic engineer, I want reliable calculations so that I minimize manual effort
- US-2: As a hydraulic engineer, I want a system that can propose plausible locations for inlets to facilitate schematic design.

## Acceptance Criteria
- Given US-1
- When given a valid elevation model and shapefile
- Then the system returns a valid shapefile or DXF

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
