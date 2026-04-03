# Delineator

Core watershed delineation library migrated from the standalone `delineator` repository.

This package owns:

- LandXML parsing
- TIN rasterization
- hydrology preprocessing and flow modeling
- watershed delineation
- boundary and report generation

The CLI interface lives in [`/home/sam/Projects/Curatores Viarum/apps/delineator-cli`](/home/sam/Projects/Curatores%20Viarum/apps/delineator-cli).

Source-era documentation was preserved in:

- [`README.source.md`](/home/sam/Projects/Curatores%20Viarum/packages/delineator/README.source.md)
- [`docs/`](/home/sam/Projects/Curatores%20Viarum/packages/delineator/docs)
- [`perf/`](/home/sam/Projects/Curatores%20Viarum/packages/delineator/perf)

## Migration Notes

- The original package-level CLI entrypoint was moved out of the library into the monorepo app.
- No formal external schema was present in the source repository, so no schema extraction was added yet.
