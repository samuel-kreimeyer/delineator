# AGENTS.md

## Purpose

Core watershed delineation logic for LandXML-driven drainage workflows.

## Invariants

- Keep reusable hydrology, parsing, geometry, and reporting logic in this package.
- Do not reintroduce CLI-specific argument parsing here.
- Preserve coordinate-system handling and output-unit behavior explicitly.

## Local Rules

- The app in `apps/delineator-cli` is the command-line interface for this package.
- Keep algorithmic changes covered by fixtures or regression-oriented tests as they are added.
- Preserve source documentation and reference material unless it is clearly obsolete.

## Commands

- Syntax check: `python -m compileall packages/delineator/src`
