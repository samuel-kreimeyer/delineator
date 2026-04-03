# AGENTS.md

## Purpose

CLI interface for the `delineator` watershed delineation package.

## Depends On

- `packages/delineator`

## Local Rules

- Keep argument parsing and exit-code behavior here.
- Do not move core hydrology or parsing logic into this app.
- Preserve the public CLI command name `delineator` unless there is a deliberate breaking-change decision.

## Commands

- Syntax check: `python -m compileall apps/delineator-cli/src`
