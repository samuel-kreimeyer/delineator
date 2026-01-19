#!/usr/bin/env python3
"""Run delineator and compare outputs against a reference folder."""

from __future__ import annotations

import argparse
import json
import hashlib
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import rasterio


def _hash_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _round_floats(value): 
    if isinstance(value, float):
        return round(value, 2)
    if isinstance(value, list):
        return [_round_floats(item) for item in value]
    if isinstance(value, dict):
        return {key: _round_floats(val) for key, val in value.items()}
    return value


def _normalize_report(data: dict) -> dict:
    metadata = data.get("metadata")
    if isinstance(metadata, dict):
        metadata = dict(metadata)
        metadata.pop("generated_at", None)
        data = dict(data)
        data["metadata"] = metadata
    return _round_floats(data)


def _compare_json(reference: Path, candidate: Path) -> tuple[bool, str]:
    try:
        ref_data = json.loads(reference.read_text(encoding="utf-8"))
        cand_data = json.loads(candidate.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return False, f"json decode error: {exc}"

    ref_norm = _normalize_report(ref_data)
    cand_norm = _normalize_report(cand_data)

    if ref_norm != cand_norm:
        return False, "json differs (excluding metadata.generated_at)"

    return True, "ok"


def _compare_raster(reference: Path, candidate: Path) -> tuple[bool, str]:
    try:
        with rasterio.open(reference) as ref_src, rasterio.open(candidate) as cand_src:
            if ref_src.profile != cand_src.profile:
                return False, "raster profile differs"
            ref_data = ref_src.read(1)
            cand_data = cand_src.read(1)
    except Exception as exc:
        return False, f"raster read error: {exc}"

    if not np.array_equal(ref_data, cand_data):
        return False, "raster data differs"

    return True, "ok"


def _compare_bytes(reference: Path, candidate: Path) -> tuple[bool, str]:
    if _hash_file(reference) != _hash_file(candidate):
        return False, "file differs"
    return True, "ok"


def _iter_reference_files(reference_dir: Path) -> Iterable[Path]:
    for path in sorted(reference_dir.iterdir()):
        if path.is_file():
            yield path


def run_delineator(args: argparse.Namespace, extra_args: list[str]) -> float:
    command = [
        sys.executable,
        "-m",
        "delineator",
        str(args.input_file),
        "-o",
        str(args.output_dir),
    ]

    if args.use_cogo_points:
        command.append("--use-cogo-points")
    if args.study_points_file:
        command.extend(["--study-points", str(args.study_points_file)])

    command.extend(extra_args)

    start = time.perf_counter()
    subprocess.run(command, check=True)
    return time.perf_counter() - start


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run delineator and compare outputs against a reference folder."
    )
    parser.add_argument("--input", dest="input_file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--use-cogo-points", action="store_true")
    parser.add_argument("--study-points", dest="study_points_file", type=Path)
    parser.add_argument("--skip-run", action="store_true")

    args, extra_args = parser.parse_known_args()

    if not args.skip_run:
        duration = run_delineator(args, extra_args)
        print(f"run: {duration:.2f}s")

    if not args.reference_dir.exists():
        print(f"reference dir not found: {args.reference_dir}")
        return 2
    if not args.output_dir.exists():
        print(f"output dir not found: {args.output_dir}")
        return 2

    failures = 0
    for reference in _iter_reference_files(args.reference_dir):
        candidate = args.output_dir / reference.name
        if not candidate.exists():
            print(f"missing: {candidate}")
            failures += 1
            continue

        ext = reference.suffix.lower()
        if ext == ".tif":
            ok, msg = _compare_raster(reference, candidate)
        elif ext == ".json":
            ok, msg = _compare_json(reference, candidate)
        else:
            ok, msg = _compare_bytes(reference, candidate)

        status = "ok" if ok else "fail"
        print(f"{status}: {reference.name} ({msg})")
        if not ok:
            failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
