#!/usr/bin/env python3
"""Merge tabular text/gzip files while retaining only the first header."""

import argparse
import gzip


def open_text(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


parser = argparse.ArgumentParser()
parser.add_argument("--output", required=True)
parser.add_argument("inputs", nargs="+")
args = parser.parse_args()

header = None
with open_text(args.output, "wt") as out:
    for path in args.inputs:
        with open_text(path, "rt") as src:
            current = src.readline()
            if not current:
                continue
            if header is None:
                header = current
                out.write(current)
            elif current.rstrip("\n") != header.rstrip("\n"):
                raise SystemExit(f"Header mismatch: {path}")
            for line in src:
                out.write(line)

