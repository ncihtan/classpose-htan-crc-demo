#!/usr/bin/env python3
"""Repo checks. Run from the repo root before each commit.

Exits 1 if any check fails. Prints file:line for possible secrets, never the value.
"""
import ast
import os
import re
import subprocess
import sys

MAX_FILE_MB = 50
MAX_REPO_MB = 200
CODE_DIRS = ("pipeline/", "analysis/", "archive/", "patches/", "pannuke/")
PERSONAL = re.compile(r"/Users/ataylor|learn-classpose")
IMAGE_EXT = re.compile(
    r"\.(tif|tiff|svs|qptiff|czi|nd2|ims|vsi|scn|mrxs|zarr|h5|h5ad|hdf5|npy|npz)$", re.I
)
SECRET = re.compile(
    r"AKIA[0-9A-Z]{16}"
    r"|-----BEGIN [A-Z ]*PRIVATE KEY-----"
    r"|gh[pousr]_[A-Za-z0-9]{30,}"
    r"|github_pat_[A-Za-z0-9_]{20,}"
    r"|eyJ[A-Za-z0-9_-]{10,}\.eyJ[A-Za-z0-9_-]{10,}"
    r"|(TOWER_ACCESS_TOKEN|SYNAPSE_AUTH_TOKEN|AWS_SECRET_ACCESS_KEY|accessToken|authToken|apiKey|password)"
    r"[\"' ]*[:=][ \"']*[A-Za-z0-9/+_.=-]{16,}"
    r"|cloud\.seqera\.io/[^ ]*(token|params)"
    r"|\b\d{12}\b(?=.*(arn:|amazonaws|account))"
)

failures = []


def tracked_files():
    out = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        capture_output=True, text=True, check=True,
    ).stdout
    return [f for f in out.split("\0") if f and os.path.isfile(f)]


def text_lines(path):
    try:
        with open(path, "r", encoding="utf-8", errors="strict") as fh:
            return fh.read().splitlines()
    except (UnicodeDecodeError, OSError):
        return None


files = tracked_files()
total = 0
for f in files:
    size = os.path.getsize(f)
    total += size
    if size > MAX_FILE_MB * 1024 * 1024:
        failures.append(f"file over {MAX_FILE_MB} MB: {f}")
    if IMAGE_EXT.search(f):
        failures.append(f"image or array data file: {f}")
    lines = text_lines(f)
    if lines is None:
        continue
    for n, line in enumerate(lines, 1):
        if SECRET.search(line):
            failures.append(f"possible secret: {f}:{n}")
        if f.startswith(CODE_DIRS) and PERSONAL.search(line):
            failures.append(f"personal or SOURCE path in code: {f}:{n}")
    if f.endswith(".py"):
        try:
            ast.parse("\n".join(lines), filename=f)
        except SyntaxError as err:
            failures.append(f"Python syntax error: {f}:{err.lineno}")

if total > MAX_REPO_MB * 1024 * 1024:
    failures.append(f"repo files total {total / 1048576:.0f} MB (limit {MAX_REPO_MB} MB)")

for f in failures:
    print("FAIL", f)
print(f"{len(files)} files, {total / 1048576:.1f} MB, {len(failures)} failures")
sys.exit(1 if failures else 0)
