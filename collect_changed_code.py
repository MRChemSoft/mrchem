#!/usr/bin/env python3

"""
collect_changed_code.py

Purpose
-------
Collect all lines from files you've **changed** in a Git repo and write them
into a single text file with line numbers, grouped by file.

By default, it gathers:
- Unstaged changes (working tree vs HEAD)
- Staged changes (index vs HEAD)
- Untracked files (not ignored)

You can also compare against a specific ref/branch with --ref (e.g., main).

Usage
-----
# From the repo root:
python3 collect_changed_code.py

# Compare against a branch/ref:
python3 collect_changed_code.py --ref main

# Save to a custom outfile:
python3 collect_changed_code.py --outfile changed_code.txt

# Only include specific extensions:
python3 collect_changed_code.py --ext .py,.js,.ts

# Remove the default 5 MB per-file size guard:
python3 collect_changed_code.py --no-size-limit
"""
import argparse
import os
import subprocess
import sys
import time
from typing import List, Set, Tuple


def run(cmd: List[str]) -> str:
    """Run a command and return stdout (text). Raise on non-zero exit."""
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Command {' '.join(cmd)} failed:\n{r.stderr.strip()}")
    return r.stdout


def get_repo_root() -> str:
    """Return absolute path to the repository root."""
    root = run(["git", "rev-parse", "--show-toplevel"]).strip()
    return root


def list_changed_files(ref: str = None) -> List[str]:
    """
    Return a sorted list of changed/untracked file paths (relative to repo root).

    If ref is provided, returns files changed between ref...HEAD.
    Otherwise, combines:
      - Unstaged changes (git diff --name-only)
      - Staged changes (git diff --name-only --cached)
      - Untracked files (git ls-files --others --exclude-standard)
    """
    files: Set[str] = set()
    if ref:
        out = run(["git", "diff", "--name-only", f"{ref}...HEAD"])
        files.update([line.strip() for line in out.splitlines() if line.strip()])
    else:
        for cmd in [
            ["git", "diff", "--name-only"],
            ["git", "diff", "--name-only", "--cached"],
            ["git", "ls-files", "--others", "--exclude-standard"],
        ]:
            out = run(cmd)
            files.update([line.strip() for line in out.splitlines() if line.strip()])

    # Filter out paths that no longer exist (deleted, etc.)
    existing = [p for p in files if os.path.exists(p)]
    existing.sort()
    return existing


def looks_text(sample: bytes) -> bool:
    """Heuristic to check if a byte sample is text-like (not binary)."""
    if b"\x00" in sample:
        return False
    if not sample:
        return True
    # Allow typical printable ASCII and common whitespace/newlines
    text_like = sum((32 <= b <= 126) or b in (9, 10, 13) for b in sample)
    ratio = text_like / len(sample)
    return ratio >= 0.70


def is_binary_file(path: str, sample_size: int = 8192) -> bool:
    try:
        with open(path, "rb") as f:
            sample = f.read(sample_size)
        return not looks_text(sample)
    except Exception:
        # If we can't read or something goes wrong, assume binary to be safe.
        return True


def should_include(path: str, exts: List[str]) -> bool:
    if not os.path.isfile(path):
        return False
    if exts:
        plower = path.lower()
        return any(plower.endswith(e.lower().strip()) for e in exts if e.strip())
    return True


def human_size(bytes_num: int) -> str:
    for unit in ["B", "KB", "MB", "GB"]:
        if bytes_num < 1024.0:
            return f"{bytes_num:.1f} {unit}"
        bytes_num /= 1024.0
    return f"{bytes_num:.1f} TB"


def collect_lines(
    files: List[str],
    outfile: str,
    exts: List[str],
    no_size_limit: bool,
    max_mb: float = 5.0,
) -> Tuple[int, int, int]:
    """
    Write all lines from target files into outfile.
    Returns (num_files_included, total_lines_written, num_files_skipped).
    """
    repo_root = get_repo_root()
    abs_out = os.path.abspath(outfile)

    included = 0
    skipped = 0
    total_lines = 0
    max_bytes = int(max_mb * 1024 * 1024)

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(outfile, "w", encoding="utf-8", errors="replace") as out:
        out.write(f"# Changed Code Dump\n")
        out.write(f"# Repo: {repo_root}\n")
        out.write(f"# Generated: {timestamp}\n\n")

        # Index
        out.write("## Files Included\n")

        index_entries = []
        file_infos = []  # (path, is_binary, size_bytes, line_count)
        for rel in files:
            # Respect extension filter
            if not should_include(rel, exts):
                skipped += 1
                continue

            size_bytes = 0
            try:
                size_bytes = os.path.getsize(rel)
            except OSError:
                skipped += 1
                continue

            if not no_size_limit and size_bytes > max_bytes:
                file_infos.append((rel, False, size_bytes, -1))
                skipped += 1
                continue

            binary = is_binary_file(rel)
            if binary:
                file_infos.append((rel, True, size_bytes, 0))
                skipped += 1
                continue

            # Count lines to put in index (we'll read again later)
            try:
                with open(rel, "r", encoding="utf-8", errors="replace") as f:
                    lines = f.read().splitlines()
                file_infos.append((rel, False, size_bytes, len(lines)))
            except Exception:
                file_infos.append((rel, False, size_bytes, -1))
                skipped += 1

        for rel, binary, size_bytes, line_count in file_infos:
            if binary or line_count <= 0:
                continue
            index_entries.append(f"- {rel} ({line_count} lines, {human_size(size_bytes)})")

        if index_entries:
            out.write("\n".join(index_entries) + "\n\n")
        else:
            out.write("(no text files to include)\n\n")

        # Body
        for rel, binary, size_bytes, line_count in file_infos:
            if binary or line_count <= 0:
                continue

            included += 1
            out.write(f"=== FILE: {rel} | {human_size(size_bytes)} ===\n")
            try:
                with open(rel, "r", encoding="utf-8", errors="replace") as f:
                    for i, line in enumerate(f.read().splitlines(), 1):
                        out.write(f"{i:6d} | {line}\n")
                        total_lines += 1
            except Exception as e:
                out.write(f"[ERROR] Could not read file: {e}\n")
            out.write("\n")

        out.write(f"# Summary: {included} files, {total_lines} lines written. Skipped {skipped}.\n")

    return included, total_lines, skipped


def main():
    parser = argparse.ArgumentParser(description="Collect lines from changed files into a single text file.")
    parser.add_argument("--ref", help="Compare against this ref/branch instead of using working-tree status (e.g., 'main').")
    parser.add_argument("--outfile", default="changed_code.txt", help="Path to output text file (default: changed_code.txt).")
    parser.add_argument(
        "--ext",
        help="Comma-separated list of file extensions to include (e.g., .py,.js,.ts). Default: include all.",
    )
    parser.add_argument("--no-size-limit", action="store_true", help="Disable the 5 MB per-file size limit.")
    args = parser.parse_args()

    # Ensure we're inside a git repo
    try:
        repo_root = get_repo_root()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    os.chdir(repo_root)

    try:
        files = list_changed_files(args.ref)
    except Exception as e:
        print(f"Error listing changed files: {e}", file=sys.stderr)
        sys.exit(2)

    exts = None
    if args.ext:
        exts = [e.strip() for e in args.ext.split(",") if e.strip()]

    included, total_lines, skipped = collect_lines(
        files=files,
        outfile=args.outfile,
        exts=exts,
        no_size_limit=args.no_size_limit,
    )

    print(f"Wrote {total_lines} lines from {included} files to {os.path.abspath(args.outfile)} (skipped {skipped}).")


if __name__ == "__main__":
    main()

