#!/usr/bin/env python3
"""
Dump selected files' full contents into a single text (and optional Markdown) file.

Usage:
  python dump_mrchem_files.py [--markdown]

- Always writes: mrchem_selection_dump_YYYYMMDD_HHMMSS.txt
- If --markdown is given, also writes: mrchem_selection_dump_YYYYMMDD_HHMMSS.md
"""

from __future__ import annotations
from datetime import datetime
from pathlib import Path
import argparse
import sys

def infer_lang(path: str) -> str:
    """
    Infer a Markdown code fence language from the file extension.
    Defaults to empty string (no tag) if not recognized.
    """
    ext = Path(path).suffix.lower()
    return {
        ".h": "cpp",
        ".hpp": "cpp",
        ".hh": "cpp",
        ".c": "c",
        ".cc": "cpp",
        ".cpp": "cpp",
        ".cxx": "cpp",
        ".cu": "cuda",
        ".py": "python",
        ".cmake": "cmake",
        ".in": "",  # template; leave untagged
        ".json": "json",
        ".txt": "",
        ".md": "markdown",
        ".toml": "toml",
        ".yml": "yaml",
        ".yaml": "yaml",
    }.get(ext, "")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--markdown", action="store_true", help="Also write a Markdown version with code fences.")
    args = parser.parse_args()

    # Exact list (in order) as requested
    files = [
        # Core DFT backend selection & wiring
        "src/mrdft/Factory.h",
        "src/mrdft/Factory.cpp",
        "src/mrdft/LibXCBackend.h",
        "src/mrdft/LibXCBackend.cpp",
        "src/mrdft/xc_backend.h",
        "src/mrdft/xc_backend_xcfun.cpp",

        # Program environment & JSON intake
        "src/mrenv.h",
        "src/mrenv.cpp",
        "src/utils/json_utils.h",
        "src/utils/json_utils.cpp",  # optional; included if it exists

        # Build system
        "CMakeLists.txt",
        "src/mrdft/CMakeLists.txt",
        "external/upstream/fetch_libxc.cmake",
        # Other CMake that may gate/define MRCHEM_ENABLE_LIBXC
        "cmake/custom/main.cmake",

        # Optional but helpful (Python front-end)
        "python/mrchem/cli.py",
        "python/mrchem/config.py.in",
    ]

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_txt = f"mrchem_selection_dump_{ts}.txt"
    out_md = f"mrchem_selection_dump_{ts}.md"

    missing_any = False

    # Write TXT bundle
    with open(out_txt, "w", encoding="utf-8", errors="replace") as o:
        for f in files:
            p = Path(f)
            if p.is_file():
                o.write(f"===== BEGIN {f} =====\n")
                try:
                    o.write(p.read_text(encoding="utf-8"))
                except UnicodeDecodeError:
                    o.write(p.read_bytes().decode("utf-8", errors="replace"))
                o.write(f"\n===== END {f} =====\n\n")
            else:
                o.write(f"===== MISSING {f} =====\n\n")
                missing_any = True

    # Optional Markdown bundle with code fences
    if args.markdown:
        with open(out_md, "w", encoding="utf-8", errors="replace") as m:
            m.write(f"# MRChem selected files (dumped {ts})\n\n")
            for f in files:
                p = Path(f)
                m.write(f"## {f}\n\n")
                if p.is_file():
                    lang = infer_lang(f)
                    fence = lang if lang else ""
                    m.write(f"```{fence}\n")
                    try:
                        m.write(p.read_text(encoding="utf-8"))
                    except UnicodeDecodeError:
                        m.write(p.read_bytes().decode("utf-8", errors="replace"))
                    m.write("\n```\n\n")
                else:
                    m.write("_Missing file_\n\n")
                    missing_any = True

    print(f"Wrote {out_txt}")
    if args.markdown:
        print(f"Wrote {out_md}")
    if missing_any:
        print("Note: Some files were missing. See markers in the output.", file=sys.stderr)

if __name__ == "__main__":
    main()

