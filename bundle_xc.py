from datetime import datetime
from pathlib import Path

out = f"xc_bundle_{datetime.now():%Y%m%d_%H%M%S}.txt"
files = [
    # 1) The XC “chokepoint”
    "src/mrdft/xc_utils.h",
    "src/mrdft/xc_utils.cpp",

    # 2) Where functionals are created/parsed
    "src/mrdft/Functional.h",
    "src/mrdft/Functional.cpp",
    "src/mrdft/Factory.h",
    "src/mrdft/Factory.cpp",

    # 3) Callers that pass rho/sigma/tau/spin into the XC layer
    "src/mrdft/LDA.h",
    "src/mrdft/LDA.cpp",
    "src/mrdft/GGA.h",
    "src/mrdft/GGA.cpp",
    "src/mrdft/SpinLDA.h",
    "src/mrdft/SpinLDA.cpp",
    "src/mrdft/SpinGGA.h",
    "src/mrdft/SpinGGA.cpp",

    # 4) Operators / properties that might call XC directly
    "src/qmoperators/xc_operator_lda.cpp",
    "src/qmoperators/xc_operator_blyp.cpp",
    "src/qmoperators/xc_hessian_lda.cpp",
    "src/qmoperators/xc_hessian_pbe.cpp",
    "src/surface_forces/xcStress.h",
    "src/surface_forces/xcStress.cpp",
    "src/qmoperators/two_electron/XCPotential.h",
    "src/qmoperators/two_electron/XCPotential.cpp",
    "src/qmoperators/two_electron/XCPotentialD1.h",
    "src/qmoperators/two_electron/XCPotentialD1.cpp",
    "src/qmoperators/two_electron/XCPotentialD2.h",
    "src/qmoperators/two_electron/XCPotentialD2.cpp",

    # 5) Build system bits I’ll patch
    "CMakeLists.txt",
    "src/mrdft/CMakeLists.txt",
    "external/upstream/fetch_xcfun.cmake",
    "cmake/custom/main.cmake",
]

with open(out, "w", encoding="utf-8", errors="replace") as o:
    for f in files:
        p = Path(f)
        if p.is_file():
            o.write(f"===== BEGIN {f} =====\n")
            try:
                o.write(p.read_text(encoding="utf-8"))
            except UnicodeDecodeError:
                # fallback for odd encodings
                o.write(p.read_bytes().decode("utf-8", errors="replace"))
            o.write(f"\n===== END {f} =====\n\n")
        else:
            o.write(f"===== MISSING {f} =====\n\n")

print(f"Wrote {out}")

