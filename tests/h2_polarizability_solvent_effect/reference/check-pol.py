#!/usr/bin/env python3

from json import load
from pathlib import Path

dipoles = []
ef_strength = 0.01

curdir = Path(__file__).parent.resolve()

for out in ["h2_m", "h2_p"]:
    with (curdir / f"{out}.json").open("r") as fh:
        data = load(fh)
        dipoles.append(
            data["output"]["properties"]["dipole_moment"]["dip-1"]["vector_el"][-1]
        )

fd_pol = (dipoles[1] - dipoles[0]) / (2 * ef_strength)

with (curdir / "h2.json").open("r") as fh:
    data = load(fh)
    pol = data["output"]["properties"]["polarizability"]["pol-0.000000"]["tensor"][-1]

print(f"Finite difference: {fd_pol}")
print(f"Analytical: {pol}")
