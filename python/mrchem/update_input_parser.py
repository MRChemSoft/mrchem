#!/usr/bin/env python

#
# MRChem, a numerical real-space code for molecular electronic structure
# calculations within the self-consistent field (SCF) approximations of quantum
# chemistry (Hartree-Fock and Density Functional Theory).
# Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
#
# This file is part of MRChem.
#
# MRChem is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MRChem is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
#
# For information on the complete list of contributors to MRChem, see:
# <https://mrchem.readthedocs.io/>
#


from pathlib import Path
from ruamel.yaml import YAML
import subprocess
import os
import shutil
import argparse

from physical_constants import MRChemPhysConstants

root = Path.cwd().parent
target = root.joinpath("mrchem", "input_parser")

yaml = YAML()
yaml.indent(sequence=4, mapping=2, offset=2)


def update_constants():
    pc = MRChemPhysConstants()
    f_template = root.joinpath("template.yml")
    template = yaml.load(f_template)

    new = {
        "keywords": template["keywords"],
        "sections": [
            section
            for section in template["sections"]
            if section["name"] != "Constants"
        ],
    }

    # Build new constants section
    constants = {
        "name": "Constants",
        "docstring": "Physical and mathematical constants used by MRChem",
        "keywords": [
            {"name": name, "default": value, "type": "float", "docstring": docstring}
            for name, _, value, docstring in pc.data
        ],
    }

    # Construct the full template file and dump
    new["sections"].append(constants)
    yaml.dump(new, f_template)


def run_parselglossy():
    os.chdir(root)
    cmd = [
        "parselglossy",
        "generate",
        "--template",
        "template.yml",
        "--docfile",
        "user_ref.rst",
        "--doc-header",
        "User input reference",
        "--target",
        target,
    ]

    subprocess.call(cmd)


def update_doc():
    src = target.joinpath("docs/user_ref.rst")
    dst = root.parent.joinpath("doc/users/user_ref.rst")
    shutil.copyfile(src, dst)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CLI for synchronizing the template with current constants, plus updating the input parser and user ref."
    )
    parser.add_argument(
        "-sT",
        "--skip-template",
        action="store_true",
        help="Do not update constants section in template.yml",
    )
    parser.add_argument(
        "-sP",
        "--skip-parselglossy",
        action="store_true",
        help="Do not update input parser with parselglossy",
    )
    parser.add_argument(
        "-sD",
        "--skip-doc",
        action="store_true",
        help="Do not update the user reference file",
    )
    parser.add_argument("-D", "--dryrun", action="store_true", help="Skip all actions")
    args = parser.parse_args()

    if args.dryrun:
        args.skip_template = True
        args.skip_parselglossy = True
        args.skip_doc = True

    if not args.skip_template:
        print(f'{"Updating template":20} ... ', end="")
        update_constants()
        print("done")

    if not args.skip_parselglossy:
        print(f'{"Running parselglossy":20} ... ', end="")
        run_parselglossy()
        print("done")

    if not args.skip_doc:
        print(f'{"Updating user ref":20} ... ', end="")
        update_doc()
        print("done")
