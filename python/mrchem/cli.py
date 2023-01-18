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

import argparse

from .config import MRCHEM_EXECUTABLE, MRCHEM_VERSION, MRCHEM_MODULE


def cli():
    cli = argparse.ArgumentParser(description="Front-end CLI for MRChem")
    cli.add_argument("-v", "--version", action="version", version=MRCHEM_VERSION)
    cli.add_argument(
        "--launcher",
        action="store",
        dest="launcher",
        type=str,
        default="",
        help="Set program launcher string",
    )
    cli.add_argument(
        "--executable",
        "-x",
        action="store",
        dest="executable",
        type=str,
        default=MRCHEM_EXECUTABLE,
        help="Set executable name",
    )
    cli.add_argument(
        "--dryrun",
        "-D",
        action="store_true",
        dest="dryrun",
        default=False,
        help="Only process input",
    )
    cli.add_argument(
        "--stdout",
        action="store_true",
        dest="stdout",
        default=False,
        help="Print to stdout",
    )
    cli.add_argument(
        "--json",
        "-j",
        action="store_true",
        dest="inp_json",
        default=False,
        help="Input file is in json format",
    )
    cli.add_argument("--module", "-m", action="version", version=str(MRCHEM_MODULE))
    cli.add_argument("inp_name", type=str, help="Input file in getkw format")

    args = cli.parse_args()

    return args
