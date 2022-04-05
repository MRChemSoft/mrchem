#
# MRChem, a numerical real-space code for molecular electronic structure
# calculations within the self-consistent field (SCF) approximations of quantum
# chemistry (Hartree-Fock and Density Functional Theory).
# Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

import math
from types import SimpleNamespace

from qcelemental.physical_constants.context import PhysicalConstantsContext


class MRChemPhysConstants:
    """Class holding physical and math constants relevant to MRChem.
    The constants are fetched from QCElemental, and new constants are
    derived from those in QCElemental (as far as possible).

    Source in ascii:
    https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    """

    # Explicitly defined constants here
    # (not present in qcelemental or cannot be derived)
    PI = 3.1415926535897932384
    HARTREE2SIMAGNETIZABILITY = 78.9451185  # TODO: We should be able to derive this one

    def __init__(self):
        """Here we extract those constants that we need to run any MRChem calculation.
        Those that are not directly available in QCElemental, we compute by via the existing
        constants in QCElementlal.
        """
        self.context = "CODATA2018"
        self.nist = "https://physics.nist.gov/cuu/Constants/Table/allascii.txt"
        qce = PhysicalConstantsContext(context=self.context)

        # Define new constants needed by MRChem
        # We follow the QCElemental 4-tuple format
        # (callname: str, units: str, value: float, description: str)
        mrchem_constants = [
            ('pi',                              '',       self.PI,                                                                                         'Pi'),
            ('pi_sqrt',                         '',       math.sqrt(self.PI),                                                                              'Square root of pi'),
            ('hartree2simagnetizability',       'J T^-2', self.HARTREE2SIMAGNETIZABILITY,                                                                  'Atomic units to J/T^2 (magnetizability)'),
            ('atomic_unit_of_bohr_magneton',    '',       qce.Bohr_magneton / qce.atomic_unit_of_magnetizability / qce.atomic_unit_of_mag_flux_density,    'Bohr magneton in atomic units'),
            ('atomic_unit_of_nuclear_magneton', '',       qce.nuclear_magneton / qce.atomic_unit_of_magnetizability / qce.atomic_unit_of_mag_flux_density, 'Nuclear magneton in atomic units'),
            ('angstrom2bohrs',                  'Å',      1.0 / qce.bohr2angstroms,                                                                        'Angstrom -> Bohr conversion factor')
        ]

        # Add the following pre-defined constants to our mrchem subset
        # NOTE: when using the get method, we use the NIST access names
        # which may contain spaces. See web page for reference.
        names = [
            "hartree2kJmol",
            "hartree2kcalmol",
            "hartree2ev",
            "hartree2wavenumbers",
            "fine-structure constant",
            "c_au",
            "electron g factor",
            "dipmom_au2debye"
        ]

        # Append the data to our list of constants
        for name in names:
            datum = qce.get(name, return_tuple=True)
            constant = (datum.label, datum.units, datum.data, datum.comment)
            mrchem_constants.append(constant)

        # Set our constants to instance attributes
        for ident, _, value, _ in sorted(mrchem_constants, key=lambda x: x[0]):
            key = ident.lower().translate(qce._transtable)
            self.__setattr__(key, float(value))

    def print_constants_for_tests(self, varname='testConstants'):
        """Helper function for printing constants for copy/pasting into the c++ code.
        We need to store the constants internally for the tests to pass."""
        for key, value in self.__dict__.items():
            print(f'{varname}["{key}"] = {value};')

if __name__ == '__main__':
    MRChemPhysConstants().print_constants_for_tests()
