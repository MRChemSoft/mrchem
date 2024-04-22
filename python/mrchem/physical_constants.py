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

import math
from pydoc import doc
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
        self.qce = PhysicalConstantsContext(context=self.context)
        self.data = []

        ########################
        # Define our constants
        ########################

        # Convert from au to SI units for magnetizability
        name = "hartree2simagnetizability"
        unit = "J T^-2"
        value = self.HARTREE2SIMAGNETIZABILITY
        docstring = f"| Conversion factor for magnetizability from atomic units to SI units  (unit: {unit}). Affected code: Printed value of the magnetizability property."
        self.add_constant(name, unit, value, docstring)

        # Speed of light in atomic units
        name = "light_speed"
        unit = "au"
        value = self.qce.c_au
        docstring = f"| Speed of light in atomic units  (unit: {unit}). Affected code: Relativistic Hamiltonians (ZORA, etc.)"
        self.add_constant(name, unit, value, docstring)

        # Convert from Angstrom to Bohr
        name = "angstrom2bohrs"
        unit = "Å^-1"
        value = 1.0 / self.qce.bohr2angstroms
        docstring = f"| Conversion factor for Cartesian coordinates from Angstrom to Bohr  (unit: {unit}). Affected code: Parsing of input coordinates, printed coordinates"
        self.add_constant(name, unit, value, docstring)

        # Convert from Hartree to kJ/mol
        name = "hartree2kjmol"
        unit = "kJ mol^-1"
        value = self.qce.hartree2kJmol
        docstring = f"| Conversion factor from Hartree to kJ/mol  (unit: {unit}). Affected code: Printed value of energies."
        self.add_constant(name, unit, value, docstring)

        # Convert from Hartree to kcal/mol
        name = "hartree2kcalmol"
        unit = "kcal mol^-1"
        value = self.qce.hartree2kcalmol
        docstring = f"| Conversion factor from Hartree to kcal/mol  (unit: {unit}). Affected code: Printed value of energies."
        self.add_constant(name, unit, value, docstring)

        # Convert from Hartree to eV
        name = "hartree2ev"
        unit = "ev"
        value = self.qce.hartree2ev
        docstring = (
            f"| Conversion factor from Hartree to eV  (unit: {unit}). Affected code: Printed value of energies."
            ""
        )
        self.add_constant(name, unit, value, docstring)

        # Convert from Hartree to cm-1
        name = "hartree2wavenumbers"
        unit = "cm^-1"
        value = self.qce.hartree2wavenumbers
        docstring = f"| Conversion factor from Hartree to wavenumbers (unit: {unit}). Affected code: Printed value of frequencies."
        self.add_constant(name, unit, value, docstring)

        # Fine-structure constant in atomic units
        name = "fine_structure_constant"
        unit = "au"
        value = self.qce.fine_structure_constant
        docstring = f"| Fine-structure constant in atomic units (unit: {unit}). Affected code: Certain magnetic interaction operators."
        self.add_constant(name, unit, value, docstring)

        # Electron g factor in atomic units
        name = "electron_g_factor"
        unit = "au"
        value = self.qce.electron_g_factor
        docstring = f"| Electron g factor in atomic units (unit: {unit}). Affected code: Certain magnetic interaction operators."
        self.add_constant(name, unit, value, docstring)

        # Convert from atomic units to Debye
        name = "dipmom_au2debye"
        unit = "?"
        value = self.qce.dipmom_au2debye
        docstring = f"| Conversion factor for dipoles from atomic units to Debye (unit: {unit}). Affected code: Printed value of dipole moments."
        self.add_constant(name, unit, value, docstring)

        # Boltzmann constant in J K^{-1}
        name = "boltzmann_constant"
        unit = "J K^-1"
        value = self.qce.kb
        docstring = f"| Boltzmann constant in (unit: {unit}). Affected code: Value of the Debye-Huckel screening parameter in the Poisson-Boltzmann equation."
        self.add_constant(name, unit, value, docstring)

        # elementary charge in C
        name = "elementary_charge"
        unit = "C"
        value = self.qce.get("elementary charge")
        docstring = f"| Elementary charge in (unit: {unit}). Affected code: Value of the Debye-Huckel screening parameter in the Poisson-Boltzmann equation."
        self.add_constant(name, unit, value, docstring)

        # Permittivity of free space in F m^-1
        name = "e0"
        unit = "F m^-1"
        value = self.qce.e0
        docstring = f"| Permittivity of free space (unit: {unit}). Affected code: Value of the Debye-Huckel screening parameter in the Poisson-Boltzmann equation."
        self.add_constant(name, unit, value, docstring)

        # Avogadro constant in mol^-1
        name = "N_a"
        unit = "mol^-1"
        value = self.qce.na
        docstring = f"| Avogadro constant (unit: {unit}). Affected code: Value of the Debye-Huckel screening parameter in the Poisson-Boltzmann equation."
        self.add_constant(name, unit, value, docstring)

        # meter to Bohr radius
        name = "meter2bohr"
        unit = "m^-1"
        value = 1 / self.qce.bohr2m
        docstring = f"| conversion factor from meter to Bohr radius (unit: {unit}). Affected code: Value of the Debye-Huckel screening parameter in the Poisson-Boltzmann equation."
        self.add_constant(name, unit, value, docstring)

        # Set our constants to instance attributes for easy access
        for name, _, value, _ in self.data:
            self.__setattr__(name.lower(), float(value))

    def add_constant(self, name, unit, value, docstring):
        """Add constant to `data` instance attribute."""
        self.data.append((name, unit, value, docstring))

    def print_constants_for_tests(self, varname="testConstants"):
        """Helper function for printing constants for copy/pasting into the c++ code.
        We need to store the constants internally for the unit tests to pass."""
        for name, _, value, _ in sorted(self.data, key=lambda x: x[0]):
            print(f'{{"{name}", {value}}},')


if __name__ == "__main__":
    MRChemPhysConstants().print_constants_for_tests()
