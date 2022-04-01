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

from qcelemental.physical_constants.context import PhysicalConstantsContext
from qcelemental.datum import Datum


class MRChemPhysConstants(PhysicalConstantsContext):
    """Wrapper over the PhysicalConstantsContext class from QCElemental.
    Subclassing it here allows for some customization, and it ensures that
    when imported the same CODATA source is used automatically (we use 2018).
    """
    def __init__(self, context="CODATA2018"):
        super().__init__(context)

        # Define custom shorthands. Each tuple is organized as follows (to be compatible with the Datum object):
        # (callname, units, value, description)
        customs = [
            ('pi',                              '',       3.1415926535897932384,                                                                              'Pi'),
            ('pi_sqrt',                         '',       1.7724538509055160273,                                                                              'Square root of pi'),
            ('hartree2simagnetizability',       'J T^-2', 78.9451185,                                                                                         'Atomic units to J/T^2 (magnetizability)'),
            ('atomic_unit_of_bohr_magneton',    '',       self.Bohr_magneton / self.atomic_unit_of_magnetizability / self.atomic_unit_of_mag_flux_density,    'Bohr magneton in atomic units'),
            ('atomic_unit_of_nuclear_magneton', '',       self.nuclear_magneton / self.atomic_unit_of_magnetizability / self.atomic_unit_of_mag_flux_density, 'Nuclear magneton in atomic units'),
            ('angstrom2bohrs',                  'Å',      1.0 / self.bohr2angstroms,                                                                           'Angstrom -> Bohr conversion factor')
        ]

        # Add aliases to internally stored constants and make them callable
        for ident, units, value, comment in customs:
            self.pc[ident.lower()] = Datum(ident, units, value, comment=comment)
            self.__setattr__(ident, value)

    def to_dict(self):
        """Generate a dictionary storing the constants, which can be passed to the c++ program.
        
        The transtable is used to make the names Python-friendly."""
        return {
            qca.label.lower().translate(self._transtable): float(qca.data) for qca in self.pc.values()
        }

    def print_subset_for_unit_tests(self, varname="testConstants"):
        """Helper function for printing a subset of the constants."""
        subset = [
            "pi",
            "pi_sqrt",
            "electron_g_factor",
            "fine_structure_constant",
            "hartree2kJmol",
            "hartree2kcalmol",
            "hartree2ev",
            "hartree2simagnetizability"
        ]

        content = [
            f'{varname}["{c.lower()}"] = {self.__getattribute__(c)};' for c in subset
        ]

        print('\n'.join(content))


if __name__ == '__main__':
    c = MRChemPhysConstants()
    c.print_subset_for_unit_tests()