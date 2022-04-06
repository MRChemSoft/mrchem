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
import itertools
import re

from .periodictable import PeriodicTable, PeriodicTableByZ


class MoleculeValidator:
    """Sanity check routines for the user's Molecule section."""

    # Thresholds
    THRESHOLD_NUCLEAR_SINGULARITY_ERROR = 1.0e-6
    THRESHOLD_NUCLEAR_SINGULARITY_WARNING = 1.0e-3

    # Unit strings
    UNIT_ANGSTROM = 'angstrom'
    UNIT_BOHR = 'bohr'

    # Error/warning messages
    ERROR_MESSAGE_ATOMIC_COORDINATES =    lambda self, details: f'ABORT: INVALID ATOMIC COORDINATES: {details}'
    ERROR_MESSAGE_ATOMIC_SYMBOLS =        lambda self, details: f'ABORT: INVALID ATOMIC SYMBOLS: {details}'
    ERROR_MESSAGE_CAVITY_COORDINATES =    lambda self, details: f'ABORT: INVALID CAVITY COORDINATES: {details}'
    ERROR_MESSAGE_CAVITY_RADII =          lambda self, details: f'ABORT: INVALID CAVITY RADII: {details}'

    ERROR_MESSAGE_NUCLEAR_SINGULARITY =   lambda self, details: f'ABORT: SOME ATOMS TOO CLOSE (norm < {MoleculeValidator.THRESHOLD_NUCLEAR_SINGULARITY_ERROR}):\n{details}'
    WARNING_MESSAGE_NUCLEAR_SINGULARITY = lambda self, details: f'WARNING: SOME ATOMS VERY CLOSE (norm < {MoleculeValidator.THRESHOLD_NUCLEAR_SINGULARITY_WARNING}):\n{details}'

    ERROR_INCOMPATIBLE_MULTIPLICITY =     lambda self, details: f'ABORT: INCOMPATIBLE MULTIPLICITY: {details}'
    ERROR_UNPHYSICAL_MULTIPLICITY =       lambda self, details: f'ABORT: UNPHYSICAL MULTIPLICITY: {details}'
    ERROR_UNPHYSICAL_CHARGE =             lambda self, details: f'ABORT: UNPHYSICAL CHARGE: {details}'

    ERROR_RESTRICTED_OPEN_SHELL =         'ABORT: Restricted open-shell not implemented'


    def __init__(self, user_dict, origin):
        """
        Raises RunTimeError with helpful messages when invalid format
        or unphysical input value are detected:

        Atomic coordinates
        - Correct XYZ format checked with regexes (both atomic symbols 
          and numbers are valid atom identifiers, and can be used 
          interchancably in the same input)
        - Nuclear singularities
        - Atomic symbols checked against periodic table 
        
        Cavity spheres
        - Correct format checked with regexes
        - Negative radii not allowed

        Electrons, charge and multiplicity
        - charge > Z not allowed
        - n_unpaired > n_electrons not allowed
        - multiplicity and n_electrons both odd/even not allowed
        - restricted open-shell not allowed

        Unit conversions
        - convert to bohr if user specified angstrom

        Parameters
        ----------
        user_dict: dict, dictionary of user input
        origin: List[float], origin to be used in the calcuation
        """
        self.user_dict = user_dict
        self.origin = origin
        self.unit = user_dict['world_unit']
        self.pc = user_dict['Constants']

        # Molecule related data
        self.user_mol = user_dict['Molecule']
        self.charge = self.user_mol['charge']
        self.mult = self.user_mol['multiplicity']
        self.do_translate = self.user_mol['translate']
        self.coords_raw = self.user_mol['coords']

        # Cavity related data
        self.cavity_dict = user_dict['Environment']['Cavity']
        self.spheres_raw = self.cavity_dict['spheres']
        self.cavity_width = self.cavity_dict['cavity_width']
        self.has_sphere_input = len(self.spheres_raw.strip().splitlines()) > 0

        # Validate atomic coordinates
        self.atomic_symbols, self.atomic_coords = self.validate_atomic_coordinates()

        # Translate center of mass if requested
        # We must test for translation before validating the cavity, 
        # in case the nuclear coordinates are to be used for the
        # sphere centers
        if self.do_translate:
            self.atomic_coords = self.translate_com_to_origin()

        # Validate cavity spheres
        self.cavity_radii, self.cavity_coords = self.validate_cavity()

        # Perform some sanity checks
        self.check_for_nuclear_singularities()
        self.check_for_invalid_electronic_configuration()

        # Convert to bohrs if user gave angstroms
        if self.unit == self.UNIT_ANGSTROM:
            self.atomic_coords = self.ang2bohr_array(self.atomic_coords)
            self.cavity_coords = self.ang2bohr_array(self.cavity_coords)
            self.cavity_radii = self.ang2bohr_vector(self.cavity_radii)
            self.cavity_width *= self.pc['angstrom2bohrs']

    def get_coords_in_program_syntax(self):
        """Convert nuclear coordinates from JSON syntax to program syntax."""
        return [
            {
                "atom": label,
                "xyz": coord
            } for label, coord in zip(self.atomic_symbols, self.atomic_coords)
        ]

    def get_cavity_in_program_syntax(self):
        """Convert cavity spheres from JSON syntax to program syntax."""
        # Use sphere coordinates and radii if given
        if self.has_sphere_input:
            return [
                {
                    "center": center,
                    "radius": radius
                } for center, radius in zip(self.cavity_coords, self.cavity_radii)
            ]
        # If not build spheres from nuclear coordinates and default radii
        else:
            return [
                {
                    "center": coord,
                    "radius": PeriodicTable[label.lower()].radius
                } for label, coord in zip(self.atomic_symbols, self.atomic_coords)
            ]

    def validate_atomic_coordinates(self):
        """Parse the $coords block and ensure correct formatting."""
        # Regex components
        line_start = r'^'
        line_end = r'$'
        symbol = r'[a-zA-Z]{1,3}'
        decimal = r'[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)'
        integer = r'[0-9]+'
        one_or_more_whitespace = r'[\s]+'
        zero_or_more_whitespace = r'[\s]*'

        # Build regex
        atom_with_symbol = line_start + zero_or_more_whitespace + symbol + (one_or_more_whitespace + decimal)*3 + zero_or_more_whitespace + line_end
        atom_with_number = line_start + zero_or_more_whitespace + integer + (one_or_more_whitespace + decimal)*3 + zero_or_more_whitespace + line_end

        p_with_symbol = re.compile(atom_with_symbol)
        p_with_number = re.compile(atom_with_number)

        # Parse coordinates
        coords = []
        labels = []
        bad_atoms = []
        for atom in self.coords_raw.strip().splitlines():
            match_symbol = p_with_symbol.match(atom)
            match_number = p_with_number.match(atom)
            if match_symbol:
                g = match_symbol.group()
                labels.append(g.split()[0].strip().lower())
                coords.append([float(c.strip()) for c in g.split()[1:]])
            elif match_number:
                g = match_number.group()
                Z = int(g.split()[0].strip())
                labels.append(PeriodicTableByZ[Z].symbol.lower())
                coords.append([float(c.strip()) for c in g.split()[1:]])
            else:
                bad_atoms.append(atom)

        if bad_atoms:
            newline = '\n'
            raise RuntimeError(self.ERROR_MESSAGE_ATOMIC_COORDINATES(
                f'One or more atomic coordinates had an invalid input format:\n{newline.join(bad_atoms)}'
            ))

        # Check that the atomic symbols represent valid elements
        fltr = filter(lambda x: x not in PeriodicTable, labels)
        if any(list(fltr)):
            newline = '\n'
            raise RuntimeError(self.ERROR_MESSAGE_ATOMIC_SYMBOLS(
                f'One or more invalid atomic symbols:\n{newline.join(set(fltr))}'
            ))

        return labels, coords


    def validate_cavity(self):
        """Parse the $spheres block and ensure correct formatting."""
        # Regex components
        line_start = r'^'
        line_end = r'$'
        decimal = r'[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)'
        one_or_more_whitespace = r'[\s]+'
        zero_or_more_whitespace = r'[\s]*'

        # Build regex
        valid_sphere = line_start + zero_or_more_whitespace + decimal + (one_or_more_whitespace + decimal)*3 + zero_or_more_whitespace + line_end
        p = re.compile(valid_sphere)

        # Parse spheres
        coords = []
        radii = []
        bad_spheres = []

        for sphere in self.spheres_raw.strip().splitlines():
            match = p.match(sphere)
            if match:
                g = match.group()
                radii.append(float(g.split()[-1].strip()))
                coords.append([float(c.strip()) for c in g.split()[:-1]])
            else:
                bad_spheres.append(sphere)

        if bad_spheres:
            newline = '\n'
            raise RuntimeError(self.ERROR_MESSAGE_CAVITY_COORDINATES(
                f'One or more cavity spheres had an invalid input format:\n{newline.join(bad_spheres)}'
            ))

        # Check for negative radii
        if any([r < 0 for r in radii]):
            raise RuntimeError(self.ERROR_MESSAGE_CAVITY_RADII(
                'Cavity radii cannot be negative'
            ))

        return radii, coords

    def check_for_nuclear_singularities(self):
        """Check for singularities in the nuclear positions."""
        # Bad pairs will be stored here
        error_pairs = []
        warning_pairs = []

        # Loop over all unique atom pairs and compute euclidian distance
        for (ca, la), (cb, lb) in itertools.combinations(zip(self.atomic_coords, self.atomic_symbols), 2):
            pair_label = f'{la}: {ca}\n{lb}: {cb}'
            R = self.euclidian_distance(ca, cb)

            # Compare distance to internal thresholds
            if R < self.THRESHOLD_NUCLEAR_SINGULARITY_ERROR:
                error_pairs.append(pair_label)
            elif R < self.THRESHOLD_NUCLEAR_SINGULARITY_WARNING:
                warning_pairs.append(pair_label)

        # Print warnings and raise exception if necessary
        if warning_pairs:
            msg = self.WARNING_MESSAGE_NUCLEAR_SINGULARITY('\n\n'.join(warning_pairs))
            print(msg)

        if error_pairs:
            msg = self.ERROR_MESSAGE_NUCLEAR_SINGULARITY('\n\n'.join(error_pairs))
            raise RuntimeError(msg)
                
    def check_for_invalid_electronic_configuration(self):
        """Check that the number of electrons and spin multiplicity are compatible.
        Also check for restricted open-shell calculation."""
        restricted = self.user_dict['WaveFunction']['restricted']
        Z = sum([PeriodicTable[atom.lower()].Z for atom in self.atomic_symbols])
        n_electrons = Z - self.charge
        n_unpaired = self.mult - 1

        # Helper function
        parity = lambda n: 'even' if n % 2 == 0 else 'odd'

        # Check for impossible charge
        if self.charge > Z:
            raise RuntimeError(self.ERROR_UNPHYSICAL_CHARGE(
                f'The specified charge ({self.charge}) cannot be larger than the nuclear charge ({Z})'
            ))

        # Check for unphysical multiplicity
        elif n_unpaired > n_electrons:
            raise RuntimeError(self.ERROR_UNPHYSICAL_MULTIPLICITY(
                f"The specified multiplicity requires more unpaired electrons ({self.mult - 1}) than are available ({n_electrons}))."
                ))

        # Check for invalid spin multiplicity
        elif parity(n_electrons) == parity(self.mult): 
            raise RuntimeError(self.ERROR_INCOMPATIBLE_MULTIPLICITY(
                f"The specified multiplicity ({parity(self.mult)}) is not compatible with the number of electrons ({parity(n_electrons)})"
                ))

        # Check for restricted open-shell
        elif restricted and n_unpaired > 0: 
            raise RuntimeError(self.ERROR_RESTRICTED_OPEN_SHELL)

    def translate_com_to_origin(self):
        """Translate center of mass to the origin."""
        masses = [PeriodicTable[label.lower()].mass for label in self.atomic_symbols]
        M = sum(masses)
        return [masses[i] * (self.atomic_coords[i] - self.origin[i]) / M for i in range(3)]

    @staticmethod
    def euclidian_distance(a, b):
        """Helper function for the nuclear singularies validation.
        Computes the euclidian distance between two vectors, a and b."""
        squared_deviations = [(a[i] - b[i])**2 for i in range(3)]
        return math.sqrt(sum(squared_deviations))

    def ang2bohr_array(self, coords):
        """Convert List[List[float]] from angstrom to bohr."""
        return [
            [c * self.pc['angstrom2bohrs'] for c in element] for element in coords
        ]

    def ang2bohr_vector(self, vec):
        """Convert List[float] from angstrom to bohr"""
        return [el * self.pc['angstrom2bohrs'] for el in vec]
