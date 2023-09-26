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

import itertools
import math
import re
from copy import deepcopy

from .periodictable import PeriodicTable, PeriodicTableByZ


class MoleculeValidator:
    """Sanity check routines for the user's Molecule section."""

    # Thresholds
    THRESHOLD_NUCLEAR_SINGULARITY_ERROR = 1.0e-6
    THRESHOLD_NUCLEAR_SINGULARITY_WARNING = 1.0e-3

    # Unit strings
    UNIT_ANGSTROM = "angstrom"
    UNIT_BOHR = "bohr"

    # Error/warning messages
    ERROR_MESSAGE_ATOMIC_COORDINATES = (
        lambda self, details: f"ABORT: INVALID ATOMIC COORDINATES: {details}"
    )
    ERROR_MESSAGE_ATOMIC_SYMBOLS = (
        lambda self, details: f"ABORT: INVALID ATOMIC SYMBOLS: {details}"
    )
    ERROR_MESSAGE_CAVITY_COORDINATES = (
        lambda self, details: f"ABORT: INVALID CAVITY COORDINATES: {details}"
    )
    ERROR_MESSAGE_CAVITY_RADII = (
        lambda self, details: f"ABORT: INVALID CAVITY RADII: {details}"
    )
    ERROR_MESSAGE_CAVITY_ALPHAS = (
        lambda self, details: f"ABORT: INVALID CAVITY SCALING FACTORS: {details}"
    )

    ERROR_MESSAGE_NUCLEAR_SINGULARITY = (
        lambda self, details: f"ABORT: SOME ATOMS TOO CLOSE (norm < {MoleculeValidator.THRESHOLD_NUCLEAR_SINGULARITY_ERROR}):\n{details}"
    )
    WARNING_MESSAGE_NUCLEAR_SINGULARITY = (
        lambda self, details: f"WARNING: SOME ATOMS VERY CLOSE (norm < {MoleculeValidator.THRESHOLD_NUCLEAR_SINGULARITY_WARNING}):\n{details}"
    )

    ERROR_INCOMPATIBLE_MULTIPLICITY = (
        lambda self, details: f"ABORT: INCOMPATIBLE MULTIPLICITY: {details}"
    )
    ERROR_UNPHYSICAL_MULTIPLICITY = (
        lambda self, details: f"ABORT: UNPHYSICAL MULTIPLICITY: {details}"
    )
    ERROR_UNPHYSICAL_CHARGE = (
        lambda self, details: f"ABORT: UNPHYSICAL CHARGE: {details}"
    )

    ERROR_RESTRICTED_OPEN_SHELL = "ABORT: Restricted open-shell not implemented"

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
        self.unit = user_dict["world_unit"]
        self.pc = user_dict["Constants"]

        # Molecule related data
        self.user_mol = user_dict["Molecule"]
        self.charge = self.user_mol["charge"]
        self.mult = self.user_mol["multiplicity"]
        self.do_translate = self.user_mol["translate"]
        self.coords_raw = self.user_mol["coords"]

        # Cavity related data
        self.cavity_dict = user_dict["PCM"]["Cavity"]
        self.cavity_mode = self.cavity_dict["mode"]
        self.spheres_raw = self.cavity_dict["spheres"]
        self.cavity_alpha = self.cavity_dict["alpha"]
        self.cavity_beta = self.cavity_dict["beta"]
        self.cavity_sigma = self.cavity_dict["sigma"]

        # Validate atomic coordinates
        (
            self.atomic_symbols,
            self.atomic_coords,
            self.atomic_rms,
        ) = self.validate_atomic_coordinates()
        self.n_atoms = len(self.atomic_coords)

        # Translate center of mass if requested
        # We must test for translation before validating the cavity,
        # in case the nuclear coordinates are to be used for the
        # sphere centers
        if self.do_translate:
            self.translate_com_to_origin()

        # Validate cavity spheres
        (
            self.cavity_radii,
            self.cavity_coords,
            self.cavity_alphas,
            self.cavity_betas,
            self.cavity_sigmas,
        ) = self.validate_cavity()

        # Perform some sanity checks
        self.check_for_nuclear_singularities()
        self.check_for_invalid_electronic_configuration()

        # Convert to bohrs if user gave angstroms
        if self.unit == self.UNIT_ANGSTROM:
            self.atomic_coords = self.ang2bohr_array(self.atomic_coords)
            self.cavity_coords = self.ang2bohr_array(self.cavity_coords)
            self.cavity_radii = self.ang2bohr_vector(self.cavity_radii)
            self.cavity_sigmas = self.ang2bohr_vector(self.cavity_sigmas)

    def get_coords_in_program_syntax(self):
        """Convert nuclear coordinates from JSON syntax to program syntax."""
        return [
            {"atom": label, "xyz": coord, "r_rms": rms}
            for label, coord, rms in zip(
                self.atomic_symbols, self.atomic_coords, self.atomic_rms
            )
        ]

    def get_cavity_in_program_syntax(self):
        """Convert cavity spheres from JSON syntax to program syntax."""
        return [
            {
                "center": center,
                "radius": radius,
                "alpha": alpha,
                "beta": beta,
                "sigma": sigma,
            }
            for center, radius, alpha, beta, sigma in zip(
                self.cavity_coords,
                self.cavity_radii,
                self.cavity_alphas,
                self.cavity_betas,
                self.cavity_sigmas,
            )
        ]

    def validate_atomic_coordinates(self):
        """Parse the $coords block and ensure correct formatting."""
        # Regex components
        line_start = r"^"
        line_end = r"$"
        symbol = r"[a-zA-Z]{1,3}"
        decimal = r"[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)"
        integer = r"[0-9]+"
        one_or_more_whitespace = r"[\s]+"
        zero_or_more_whitespace = r"[\s]*"

        # Build regex
        atom_with_symbol = (
            line_start
            + zero_or_more_whitespace
            + symbol
            + (one_or_more_whitespace + decimal) * 3
            + zero_or_more_whitespace
            + line_end
        )
        atom_with_number = (
            line_start
            + zero_or_more_whitespace
            + integer
            + (one_or_more_whitespace + decimal) * 3
            + zero_or_more_whitespace
            + line_end
        )

        p_with_symbol = re.compile(atom_with_symbol)
        p_with_number = re.compile(atom_with_number)

        lines = [
            x.strip() for x in self.coords_raw.strip().splitlines() if x != ""
        ]
        # Parse coordinates
        coords = []
        labels = []
        radii = []
        bad_atoms = []
        for atom in lines:
            match_symbol = p_with_symbol.match(atom)
            match_number = p_with_number.match(atom)
            if match_symbol:
                g = match_symbol.group()
                symbol = g.split()[0].strip().lower()
                labels.append(symbol)
                radii.append(float(PeriodicTable[symbol].r_rms))
                coords.append([float(c.strip()) for c in g.split()[1:]])
            elif match_number:
                g = match_number.group()
                Z = int(g.split()[0].strip())
                labels.append(PeriodicTableByZ[Z].symbol.lower())
                radii.append(float(PeriodicTableByZ[Z].r_rms))
                coords.append([float(c.strip()) for c in g.split()[1:]])
            else:
                bad_atoms.append(atom)

        if bad_atoms:
            newline = "\n"
            raise RuntimeError(
                self.ERROR_MESSAGE_ATOMIC_COORDINATES(
                    f"One or more atomic coordinates had an invalid input format:\n{newline.join(bad_atoms)}"
                )
            )

        # Check that the atomic symbols represent valid elements
        fltr = filter(lambda x: x not in PeriodicTable, labels)
        if any(list(fltr)):
            newline = "\n"
            raise RuntimeError(
                self.ERROR_MESSAGE_ATOMIC_SYMBOLS(
                    f"One or more invalid atomic symbols:\n{newline.join(set(fltr))}"
                )
            )

        return labels, coords, radii

    def validate_cavity(self):
        """Parse the $spheres block and ensure correct formatting.

        - Each sphere "object" in the JSON will contain:

          * Position in Cartesian coordinates x_i, y_i, z_i
          * **Unscaled** radius R_i^vdW
          * Radius scaling factor alpha_i
          * Width scaling factor beta_i
          * Width sigma_i

          Inside MRChem, the radius to be used will be computed as
          R_i = alpha_i*R_i^vdW + beta_i*sigma_i.
        - In `atoms` mode, we first create the list of spheres from the atomic
          coordinates and built-in radii, alpha (1.1), beta (0.5), and
          sigma (0.2) values. We then read the `$spheres`/`$end` blob of text
          and amend this list. Spheres can be given in two syntaxes:

          * **substitution**: `index radius [alpha] [beta] [sigma]`, and
          * **addition**: `x y z radius [alpha] [beta] [sigma]`

          These lines are parsed by two regexes. The "substitution" regex (`p`
          in the code) has higher priority: it is run and post-processed first.
          Note that the `index` in "substitution" syntax **starts from 1**, but
          it's easy to go to 0-based indexing.
        - In `explicit` mode, we start from an empty list of spheres. We then
          read the `$spheres`/`$end` blob of text and amend this list. Spheres can
          be given in two syntaxes:

          * **centered**: `index radius [alpha] [beta] [sigma]`, and
          * **free**: `x y z radius [alpha] [beta] [sigma]`

          These lines are parsed by the same two regexes, with same priority and
          ordering considerations, as above.
        - "Centered" and "substituted" spheres can be aware of their parent
          atom (I believe this needs some change in the C++ code though!) which
          means we can track their contribution to the molecular gradient. "Added"
          and "free" spheres are always treated as not atom-centered, so they will
          never contribute to the molecular gradient.
        - After discussions with @ilfreddy, I decided that it's less surprising
          for the user if the default values of `alpha`, `beta`, and `sigma`
          are **always** applied. If special behavior is desired, then it needs to
          be requested *explicitly* in the input.
        """

        # Regex components
        integer = r"[0-9]+"
        decimal = r"[+-]?[0-9]+\.?[0-9]*|\.[0-9]+"
        positive_decimal = r"[0-9]+\.?[0-9]*|\.[0-9]+"

        # Build regexes
        valid_atom = rf"^(?P<index>{integer})\s+(?P<radius>{positive_decimal})(?:\s+)?(?P<alpha>{positive_decimal})?(?:\s+)?(?P<beta>{positive_decimal})?(?:\s+)?(?P<sigma>{positive_decimal})?$"
        p = re.compile(valid_atom)
        valid_sphere = rf"^(?P<X>{decimal})\s+(?P<Y>{decimal})\s+(?P<Z>{decimal})\s+(?P<radius>{positive_decimal})(?:\s+)?(?P<alpha>{positive_decimal})?(?:\s+)?(?P<beta>{positive_decimal})?(?:\s+)?(?P<sigma>{positive_decimal})?$"
        q = re.compile(valid_sphere)

        # the centers of the spheres are the same as the atoms
        if self.cavity_mode == "atoms":

            def radiusNotFound(x):
                # This raises an exception in the list comprehension if the radius is not valid.
                raise ValueError(
                    f"The vdw-radius of element {x} is not defined in the mantina set"
                )

            coords = deepcopy(self.atomic_coords)
            # fetches mantina radii from the template.yml file for each atom x.
            # If the radius is negative, it means it is not defined in the
            # mantina set and it raises an exception.
            radii = [
                self.user_dict["Elements"][x.lower()]["vdw-radius"]
                if (self.user_dict["Elements"][x.lower()]["vdw-radius"] > 0.0)
                else radiusNotFound(x)
                for x in self.atomic_symbols
            ]
            alphas = [self.cavity_alpha] * len(radii)
            betas = [self.cavity_beta] * len(radii)
            sigmas = [self.cavity_sigma] * len(radii)

        else:
            coords = []
            radii = []
            alphas = []
            betas = []
            sigmas = []

        # Parse spheres
        bad_spheres = []

        lines = [
            x.strip() for x in self.spheres_raw.strip().splitlines() if x != ""
        ]
        for sphere in lines:
            p_match = p.match(sphere)
            q_match = q.match(sphere)
            if p_match:
                # the indexing of the atoms is 0-based in the user input!
                index = int(p_match.group("index"))

                if self.cavity_mode == "atoms":
                    radii[index] = float(p_match.group("radius"))

                    if p_match.group("alpha"):
                        alphas[index] = float(p_match.group("alpha"))

                    if p_match.group("beta"):
                        betas[index] = float(p_match.group("beta"))

                    if p_match.group("sigma"):
                        sigmas[index] = float(p_match.group("sigma"))
                else:
                    coords.append(self.atomic_coords[index])
                    radii.append(float(p_match.group("radius")))
                    alphas.append(
                        (
                            float(p_match.group("alpha"))
                            if p_match.group("alpha")
                            else self.cavity_alpha
                        )
                    )
                    betas.append(
                        (
                            float(p_match.group("beta"))
                            if p_match.group("beta")
                            else self.cavity_beta
                        )
                    )
                    sigmas.append(
                        (
                            float(p_match.group("sigma"))
                            if p_match.group("sigma")
                            else self.cavity_sigma
                        )
                    )
            elif q_match:
                coords.append(
                    [
                        float(q_match.group("X")),
                        float(q_match.group("Y")),
                        float(q_match.group("Z")),
                    ]
                )
                radii.append(float(q_match.group("radius")))
                alphas.append(
                    (
                        float(q_match.group("alpha"))
                        if q_match.group("alpha")
                        else self.cavity_alpha
                    )
                )
                betas.append(
                    (
                        float(q_match.group("beta"))
                        if q_match.group("beta")
                        else self.cavity_beta
                    )
                )
                sigmas.append(
                    (
                        float(q_match.group("sigma"))
                        if q_match.group("sigma")
                        else self.cavity_sigma
                    )
                )
            else:
                bad_spheres.append(sphere)

        if bad_spheres:
            newline = "\n"
            raise RuntimeError(
                self.ERROR_MESSAGE_CAVITY_COORDINATES(
                    f"One or more cavity spheres had an invalid input format:\n{newline.join(bad_spheres)}"
                )
            )

        # Check for negative or zero radii
        invalid_radii = {
            i: r
            for i, r in enumerate(radii)
            if ((r < 0) or math.isclose(r, 0.0))
        }
        if invalid_radii:
            invalid = "\n".join(
                [
                    f"Sphere {i} has invalid radius {r}"
                    for i, r in invalid_radii.items()
                ]
            )
            raise RuntimeError(
                self.ERROR_MESSAGE_CAVITY_RADII(
                    f"Cavity radii cannot be negative or zero:\n{invalid}"
                )
            )

        # Check for negative or zero scaling factors
        invalid_alphas = {
            i: a
            for i, a in enumerate(alphas)
            if ((a < 0) or math.isclose(a, 0.0))
        }
        if invalid_alphas:
            invalid = "\n".join(
                [
                    f"Sphere {i} has invalid radius {a}"
                    for i, a in invalid_alphas.items()
                ]
            )
            raise RuntimeError(
                self.ERROR_MESSAGE_CAVITY_ALPHAS(
                    f"Cavity radii scaling factors cannot be negative or zero:\n{invalid}"
                )
            )

        return radii, coords, alphas, betas, sigmas

    def check_for_nuclear_singularities(self):
        """Check for singularities in the nuclear positions."""
        # Bad pairs will be stored here
        error_pairs = []
        warning_pairs = []

        # Loop over all unique atom pairs and compute euclidian distance
        for (ca, la), (cb, lb) in itertools.combinations(
            zip(self.atomic_coords, self.atomic_symbols), 2
        ):
            pair_label = f"{la}: {ca}\n{lb}: {cb}"
            R = self.euclidian_distance(ca, cb)

            # Compare distance to internal thresholds
            if R < self.THRESHOLD_NUCLEAR_SINGULARITY_ERROR:
                error_pairs.append(pair_label)
            elif R < self.THRESHOLD_NUCLEAR_SINGULARITY_WARNING:
                warning_pairs.append(pair_label)

        # Print warnings and raise exception if necessary
        if warning_pairs:
            msg = self.WARNING_MESSAGE_NUCLEAR_SINGULARITY(
                "\n\n".join(warning_pairs)
            )
            print(msg)

        if error_pairs:
            msg = self.ERROR_MESSAGE_NUCLEAR_SINGULARITY(
                "\n\n".join(error_pairs)
            )
            raise RuntimeError(msg)

    def check_for_invalid_electronic_configuration(self):
        """Check that the number of electrons and spin multiplicity are compatible.
        Also check for restricted open-shell calculation."""
        restricted = self.user_dict["WaveFunction"]["restricted"]
        Z = sum([PeriodicTable[atom.lower()].Z for atom in self.atomic_symbols])
        n_electrons = Z - self.charge
        n_unpaired = self.mult - 1

        # Helper function
        parity = lambda n: "even" if n % 2 == 0 else "odd"

        # Check for impossible charge
        if self.charge > Z:
            raise RuntimeError(
                self.ERROR_UNPHYSICAL_CHARGE(
                    f"The specified charge ({self.charge}) cannot be larger than the nuclear charge ({Z})"
                )
            )

        # Check for unphysical multiplicity
        elif n_unpaired > n_electrons:
            raise RuntimeError(
                self.ERROR_UNPHYSICAL_MULTIPLICITY(
                    f"The specified multiplicity requires more unpaired electrons ({self.mult - 1}) than are available ({n_electrons}))."
                )
            )

        # Check for invalid spin multiplicity
        elif parity(n_electrons) == parity(self.mult):
            raise RuntimeError(
                self.ERROR_INCOMPATIBLE_MULTIPLICITY(
                    f"The specified multiplicity ({parity(self.mult)}) is not compatible with the number of electrons ({parity(n_electrons)})"
                )
            )

        # Check for restricted open-shell
        elif restricted and n_unpaired > 0:
            raise RuntimeError(self.ERROR_RESTRICTED_OPEN_SHELL)

    def translate_com_to_origin(self):
        """Translate center of mass to the origin (in-place)."""
        masses = [
            PeriodicTable[label.lower()].mass for label in self.atomic_symbols
        ]
        M = sum(masses)

        # Compute center of mass
        com = []
        for dim in range(3):
            component = 0.0
            for atom in range(self.n_atoms):
                component += self.atomic_coords[atom][dim] * masses[atom]
            com.append(component / M - self.origin[dim])

        # Translate coordinates
        for dim in range(3):
            for atom in range(self.n_atoms):
                self.atomic_coords[atom][dim] -= com[dim]

    @staticmethod
    def euclidian_distance(a, b):
        """Helper function for the nuclear singularies validation.
        Computes the euclidian distance between two vectors, a and b."""
        squared_deviations = [(a[i] - b[i]) ** 2 for i in range(3)]
        return math.sqrt(sum(squared_deviations))

    def ang2bohr_array(self, coords):
        """Convert List[List[float]] from angstrom to bohr."""
        return [
            [c * self.pc["angstrom2bohrs"] for c in element]
            for element in coords
        ]

    def ang2bohr_vector(self, vec):
        """Convert List[float] from angstrom to bohr"""
        return [el * self.pc["angstrom2bohrs"] for el in vec]
