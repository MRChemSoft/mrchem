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

from .helpers import *
from .CUBEparser import *
from .periodictable import PeriodicTable as PT, PeriodicTableByZ as PT_Z
from .validators import MoleculeValidator


def translate_input(user_dict):
    # get the origin in the desired units of measure
    origin = user_dict["world_origin"]
    if user_dict["world_unit"] == "angstrom":
        origin = [ANGSTROM_2_BOHR * r for r in origin]

    # prepare bits and pieces
    mol_dict = MoleculeValidator(user_dict, origin).write_mol()
    mpi_dict = write_mpi(user_dict)
    mra_dict = write_mra(user_dict, mol_dict)
    scf_dict = write_scf_calculation(user_dict, origin)
    rsp_dict = write_rsp_calculations(user_dict, mol_dict, origin)

    # piece everything together
    program_dict = {
        "schema_name": "mrchem_input",
        "schema_version": 1,
        "mpi": mpi_dict,
        "mra": mra_dict,
        "printer": user_dict["Printer"],
        "molecule": mol_dict,
        "scf_calculation": scf_dict,
        "rsp_calculations": rsp_dict,
    }
    return program_dict


def write_mpi(user_dict):
    mpi_dict = {
        "numerically_exact": user_dict["MPI"]["numerically_exact"],
        "shared_memory_size": user_dict["MPI"]["shared_memory_size"],
        "bank_size": user_dict["MPI"]["bank_size"],
    }
    return mpi_dict


def write_mra(user_dict, mol_dict):
    order = user_dict["Basis"]["order"]
    if order < 0:
        # Set polynomial order based on world_prec
        prec = user_dict["world_prec"]
        order = int(math.ceil(-1.5 * math.log10(prec)))

    min_scale = -(user_dict["world_size"] - 1)
    if min_scale > 1:
        # Compute auto box
        max_coord = 0.0  # single (coord + Z) with largest abs value
        for coord in mol_dict["coords"]:
            Z_i = max(5.0, float(PT[coord["atom"]].Z))
            max_coord = max(max_coord, abs(max(coord["xyz"], key=abs)) + 2.0 * Z_i)

        min_scale = 0
        while 2.0 ** (-min_scale) < max_coord:
            min_scale = min_scale - 1
    else:
        # Approximately scale world angstrom -> bohr
        if user_dict["world_unit"] == "angstrom":
            min_scale -= 1

    max_scale = 20
    if (max_scale - min_scale) > 30:
        max_scale = 30 + min_scale

    mra_dict = {
        "basis_type": user_dict["Basis"]["type"].lower(),
        "basis_order": order,
        "boxes": [2, 2, 2],
        "corner": [-1, -1, -1],
        "min_scale": min_scale,
        "max_scale": max_scale,
    }
    return mra_dict


############################################################
#                         MOLECULE                         #
############################################################


def write_molecule(user_dict, origin):
    # Translate into program syntax
    coords_raw = user_dict["Molecule"]["coords"]
    coords_dict = []
    for line in coords_raw.split("\n"):
        sp = line.split()
        if len(sp) > 0:

            try:
                int(sp[0])
            except:
                atom = sp[0].lower()
            else:
                atom = PT_Z[int(sp[0])].symbol.lower()

            xyz = list(map(float, sp[1:]))
            if len(xyz) != 3:
                raise RuntimeError(f"Invalid coordinate: {atom.upper()} {str(xyz)}")
            coords_dict.append({"atom": atom, "xyz": xyz})

    # Convert angstrom -> bohr
    if user_dict["world_unit"] == "angstrom":
        for coord in coords_dict:
            coord["xyz"] = [ANGSTROM_2_BOHR * r for r in coord["xyz"]]

    # Check for singularity
    for a in range(len(coords_dict)):
        for b in range(a + 1, len(coords_dict)):
            A = coords_dict[a]
            B = coords_dict[b]
            R = math.sqrt(sum([(a - b) ** 2 for a, b in zip(A["xyz"], B["xyz"])]))
            if R < 1.0e-6:
                msg = f"ABORT: Atoms are too close\n {A['atom']}: {str(A['xyz'])}\n {B['atom']}: {str(B['xyz'])}"
                raise RuntimeError(msg)
            elif R < 1.0e-3:
                msg = f"WARNING: Atoms are too close\n {A['atom']}: {str(A['xyz'])}\n {B['atom']}: {str(B['xyz'])}"
                print(msg)

    # Translate center of mass to origin
    if user_dict["Molecule"]["translate"]:
        # Calc center of mass
        M = 0.0
        CoM = [0.0, 0.0, 0.0]
        for coord in coords_dict:
            m = PT[coord["atom"]].mass
            M += m
            CoM = [m * r + o for r, o in zip(coord["xyz"], CoM)]
        CoM = [x / M - o for x, o in zip(CoM, origin)]

        # Translate coords
        for coord in coords_dict:
            coord["xyz"] = [r - o for r, o in zip(coord["xyz"], CoM)]

    # initialize the cavity
    cav_coords_dict = []
    if len(user_dict["Environment"]["Cavity"]["spheres"]) > 0:
        cav_coords_raw = user_dict["Environment"]["Cavity"]["spheres"]
        for line in cav_coords_raw.split('\n'):
            sp = line.split()
            if len(sp) > 0:
                center = list(map(float, sp[:-1]))
                radius = float(sp[-1])
                if len(center) != 3:
                    raise RuntimeError(f"Invalid coordinate: {center}")
            cav_coords_dict.append({
                "center": center,
                "radius": radius
            })
    else:
        for atom in coords_dict:
            radius = float(PT[atom["atom"]].radius)
            cav_coords_dict.append({
                "center": atom["xyz"],
                "radius": radius
            })

    if user_dict["world_unit"] == "angstrom" :
        for coord in cav_coords_dict:
            coord["center"] = [ANGSTROM_2_BOHR * r for r in coord["center"]]
            coord["radius"] *= ANGSTROM_2_BOHR
        cavity_width = user_dict["Environment"]["Cavity"]["cavity_width"]*ANGSTROM_2_BOHR
    else:
        cavity_width = user_dict["Environment"]["Cavity"]["cavity_width"]

    # Make sure the specified charge and multiplicity make sense
    mult = user_dict["Molecule"]["multiplicity"]
    charge = user_dict["Molecule"]["charge"]
    restricted = user_dict["WaveFunction"]["restricted"]

    # Collect number of electrons
    n_electrons = sum([PT[atom["atom"]].Z for atom in coords_dict]) - charge
    n_unpaired = mult - 1

    # Helper function
    parity = lambda n: 'even' if n % 2 == 0 else 'odd'

    # Check for unphysical multiplicity (not likely to be much of a problem, but still...)
    if n_unpaired > n_electrons:
        raise RuntimeError(f"The specified multiplicity requires more unpaired electrons ({mult - 1}) than are available ({n_electrons}))")

    # Check for restricted open-shell
    elif restricted and n_unpaired > 0: 
        raise RuntimeError("Restricted open-shell calculations are not available")

    # Check for invalid spin multiplicity
    elif parity(n_electrons) == parity(mult): 
        raise RuntimeError(f"The specified multiplicity ({parity(mult)}) is not compatible with the number of electrons ({parity(n_electrons)})")

    mol_dict = {
        "multiplicity": mult,
        "charge": charge,
        "coords": coords_dict,
        "cavity_coords": cav_coords_dict,
        "cavity_width": cavity_width
    }

    return mol_dict


############################################################
#                       SCF CALCULATION                    #
############################################################


def write_scf_calculation(user_dict, origin):
    wf_dict = parse_wf_method(user_dict)

    scf_dict = {}
    scf_dict["fock_operator"] = write_scf_fock(user_dict, wf_dict, origin)
    scf_dict["initial_guess"] = write_scf_guess(user_dict, wf_dict)

    path_orbitals = user_dict["SCF"]["path_orbitals"]
    if user_dict["SCF"]["write_orbitals"]:
        scf_dict["write_orbitals"] = {
            "file_phi_p": path_orbitals + "/phi_p_scf",
            "file_phi_a": path_orbitals + "/phi_a_scf",
            "file_phi_b": path_orbitals + "/phi_b_scf",
        }

    if user_dict["SCF"]["run"]:
        scf_dict["scf_solver"] = write_scf_solver(user_dict, wf_dict)

    prop_dict = write_scf_properties(user_dict, origin)
    if len(prop_dict) > 0:
        scf_dict["properties"] = prop_dict

    plot_dict = write_scf_plot(user_dict)
    if len(plot_dict) > 0:
        scf_dict["plots"] = plot_dict

    return scf_dict


############################################################
#                   RESPONSE CALCULATIONS                  #
############################################################


def write_rsp_calculations(user_dict, mol_dict, origin):
    rsp_dict = {}
    run_pol = user_dict["Properties"]["polarizability"]
    run_mag = user_dict["Properties"]["magnetizability"]
    run_nmr = user_dict["Properties"]["nmr_shielding"]
    nuc_spec = user_dict["NMRShielding"]["nuclear_specific"]

    if run_pol:
        for omega in user_dict["Polarizability"]["frequency"]:
            freq_key = f"{omega:6f}"
            pol_key = "pol-" + freq_key
            rsp_calc = write_rsp_calc(omega, user_dict, origin)
            rsp_calc["perturbation"] = {"operator": "h_e_dip", "r_O": origin}
            rsp_calc["properties"] = {}
            rsp_calc["properties"]["polarizability"] = {}
            rsp_calc["properties"]["polarizability"][pol_key] = {
                "frequency": omega,
                "precision": user_dict["world_prec"],
                "operator": "h_e_dip",
                "r_O": origin,
            }
            rsp_key = "ext_el-" + freq_key
            rsp_dict[rsp_key] = rsp_calc

    if run_mag or (run_nmr and not nuc_spec):
        omega = 0.0  # only static magnetic response
        freq_key = f"{omega:6f}"
        mag_key = "mag-" + freq_key
        rsp_calc = write_rsp_calc(0.0, user_dict, origin)
        rsp_calc["perturbation"] = {
            "operator": "h_b_dip",
            "derivative": user_dict["Derivatives"]["h_b_dip"],
            "r_O": origin,
        }
        rsp_calc["properties"] = {}
        if run_mag:
            rsp_calc["properties"]["magnetizability"] = {}
            rsp_calc["properties"]["magnetizability"][mag_key] = {
                "frequency": omega,
                "precision": user_dict["world_prec"],
                "dia_operator": "h_bb_dia",
                "para_operator": "h_b_dip",
                "derivative": user_dict["Derivatives"]["h_b_dip"],
                "r_O": origin,
            }
        if run_nmr and not nuc_spec:
            rsp_calc["properties"]["nmr_shielding"] = {}
            nuc_vec = user_dict["NMRShielding"]["nucleus_k"]
            all_nucs = nuc_vec[0] < 0
            nuclei = mol_dict["coords"]
            for k in range(len(nuclei)):
                if all_nucs or k in nuc_vec:
                    atom_key = str(k) + nuclei[k]["atom"]
                    nmr_key = "nmr-" + atom_key
                    rsp_calc["properties"]["nmr_shielding"][nmr_key] = {
                        "precision": user_dict["world_prec"],
                        "dia_operator": "h_bm_dia",
                        "para_operator": "h_m_pso",
                        "smoothing": user_dict["world_prec"],
                        "derivative": user_dict["Derivatives"]["h_m_pso"],
                        "r_O": origin,
                        "r_K": nuclei[k]["xyz"],
                    }
        rsp_key = "ext_mag-" + freq_key
        rsp_dict[rsp_key] = rsp_calc

    if run_nmr and nuc_spec:
        nuc_vec = user_dict["NMRShielding"]["nucleus_k"]
        all_nucs = nuc_vec[0] < 0
        nuclei = mol_dict["coords"]
        for k in range(len(nuclei)):
            if (all_nucs) or (k in nuc_vec):
                atom_key = str(k) + nuclei[k]["atom"]
                nmr_key = "nmr-" + atom_key
                rsp_calc = write_rsp_calc(0.0, user_dict, origin)
                rsp_calc["perturbation"] = {
                    "operator": "h_m_pso",
                    "smoothing": user_dict["world_prec"],
                    "derivative": user_dict["Derivatives"]["h_m_pso"],
                    "r_K": nuclei[k]["xyz"],
                }
                rsp_calc["properties"] = {}
                rsp_calc["properties"]["nmr_shielding"] = {}
                rsp_calc["properties"]["nmr_shielding"][nmr_key] = {
                    "precision": user_dict["world_prec"],
                    "dia_operator": "h_mb_dia",
                    "para_operator": "h_b_dip",
                    "smoothing": user_dict["world_prec"],
                    "derivative": user_dict["Derivatives"]["h_b_dip"],
                    "r_O": origin,
                    "r_K": nuclei[k]["xyz"],
                }
                rsp_key = "nuc_mag-" + atom_key
                rsp_dict[rsp_key] = rsp_calc

    return rsp_dict
