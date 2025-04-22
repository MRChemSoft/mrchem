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

from math import sqrt
from pathlib import Path

from .CUBEparser import parse_files

# yapf: disable
SHORTHAND_FUNCTIONALS = [
    'svwn3',
    'svwn5',
    'pbe',
    'pbe0',
    'bpw91',
    'bp86',
    'b3p86',
    'b3p86-g',
    'blyp',
    'b3lyp',
    'b3lyp-g',
    'olyp',
    'kt1',
    'kt2',
    'kt3'
]
# yapf: enable
"""List of recognized shorthands for functionals"""


def write_scf_fock(user_dict, wf_dict, origin):
    fock_dict = {}

    # ZORA
    if user_dict["WaveFunction"]["relativity"].lower() == "zora":
        fock_dict["zora_operator"] = {
            "include_nuclear": user_dict["ZORA"]["include_nuclear"],
            "include_coulomb": user_dict["ZORA"]["include_coulomb"],
            "include_xc": user_dict["ZORA"]["include_xc"],
            "isAZORA": False,
            "azora_potential_path": user_dict["ZORA"]["azora_potential_path"]
        }
    # AZORA
    if user_dict["WaveFunction"]["relativity"].lower() == "azora":
        fock_dict["zora_operator"] = {
            "include_nuclear": False,
            "include_coulomb": False,
            "include_xc": False,
            "isAZORA": True
        }
        if user_dict["ZORA"]["azora_potential_path"].lower() != "none":
            fock_dict["zora_operator"]["azora_potential_path"] = user_dict["ZORA"]["azora_potential_path"]

    # Kinetic
    fock_dict["kinetic_operator"] = {"derivative": user_dict["Derivatives"]["kinetic"]}

    # Nuclear
    fock_dict["nuclear_operator"] = {
        "proj_prec": user_dict["Precisions"]["nuclear_prec"],
        "smooth_prec": user_dict["Precisions"]["nuclear_prec"],
        "nuclear_model": user_dict["WaveFunction"]["nuclear_model"],
        "shared_memory": user_dict["MPI"]["share_nuclear_potential"],
    }

    # Reaction
    if user_dict["WaveFunction"]["environment"].lower() != "none":
        fock_dict["reaction_operator"] = _reaction_operator_handler(user_dict)

    # Coulomb
    if wf_dict["method_type"] in ["hartree", "hf", "dft"]:
        fock_dict["coulomb_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "shared_memory": user_dict["MPI"]["share_coulomb_potential"],
        }

    # Exchange
    if wf_dict["method_type"] in ["hf", "dft"]:
        fock_dict["exchange_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "exchange_prec": user_dict["Precisions"]["exchange_prec"],
        }

    # Exchange-Correlation
    if wf_dict["method_type"] in ["dft"]:
        func_dict = []
        for line in wf_dict["dft_funcs"].split("\n"):
            sp = line.split()
            if len(sp) > 0:
                func = sp[0].lower()
                coef = [1.0]
                if len(sp) > 1:
                    coef = list(map(float, sp[1:]))
                func_dict.append({"name": func, "coef": coef[0]})
        fock_dict["xc_operator"] = {
            "shared_memory": user_dict["MPI"]["share_xc_potential"],
            "xc_functional": {
                "spin": user_dict["DFT"]["spin"],
                "cutoff": user_dict["DFT"]["density_cutoff"],
                "functionals": func_dict,
            },
        }

    # External electric field
    if len(user_dict["ExternalFields"]["electric_field"]) > 0:
        fock_dict["external_operator"] = {
            "electric_field": user_dict["ExternalFields"]["electric_field"],
            "r_O": origin,
        }

    return fock_dict


def _reaction_operator_handler(user_dict, rsp=False):
    # convert density_type from string to integer
    if user_dict["PCM"]["SCRF"]["density_type"] == "total":
        density_type = 0
    elif user_dict["PCM"]["SCRF"]["density_type"] == "electronic":
        density_type = 1
    else:
        density_type = 2

    # reaction field operator settings common to all continuum models
    reo_dict = {
        "solver_type": "Generalized_Poisson",
        "poisson_prec": user_dict["world_prec"],
        "kain": user_dict["PCM"]["SCRF"]["kain"],
        "max_iter": user_dict["PCM"]["SCRF"]["max_iter"],
        "dynamic_thrs": user_dict["PCM"]["SCRF"]["dynamic_thrs"],
        # if doing a response calculation, then density_type is set to 1 (electronic only)
        "density_type": 1 if rsp else density_type,
        "epsilon_in": user_dict["PCM"]["Solvent"]["Permittivity"]["epsilon_in"],
        "epsilon_static": user_dict["PCM"]["Solvent"]["Permittivity"]["epsilon_out"]["static"],
        "epsilon_dynamic": user_dict["PCM"]["Solvent"]["Permittivity"]["epsilon_out"]["dynamic"],
        "nonequilibrium": user_dict["PCM"]["Solvent"]["Permittivity"]["epsilon_out"][
            "nonequilibrium"
        ],
        "formulation": user_dict["PCM"]["Solvent"]["Permittivity"]["formulation"],
        "kappa_out": 0.0,
        "ion_radius": user_dict["PCM"]["Solvent"]["DebyeHuckelScreening"]["ion_radius"],
        "ion_width": user_dict["PCM"]["Solvent"]["DebyeHuckelScreening"]["ion_width"],
        "DHS-formulation": user_dict["PCM"]["Solvent"]["DebyeHuckelScreening"]["formulation"],
    }

    # ionic solvent continuum model
    ionic_model = user_dict["WaveFunction"]["environment"].lower().split("_")[-1]
    if ionic_model in ("pb", "lpb"):
        permittivity = user_dict["PCM"]["Solvent"]["Permittivity"]["epsilon_out"]["static"]
        ionic_strength = user_dict["PCM"]["Solvent"]["DebyeHuckelScreening"]["ion_strength"]
        kappa_out = compute_kappa(user_dict["Constants"], permittivity, ionic_strength)
        reo_dict |= {
            "kappa_out": kappa_out,
            "solver_type": "Poisson-Boltzmann"
            if ionic_model == "pb"
            else "Linearized_Poisson-Boltzmann",
        }

    return reo_dict


def write_scf_guess(user_dict, wf_dict):
    guess_str = user_dict["SCF"]["guess_type"].lower()
    guess_type = guess_str.split("_")[0]
    nmix = user_dict["SCF"]["initial_mixing_steps"]
    alpha_mix = user_dict["SCF"]["initial_mixing_step_size"]
    nao_directory = user_dict["SCF"]["nao_directory"]
    zeta = 0

    scf_dict = user_dict["SCF"]

    guess_prec = scf_dict["guess_prec"]

    if guess_type == "chk":
        # At least one orbital must be present in the checkpoint folder
        chk_Phi = Path(f"{scf_dict['path_checkpoint']}/phi_scf_idx_0.meta")
        if not chk_Phi.is_file():
            print(
                f"No checkpoint guess found in {scf_dict['path_checkpoint']}, falling back to 'sad_gto' initial guess"
            )
            guess_type = "sad_gto"
        else:
            # adjust guess precision if checkpoint files are present
            guess_prec = user_dict["world_prec"]

    if guess_type in ["core", "sad"]:
        zeta_str = guess_str.split("_")[1]
        if zeta_str == "sz":
            zeta = 1
        elif zeta_str == "dz":
            zeta = 2
        elif zeta_str == "tz":
            zeta = 3
        elif zeta_str == "qz":
            zeta = 4
        elif zeta_str == "gto":
            guess_type = guess_str
        else:
            print("Invalid zeta:" + zeta_str)

    file_dict = user_dict["Files"]

    if guess_type == "cube":
        found = parse_files(user_dict)
        if not found:
            print(
                f"No CUBE guess found in any of the 'initial_guess' sub-folders, falling back to 'sad_gto' initial guess"
            )
            guess_type = "sad_gto"

    vector_dir = file_dict["cube_vectors"]
    guess_dict = {
        "zeta": zeta,
        "prec": guess_prec,
        "type": guess_type,
        "method": wf_dict["method_name"],
        "relativity": wf_dict["relativity_name"],
        "environment": wf_dict["environment_name"],
        "external_field": wf_dict["external_name"],
        "screen": scf_dict["guess_screen"],
        "localize": scf_dict["localize"],
        "rotate": scf_dict["guess_rotate"],
        "restricted": user_dict["WaveFunction"]["restricted"],
        "file_chk": f"{scf_dict['path_checkpoint']}/phi_scf",
        "file_basis": file_dict["guess_basis"],
        "file_gto_p": file_dict["guess_gto_p"],
        "file_gto_a": file_dict["guess_gto_a"],
        "file_gto_b": file_dict["guess_gto_b"],
        "file_phi_p": file_dict["guess_phi_p"] + "_scf",
        "file_phi_a": file_dict["guess_phi_a"] + "_scf",
        "file_phi_b": file_dict["guess_phi_b"] + "_scf",
        "file_CUBE_p": f"{vector_dir}CUBE_p_vector.json",
        "file_CUBE_a": f"{vector_dir}CUBE_a_vector.json",
        "file_CUBE_b": f"{vector_dir}CUBE_b_vector.json",
        "initial_mixing_steps": nmix,
        "initial_mixing_step_size": alpha_mix,
    }
    if nao_directory != "none":
        guess_dict["nao_directory"] = nao_directory
    return guess_dict


def write_scf_solver(user_dict, wf_dict):
    # SCF precisions and thresholds
    start_prec = user_dict["SCF"]["start_prec"]
    final_prec = user_dict["SCF"]["final_prec"]
    if final_prec < 0.0:
        final_prec = user_dict["world_prec"]
    if start_prec < 0.0:
        start_prec = final_prec

    scf_dict = user_dict["SCF"]
    solver_dict = {
        "method": wf_dict["method_name"],
        "relativity": wf_dict["relativity_name"],
        "environment": wf_dict["environment_name"],
        "external_field": wf_dict["external_name"],
        "kain": scf_dict["kain"],
        "max_iter": scf_dict["max_iter"],
        "rotation": scf_dict["rotation"],
        "localize": scf_dict["localize"],
        "file_chk": scf_dict["path_checkpoint"] + "/phi_scf",
        "checkpoint": scf_dict["write_checkpoint"],
        "start_prec": start_prec,
        "final_prec": final_prec,
        "energy_thrs": scf_dict["energy_thrs"],
        "orbital_thrs": scf_dict["orbital_thrs"],
        "helmholtz_prec": user_dict["Precisions"]["helmholtz_prec"],
    }

    return solver_dict


def write_scf_properties(user_dict, origin):
    prop_dict = {}
    if user_dict["Properties"]["dipole_moment"]:
        prop_dict["dipole_moment"] = {} 
        prop_dict["dipole_moment"]["dip-1"] = {
            "operator": "h_e_dip",
            "precision": user_dict["world_prec"],
            "r_O": origin,
        }
    if user_dict["Properties"]["quadrupole_moment"]:
        prop_dict["quadrupole_moment"] = {}
        prop_dict["quadrupole_moment"]["quad-1"] = {
            "operator": "h_e_quad",
            "precision": user_dict["world_prec"],
            "r_O": origin,
        }
    if user_dict["Properties"]["geometric_derivative"]:
        prop_dict["geometric_derivative"] = {}
        prop_dict["geometric_derivative"]["geom-1"] = {
            "operator": "h_nuc_grad",
            "precision": user_dict["world_prec"],
            "smoothing": user_dict["Precisions"]["nuclear_prec"],
            "method": user_dict["Forces"]["method"],
            "surface_integral_precision": user_dict["Forces"]["surface_integral_precision"],
            "radius_factor": user_dict["Forces"]["radius_factor"],
        }
    if user_dict["Properties"]["hirshfeld_charges"]:
        prop_dict["hirshfeld_charges"] = {}
        prop_dict["hirshfeld_charges"]["hirshfeld-1"] = {
            'precision': user_dict["world_prec"]
        }

    return prop_dict


def write_scf_plot(user_dict):
    plot_dict = {}
    if user_dict["Properties"]["plot_density"] or len(
        user_dict["Properties"]["plot_orbitals"]
    ):
        plot_dict["orbitals"] = user_dict["Properties"]["plot_orbitals"]
        plot_dict["density"] = user_dict["Properties"]["plot_density"]
        plot_dict["plotter"] = user_dict["Plotter"]
        if user_dict["world_unit"] == "angstrom":
            plot_dict["plotter"] = {
                k: [
                    user_dict["Constants"]["angstrom2bohrs"] * r
                    for r in plot_dict["plotter"][k]
                ]
                for k in plot_dict["plotter"].keys()
            }
    return plot_dict


def write_rsp_calc(omega, user_dict, origin):
    wf_dict = parse_wf_method(user_dict)
    if not wf_dict["relativity_name"] in ["None", "Off"]:
        raise RuntimeError(
            "Linear response not available: " + wf_dict["relativity_name"]
        )

    rsp_dict = user_dict["Response"]
    file_dict = user_dict["Files"]

    rsp_calc = {}
    rsp_calc["frequency"] = omega
    rsp_calc["dynamic"] = omega > 1.0e-12
    rsp_calc["fock_operator"] = write_rsp_fock(user_dict, wf_dict)
    rsp_calc["unperturbed"] = {
        "precision": user_dict["world_prec"],
        "localize": rsp_dict["localize"],
        "fock_operator": write_scf_fock(user_dict, wf_dict, origin),
    }

    guess_str = rsp_dict["guess_type"].lower()
    user_guess_type = guess_str.split("_")[0]
    user_guess_prec = rsp_dict["guess_prec"]

    vector_dir = file_dict["cube_vectors"]

    rsp_calc["components"] = []
    for dir in [0, 1, 2]:
        rsp_comp = {}

        program_guess_type = user_guess_type
        program_guess_prec = user_guess_prec

        # check that initial guess files exist
        if user_guess_type == "chk":
            chk_X = Path(f"{rsp_dict['path_checkpoint']}/X_rsp_{dir:d}")
            chk_Y = Path(f"{rsp_dict['path_checkpoint']}/Y_rsp_{dir:d}")
            if not (chk_X.is_file() and chk_Y.is_file()):
                print(
                    f"No checkpoint guess found in {rsp_dict['path_checkpoint']} for direction {dir:d}, falling back to zero initial guess"
                )
                program_guess_type = "none"
            else:
                # adjust guess precision if checkpoint files are present
                program_guess_prec = user_dict["world_prec"]
        elif user_guess_type == "cube":
            found = parse_files(user_dict, dir)
            if not found:
                print(
                    f"No CUBE guess found in any of the 'initial_guess' sub-folders for direction {dir:d}, falling back to zero initial guess"
                )
                program_guess_type = "none"
        else:
            # do no checks on other types of guess
            pass

        rsp_comp["initial_guess"] = {
            "prec": program_guess_prec,
            "type": program_guess_type,
            "file_chk_x": f"{rsp_dict['path_checkpoint']}/X_rsp_{dir:d}",
            "file_chk_y": f"{rsp_dict['path_checkpoint']}/Y_rsp_{dir:d}",
            "file_x_p": f"{file_dict['guess_x_p']}_rsp_{dir:d}",
            "file_x_a": f"{file_dict['guess_x_a']}_rsp_{dir:d}",
            "file_x_b": f"{file_dict['guess_x_b']}_rsp_{dir:d}",
            "file_y_p": f"{file_dict['guess_y_p']}_rsp_{dir:d}",
            "file_y_a": f"{file_dict['guess_y_a']}_rsp_{dir:d}",
            "file_y_b": f"{file_dict['guess_y_b']}_rsp_{dir:d}",
            "file_CUBE_x_p": f"{vector_dir}CUBE_x_p_{dir:d}_vector.json",
            "file_CUBE_x_a": f"{vector_dir}CUBE_x_a_{dir:d}_vector.json",
            "file_CUBE_x_b": f"{vector_dir}CUBE_x_b_{dir:d}_vector.json",
            "file_CUBE_y_p": f"{vector_dir}CUBE_y_p_{dir:d}_vector.json",
            "file_CUBE_y_a": f"{vector_dir}CUBE_y_a_{dir:d}_vector.json",
            "file_CUBE_y_b": f"{vector_dir}CUBE_y_b_{dir:d}_vector.json",
        }
        if rsp_dict["write_orbitals"]:
            path_orbitals = rsp_dict["path_orbitals"]
            rsp_comp["write_orbitals"] = {
                "file_x_p": f"{path_orbitals}/X_p_rsp_{dir:d}",
                "file_x_a": f"{path_orbitals}/X_a_rsp_{dir:d}",
                "file_x_b": f"{path_orbitals}/X_b_rsp_{dir:d}",
                "file_y_p": f"{path_orbitals}/Y_p_rsp_{dir:d}",
                "file_y_a": f"{path_orbitals}/Y_a_rsp_{dir:d}",
                "file_y_b": f"{path_orbitals}/Y_b_rsp_{dir:d}",
            }
        if rsp_dict["run"][dir]:
            rsp_comp["rsp_solver"] = write_rsp_solver(user_dict, wf_dict, dir)
        rsp_calc["components"].append(rsp_comp)
    return rsp_calc


def write_rsp_fock(user_dict, wf_dict):
    fock_dict = {}

    # Coulomb
    if wf_dict["method_type"] in ["hartree", "hf", "dft"]:
        fock_dict["coulomb_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "shared_memory": user_dict["MPI"]["share_coulomb_potential"],
        }

    # Exchange
    if wf_dict["method_type"] in ["hf", "dft"]:
        fock_dict["exchange_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "exchange_prec": user_dict["Precisions"]["exchange_prec"],
        }

    # Exchange-Correlation
    if wf_dict["method_type"] in ["dft"]:
        func_dict = []
        for line in wf_dict["dft_funcs"].split("\n"):
            sp = line.split()
            if len(sp) > 0:
                func = sp[0].lower()
                coef = [1.0]
                if len(sp) > 1:
                    coef = list(map(float, sp[1:]))
                func_dict.append({"name": func, "coef": coef[0]})
        fock_dict["xc_operator"] = {
            "shared_memory": user_dict["MPI"]["share_xc_potential"],
            "xc_functional": {
                "spin": user_dict["DFT"]["spin"],
                "cutoff": user_dict["DFT"]["density_cutoff"],
                "functionals": func_dict,
            },
        }

    # Reaction
    if user_dict["WaveFunction"]["environment"].lower() != "none":
        fock_dict["reaction_operator"] = _reaction_operator_handler(user_dict, rsp=True)

    return fock_dict


def write_rsp_solver(user_dict, wf_dict, d):
    # Response precisions and thresholds
    start_prec = user_dict["Response"]["start_prec"]
    final_prec = user_dict["Response"]["final_prec"]
    if final_prec < 0.0:
        final_prec = user_dict["world_prec"]
    if start_prec < 0.0:
        start_prec = final_prec

    rsp_dict = user_dict["Response"]
    solver_dict = {
        "method": wf_dict["method_name"],
        "kain": rsp_dict["kain"],
        "max_iter": rsp_dict["max_iter"],
        "file_chk_x": rsp_dict["path_checkpoint"] + "/X_rsp_" + str(d),
        "file_chk_y": rsp_dict["path_checkpoint"] + "/Y_rsp_" + str(d),
        "checkpoint": rsp_dict["write_checkpoint"],
        "start_prec": start_prec,
        "final_prec": final_prec,
        "orbital_thrs": user_dict["Response"]["orbital_thrs"],
        "property_thrs": user_dict["Response"]["property_thrs"],
        "helmholtz_prec": user_dict["Precisions"]["helmholtz_prec"],
        "orth_prec": 1.0e-14,
    }
    return solver_dict


def parse_psppar(filename):
    """
    Parses a PSPPAR pseudopotential file and returns the data as a dictionary.
    
    Args:
        filename (str): Path to the PSPPAR file.

    Returns:
        dict: Parsed pseudopotential data.
    """
    psppar_data = dict()
    psppar_data["local"] = dict()
    psppar_data["nonlocal"] = dict()
    
    with open(filename, 'r') as file:
        # Skip the first line
        next(file)
        
        # Read zion and zeff
        line = file.readline().strip()
        zion = float(line.split()[0])
        zeff = float(line.split()[1])
        psppar_data['zion'] = zion
        psppar_data['zeff'] = zeff
        
        # Skip the next line
        next(file)
        
        # Read rloc, nloc, and c coefficients
        line = file.readline().strip()
        tokens = line.split()
        rloc = float(tokens[0])
        nloc = int(tokens[1])
        c_coefficients = list(map(float, tokens[2:2+nloc]))
        
        alpha_pp = 1.0 / (sqrt(2.0) * rloc)
        psppar_data['local']['rloc'] = rloc
        psppar_data['local']['alpha_pp'] = alpha_pp
        psppar_data['local']['nloc'] = nloc
        psppar_data['local']['c'] = c_coefficients
        
        # Read nsep
        line = file.readline().strip()
        nsep = int(line.split()[0])
        if nsep < 0 or nsep > 4:
            raise ValueError("Error: nsep must be between 0 and 4.")
        psppar_data['nonlocal']['nsep'] = nsep
        
        # Initialize containers for projectors
        rl = []
        h = []
        dim_h = []

        
        for l in range(nsep):
            # Read projector radii and dimension
            line = file.readline().strip()
            tokens = line.split()
            rl_l = float(tokens[0])
            dim_h_l = int(tokens[1])
            
            rl.append(rl_l)
            dim_h.append(dim_h_l)
            
            # Initialize h matrix
            h_matrix = []
            for _ in range(dim_h_l):
                h_matrix.append([0.0] * dim_h_l)
            # convert h_matrix to list of lists

            for i in range(dim_h_l):
                h_matrix[0][i] = float(line.split()[2 + i])
                h_matrix[i][0] = h_matrix[0][i]
            
            # Fill in the upper triangle of the h matrix
            for i in range(1, dim_h_l):
                line = file.readline().strip()
                values = list(map(float, line.split()))
                for j, value in enumerate(values):
                    h_matrix[i][i + j] = value
                    h_matrix[i + j][i] = value
            
            # h_matrix = h_matrix.tolist()
            h.append(h_matrix)
            
            # Skip irrelevant SOC lines if l > 0
            if l > 0:
                for _ in range(dim_h_l):
                    file.readline()
            
        # check if file is at the end
        line = file.readline()
        if line: # try to read rcore and qcore
            words = line.split()
            if len(words) >= 2:
                rcore = float(words[0])
                qcore = float(words[1])
                psppar_data['nlcc'] = dict()
                psppar_data['nlcc']['rcore'] = rcore
                psppar_data['nlcc']['qcore'] = qcore
            
        
        psppar_data['nonlocal']['rl'] = rl
        psppar_data['nonlocal']['dim_h'] = dim_h
        psppar_data['nonlocal']['h'] = h
            
    return psppar_data


def write_pseudo_potential(user_dict, mol_dict):
    import json

    pp_dict_out = {
        "use_pp": False,
        "pp_prec": user_dict["Pseudopotential"]["pp_prec"],
    }

    # if "Pseudopotential" not in user_dict:
    #     return []

    pp_dict = json.loads(user_dict['Pseudopotential']['pp_files'])

    # loop over keys in pp_dict and convert them to lower case
    pp_dict = {k.lower(): v for k, v in pp_dict.items()}
    all_elements = {k.lower() for k in user_dict["Elements"]}

    int_keys = []
    # check if all keys in pp_dict are valid elements and collect the integer keys
    for key in pp_dict.keys():
        if key not in all_elements:
            try:
                int_keys.append(int(key))
            except ValueError:
                raise RuntimeError(
                    f"Invalid key in pseudopotential dictionary (neither a numeric value or an element): {key}"
                )

    nat = len(mol_dict["coords"])
    atomNames = []
    atomNamesSet = set()
    for atom in mol_dict["coords"]:
        atomNames.append(atom["atom"])
        atomNamesSet.add(atom["atom"])
    
    pps = []
    for i in range(nat):
        key = ''
        if i in int_keys:
            key = str(i)
        # check if the key is an element
        elif atomNames[i] in pp_dict.keys():
            key = atomNames[i]
        # if none of the above, all electron calculation for atom i
        else:
            pp = {
                "use_pp": False,
            }
            pps.append(pp)
        
        if key != "":
            if pp_dict[key].lower() == "none" or pp_dict[key].lower() == "":
                pp = {
                    "use_pp": False,
                }
            else:
                pp = {
                    "use_pp": True,
                    "file": pp_dict[key],
                }
                pp_dict_out["use_pp"] = True
                # temp_string = json.dumps(temp_pp, indent=4, cls=NumpyEncoder)
                # print("Parsed pseudopotential for element:", key, "\n \n")
                # print(temp_string)
                # print("\n \n \n")
                pp['pseuodopotential'] = parse_psppar(pp["file"])
            pps.append(pp)
    # print(pps)
    pp_dict_out["pp_list"] = pps
    return pp_dict_out


def parse_wf_method(user_dict):
    method_name = ""
    restricted = user_dict["WaveFunction"]["restricted"]
    method_type = user_dict["WaveFunction"]["method"].lower()
    dft_funcs = user_dict["DFT"]["functionals"].lower()
    if method_type in ["core"]:
        method_name = "Core Hamiltonian"
    elif method_type in ["hartree"]:
        method_name = "Hartree"
    elif method_type in ["hf", "hartree-fock", "hartreefock"]:
        method_name = "Hartree-Fock"
        method_type = "hf"
    elif method_type in ["dft"]:
        method_name = "DFT"
    elif method_type in ["lda"]:
        method_name = "DFT (SVWN5)"
        dft_funcs = "svwn5"
        method_type = "dft"
    elif method_type in SHORTHAND_FUNCTIONALS:
        method_name = "DFT (" + method_type.upper() + ")"
        dft_funcs = method_type
        method_type = "dft"
    else:
        raise RuntimeError(
            f"Invalid wavefunction method {user_dict['WaveFunction']['method']}"
        )

    # Determine relativity name label for print outs to the output file
    relativity_name = "None"
    wf_dict = dict()
    if user_dict["WaveFunction"]["relativity"].lower() in ["none"]:
        user_dict["WaveFunction"]["relativity"] = "off"
        user_dict["ZORA"]["include_nuclear"] = False
        user_dict["ZORA"]["include_coulomb"] = False
        user_dict["ZORA"]["include_xc"] = False

    if user_dict["WaveFunction"]["relativity"].lower() in ["azora"]:
        relativity_name = "AZORA"
        user_dict["WaveFunction"]["relativity"] = "azora"
        user_dict["ZORA"]["include_nuclear"] = False
        user_dict["ZORA"]["include_coulomb"] = False
        user_dict["ZORA"]["include_xc"] = False
        user_dict["ZORA"]["isAZORA"] = True

    if user_dict["WaveFunction"]["relativity"].lower() in ["nzora"]:
        user_dict["WaveFunction"]["relativity"] = "zora"
        user_dict["ZORA"]["include_nuclear"] = True
        user_dict["ZORA"]["include_coulomb"] = False
        user_dict["ZORA"]["include_xc"] = False

    if user_dict["WaveFunction"]["relativity"].lower() in ["zora"]:
        components = [
            user_dict["ZORA"]["include_nuclear"],
            user_dict["ZORA"]["include_coulomb"],
            user_dict["ZORA"]["include_xc"],
        ]
        names = ["V_nuc", "J", "V_xc"]

        if any(components):
            zora_terms = " + ".join(
                [name for name, comp in zip(names, components) if comp]
            )
            relativity_name = "ZORA (" + zora_terms + ")"
        else:
            raise RuntimeError("ZORA selected, but no ZORA potential included")

        # if user_dict["ZORA"]["include_xc"] and not restricted:
        #     raise RuntimeError(
        #         "ZORA (V_xc) not available for unrestricted wavefunctions"
        #     )

    # Determine environment name label for print outs to the output file
    environment_name = "None"
    if user_dict["WaveFunction"]["environment"].lower() == "pcm":
        environment_name = "PCM"

    # Determine external name label for print outs to the output file
    ext_dict = user_dict["ExternalFields"]
    has_external_fields = len(ext_dict["electric_field"]) > 0

    external_name = "None"
    if has_external_fields:
        # If no external fields, then the list will be empty
        # Need to catch the exception and store placeholders
        try:
            x, y, z = ext_dict["electric_field"]
        except ValueError:
            x, y, z = None, None, None  # Useless placeholders

        # Labels to aggregate
        external_name = f"Electric field ({x}, {y}, {z})"

    wf_dict["relativity_name"] = relativity_name
    wf_dict["environment_name"] = environment_name
    wf_dict["external_name"] = external_name
    wf_dict["method_name"] = method_name
    wf_dict["method_type"] = method_type
    wf_dict["dft_funcs"] = dft_funcs
    return wf_dict


def compute_kappa(constants, eps, I):
    kb = constants["boltzmann_constant"]
    e = constants["elementary_charge"]
    e_0 = constants["e0"]
    N_a = constants["N_a"]
    m2au = constants["meter2bohr"]
    T = 298.15

    numerator = e_0 * eps * kb * T
    denominator = 2.0 * (e**2) * N_a * 1000.0 * I

    debye_length = sqrt(numerator / denominator) * m2au

    return 1.0 / debye_length
