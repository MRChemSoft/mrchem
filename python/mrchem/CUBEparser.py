
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

import os
from json import dump

from .input_parser.plumbing import pyparsing as pp

global pc


def write_cube_dict(user_dict):
    file_dict = user_dict['Files']
    world_unit = user_dict['world_unit']
    pc = user_dict['Constants']

    all_path_list = []
    all_path_list.append(sort_paths(file_dict["guess_cube_p"]))
    all_path_list.append(sort_paths(file_dict["guess_cube_a"]))
    all_path_list.append(sort_paths(file_dict["guess_cube_b"]))
    all_cube_list = []
    for path_list in all_path_list:
        cube_list = []
        if (len(path_list) != 0):
            for path in path_list:
                cube_list.append(parse_cube_file(path, world_unit))
        all_cube_list.append(cube_list)

    vector_dir = file_dict["cube_vectors"]

    for index, x_list in enumerate(all_cube_list):
        sorted_list = sorted(x_list, key=lambda d: d["ORB_IDS"])
        all_cube_list[index] = sorted_list

    if (not os.path.isdir(vector_dir)) :
        os.mkdir(vector_dir)

    with open(f"{vector_dir}CUBE_p_vector.json", "w") as fd:
        dump(all_cube_list[0], fd, indent=2)

    with open(f"{vector_dir}CUBE_a_vector.json", "w") as fd:
        dump(all_cube_list[1], fd, indent=2)

    with open(f"{vector_dir}CUBE_b_vector.json", "w") as fd:
        dump(all_cube_list[2], fd, indent=2)


def sort_paths(path):
    path_l = []
    dir_path = "/".join(path.split("/")[:-1])
    directory = os.fsencode(dir_path)
    if (os.path.isdir(dir_path)):
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if ((filename.startswith(path.split("/")[-1])) and (filename.endswith(".cube"))):
                path_l.append(dir_path+"/"+filename)
    return path_l

#TODO do a sanity check on the naming of the files

def parse_cube_file(cube_path, world_unit):

    """
    Pyparsing CUBE file
    Authors: Roberto Di Remigio & Gabriel Gerez

    Parses a CUBE file following the standard outlined in https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html.
    There are two optional values in this standard, giving a set of 4 different ways of writing a CUBE file.

    Parameters
    ----------
    cube_path: path to Cube file to be parsed

    Returns
    -------
    a dictionary with information about which path the CUBE file was parsed from, the header of the CUBE file
    and vectors detailing the individual function values for each function plotted.
    """

    # non-zero unsigned integer
    nonzero_uint_t = pp.Word("123456789", pp.nums).setParseAction(pp.pyparsing_common.convertToInteger)
    # non-zero signed integer
    nonzero_int_t = pp.Word("+-123456789", pp.nums).setParseAction(lambda t: abs(int(t[0])))
    # floating point numbers, can be in scientific notation
    float_t = pp.pyparsing_common.sci_real

    # NVAL token
    nval_t = pp.Optional(~pp.LineEnd() + nonzero_uint_t, default=1)("NVAL")

    # Cartesian coordinates
    # it could be alternatively defined as: coords = pp.Group(float_t("x") + float_t("y") + float_t("z"))
    coords = pp.Group(float_t * 3)

    # row with molecular geometry
    geom_field_t = pp.Group(nonzero_uint_t("ATOMIC_NUMBER") + float_t("CHARGE") + coords("POSITION"))

    # specification of cube axes
    def axis_spec_t(d):
        return pp.Group(nonzero_uint_t("NVOXELS") + coords("VECTOR"))(f"{d.upper()}AXIS")

    before_t = pp.Group(float_t * 3)("ORIGIN") + nval_t + axis_spec_t("X") + axis_spec_t("Y") + axis_spec_t("Z")
    # the parse action flattens the list
    after_t = pp.Optional(pp.countedArray(pp.pyparsing_common.integer))("DSET_IDS").setParseAction(lambda t: t[0] if len(t) != 0 else t)

    # The molecular geometry is a variable-length list of `geom_field_t` tokens.
    # We use a modified implementation of `pyparsing`'s `countedArray` to define
    # the token for thte whole cubefile preamble, with:

    #- `before_t`: the tokens to be matched _before_ the molecular geometry:
    #  `NATOMS`, `ORIGIN`, `NVAL` (defaults to 1 if absent), and the cube axes (`XAXIS`, `YAXIS`, `ZAXIS`)
    #- `after_t`: the tokens to be matched _after_ the molecular geometry: `DSET_IDS`, if present.

    def preamble_t(pre, post):
        expr = pp.Forward()

        def count(s, l, t):
            n = t[0]
            expr << (geom_field_t * n)("GEOM")
            return n

        natoms_t = nonzero_int_t("NATOMS")
        natoms_t.addParseAction(count, callDuringTry=True)

        return natoms_t + pre + expr + post

    cube_t = preamble_t(before_t, after_t)

    # parse the whole file and extract both all the values and the header
    with open(cube_path, "r") as cube_file:
        # we skip the first two lines as they are only comments and put the rest of the file in a single string
        cube_str = "".join(cube_file.readlines()[2:])  # The \n are already included in the string, don't need to double include them

    parsed_cube = cube_t.parseString(cube_str).asDict()
    # manually extracting voxel values
    if ("DSET_IDS" not in parsed_cube.keys()):
        parsed_cube["DSET_IDS"] = []

    cube_s = cube_str.split("\n")

    all_data_list = []

    # parse through a list of lines where the header has been removed, but the ORB_IDS remain, and append each value in a new all_data_list.
    for line in (cube_s[(4 + abs(parsed_cube["NATOMS"])):]):
        all_data_list.extend(line.split())

    if (len(parsed_cube["DSET_IDS"]) != 0):
        voxel_list = all_data_list[(len(parsed_cube["DSET_IDS"]) + 1):] # remove ORB_IDS from the all_data_list
    else:
        voxel_list = all_data_list

    parsed_cube["DATA"] = [float(value) for value in voxel_list]

    # start extracting the header values


    # get cube origin data
    N_atoms = parsed_cube["NATOMS"]
    origin = parsed_cube["ORIGIN"] if (world_unit == "bohr") else [p * pc['angstrom2bohrs'] for p in parsed_cube["ORIGIN"]]

    # Set the amount of values depending on if the DSET_IDs were present or not
    if (len(parsed_cube["DSET_IDS"]) != 0):
        N_vals = len(parsed_cube["DSET_IDS"])
    else:
        N_vals = parsed_cube["NVAL"][0]

    # files are given as [phi,rho,x,y]_[p,a,b]_[rsp,scf]_idx_#_[re,im].cube
    #TODO test that the file name makes sense
    path_ids = cube_path.split("/")[-1].split("_")[4]
    if ("-" in path_ids):
        from_to = path_ids.split("-")
        orb_ids = list(range(from_to[0], (from_to[1]+1), 1))  # we work as including the from and to, so 4-7 includes 4, 5, 6 and 7.
        if (N_vals != len(orb_ids)):
            raise ValueError("Different amount of orbitals in file and amount of orbitals specified in file name.")
    else:
        orb_ids = [float(path_ids)]


    # get voxel axis data to construct a basis for the cube space.
    N_steps = [parsed_cube["XAXIS"]["NVOXELS"], parsed_cube["YAXIS"]["NVOXELS"], parsed_cube["ZAXIS"]["NVOXELS"]]
    if (world_unit == "bohr"):
        Voxel_axes = [parsed_cube["XAXIS"]["VECTOR"], parsed_cube["YAXIS"]["VECTOR"], parsed_cube["ZAXIS"]["VECTOR"]]
    else:
        X_voxel = [p * pc['angstrom2bohrs'] for p in parsed_cube["XAXIS"]["VECTOR"]]
        Y_voxel = [p * pc['angstrom2bohrs'] for p in parsed_cube["YAXIS"]["VECTOR"]]
        Z_voxel = [p * pc['angstrom2bohrs'] for p in parsed_cube["ZAXIS"]["VECTOR"]]
        Voxel_axes = [X_voxel, Y_voxel, Z_voxel]

    # get the atom coordinates
    if (type(parsed_cube["GEOM"]) != list) :
        parsed_cube["GEOM"] = [parsed_cube["GEOM"]]

    Z_n = [atom["ATOMIC_NUMBER"] for atom in parsed_cube["GEOM"]]
    atom_charges = [atom["CHARGE"] for atom in parsed_cube["GEOM"]]
    atom_coords = [atom["POSITION"] if (world_unit == "bohr") else [p * pc['angstrom2bohrs'] for p in atom["POSITION"]] for atom in parsed_cube["GEOM"]]

    # construct the CUBE vector. Indexing is CUBE_vector[MO_ID][i*N_vals[1]*N_vals[2] + j*N_vals[2] + k] where i, j and k correspond to steps in the X, Y and Z voxel axes directions respectively.
    CUBE_vector = []
    for ID in range (N_vals):
        single_function = []
        for i in range(N_steps[0]):
            for j in range(N_steps[1]):
                for k in range(N_steps[2]):
                    single_function.append(parsed_cube["DATA"][i*N_steps[1]*N_steps[2]*N_vals + j*N_steps[2]*N_vals + k*N_vals + ID])
        CUBE_vector.append(single_function)

    cube_dict= {
            "CUBE_file": cube_path,
            "Header": {
                "N_atoms": N_atoms,
                "origin": origin,
                "N_steps": N_steps,
                "Voxel_axes": Voxel_axes,
                "Z_n": Z_n,
                "atom_charges": atom_charges,
                "atom_coords": atom_coords,
                "N_vals": N_vals
            },
            "ORB_IDS": orb_ids,
            "CUBE_data": CUBE_vector
    }
    return cube_dict
