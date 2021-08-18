
#
# MRChem, a numerical real-space code for molecular electronic structure
# calculations within the self-consistent field (SCF) approximations of quantum
# chemistry (Hartree-Fock and Density Functional Theory).
# Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

import input_parser.plumbing.pyparsing as pp
from json import dump

def write_cube_dict(file_dict):
    all_path_list = sort_paths(file_dict["CUBEfiles"]) # this is a temporary fix until I have decided how to order the different files
    all_cube_list = []
    for path_list in all_path_list:
        cube_list = []
        if (len(path_list) != 0):
            for path in path_list:
                cube_list.append(parse_cube_file(path))
        all_cube_lists.append(cube_list)

    with open("CUBE_p_vector.json", "w") as fd:
        dump(all_cube_list[0], fd, indent=2)

    with open("CUBE_a_vector.json", "w") as fd:
        dump(all_cube_list[1], fd, indent=2)

    with open("CUBE_b_vector.json", "w") as fd:
        dump(all_cube_list[2], fd, indent=2)


def sort_paths(path_l):
    path_p = []
    path_a = []
    path_b = []
    if (len(path_l) != 0):
        for path in path_l:
            path_s = path.split("_")
            if ("p" == path_s[1]):
                path_p.append(path)
            elif("a" == path_s[1]):
                path_a.append(path)
            elif("b" == path_s[1]):
                path_b.append(path)
            else:
                raise ValueError("Wrong CUBE file name format")
    return [path_p, path_a, path_b]


#TODO do a sanity check on the naming of the files

def parse_cube_file(cube_path):

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

    # volumetric data
    voxel_t = pp.delimitedList(float_t, delim=pp.Empty())("DATA")


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

    cube_t = preamble_t(before_t, after_t) + voxel_t

    # parse the whole file and extract both all the values and the header
    with open(cube_path, "r") as cube_file:
        # we skip the first two lines as they are only comments and put the rest of the file in a single string
        cube_str = "\n".join(cube_file.readlines()[2:])

    parsed_cube = cube_t.parseString(cube_str).asDict()


    # start extracting the header values


    # get cube origin data
    N_atoms = parsed_cube["NATOMS"]
    origin = parsed_cube["ORIGIN"]

    # Set the amount of values depending on if the DSET_IDs were present or not
    if ("DSET_IDS" in parsed_cube.keys):
        N_vals = len(parsed_cube["DSET_IDS"])
    else:
        N_vals = parsed_cube["NVAL"]

    # files are given as [phi,rho]_[p,a,b]_[rsp,scf]_#_[re,im]_[x,y].cube
    #test that the file name makes sense
    path_ids = cube_path.split("_")[3]
    if ("-" in path_ids):
        from_to = path_ids.split("-")
        orb_ids = list(range(from_to[0], (from_to[1]+1), 1))  # we work as including the from and to, so 4-7 includes 4, 5, 6 and 7.
        if (N_vals != len(orb_ids)):
            raise ValueError("Different amount of orbitals in file and amount of orbitals specified in file name.")
    else:
        orb_ids = [float(path_ids)]


    # get voxel axis data to construct a basis for the cube space.
    N_steps = [parsed_cube["XAXIS"]["NVOXELS"], parsed_cube["YAXIS"]["NVOXELS"], parsed_cube["ZAXIS"]["NVOXELS"]]
    Voxel_axes = [parsed_cube["XAXIS"]["VECTOR"], parsed_cube["YAXIS"]["VECTOR"], parsed_cube["ZAXIS"]["VECTOR"]]

    # get the atom coordinates
    Z_n = [atom["ATOMIC_NUMBER"] for atom in parsed_cube["GEOM"]]
    atom_charges = [atom["CHARGE"] for atom in parsed_cube["GEOM"]]
    atom_coords = [atom["POSITION"] for atom in parsed_cube["GEOM"]]

    # construct the CUBE vector. Indexing is CUBE_vector[MO_ID][i*N_vals[1]*N_vals[2] + j*N_vals[2] + k] where i, j and k correspond to steps in the X, Y and Z voxel axes directions respectively.
    CUBE_vector = [ [parsed_cube["DATA"][i*N_steps[1]*N_steps[2]*N_vals + j*N_steps[2]*N_vals + k*N_vals + ID] for i in range(N_steps[0]) for j in range(N_steps[1]) for k in range(N_steps[2])]  for ID in range(N_vals)]

    cube_dict= {
            "CUBE_file": cube_path,
            "Header": {
                "comments": comments[0]+comments[1],
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
