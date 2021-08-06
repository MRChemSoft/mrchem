
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

import re
from json import dump

def write_cube_dict(file_dict):
    cube_path_list = file_dict["CUBEfiles"]
    if (len(cube_path_list) != 0):
        cube_list= []
        for cube_path in cube_path_list:
            cube_list.append(parse_cube_file(cube_path))
        with open("CUBE_p_vector.json", "w") as fd:
            dump(cube_list, fd, indent=2)


def parse_cube_file(cube_path):
    # set up the regex to fetch all voxel values
    voxel_values = re.compile("((-| )\w\.\d+(E-|E\+)\d+)")

    # parse the whole file and extract both all the values and the header
    cube_list = []
    header_list = []
    with open(cube_path, "r") as cube_file:
        for line in cube_file:
            matches = re.findall(voxel_values, line)
            if (len(matches) == 0):
                header_list.append(line)
            else:
                voxel_list = [float(match[0]) for match in matches]
                cube_list.extend(voxel_list)
    # start extracting the header values
    # get comments
    comments = header_list[:2]

    # get cube origin data
    origin_line = header_list[2].split()
    N_atoms =int(origin_line[0])
    origin = list(map(float, origin_line[1:4]))

    # get voxel axis data to construct a basis for the cube space.
    N_steps = []
    Voxel_axes = []
    for line in header_list[3:6]:
        ls = line.split()
        N_steps.append(int(ls[0]))
        Voxel_axes.append(list(map(float, ls[1:])))

    # get the atom coordinates
    Z_n = []
    atom_charges = []
    atom_coords = []
    for line in header_list[6:6+abs(N_atoms)]:
        ls = line.split()
        Z_n.append(int(ls[0]))
        atom_charges.append(float(ls[1]))
        atom_coords.append(list(map(float, ls[2:])))

    # Set the amount of values depending on if the DSET_IDs were present or not
    if (N_atoms < 0):
    # fetch all DSET_IDs, of these, only the first value is necessary, as it tells me the amount of values per voxel point. the other values might be important in the future, so I am fetching everything.
        DSETIDs = []
        for line in header_list[6+abs(N_atoms):]:
            DSETIDs.extend(list(map(int, line.split()[:])))

        N_vals = DSETIDs[0]

    else:
        N_vals = int(origin_line[-1])

    # construct the CUBE vector. Indexing is CUBE_vector[MO_ID][i*N_vals[1]*N_vals[2] + j*N_vals[2] + k] where i, j and k correspond to steps in the X, Y and Z voxel axes directions respectively.
    CUBE_vector = [ [cube_list[i*N_steps[1]*N_steps[2]*N_vals + j*N_steps[2]*N_vals + k*N_vals + ID] for i in range(N_steps[0]) for j in range(N_steps[1]) for k in range(N_steps[2])]  for ID in range(N_vals)]

    cube_dict= {
            "CUBE_file": cube_path,
            "Header": {
                "comments": comments[0]+comments[1],
                "N_atoms": abs(N_atoms),
                "origin": origin,
                "N_steps": N_steps,
                "Voxel_axes": Voxel_axes,
                "Z_n": Z_n,
                "atom_charges": atom_charges,
                "atom_coords": atom_coords,
                "N_vals": N_vals
            },
            "CUBE_data": CUBE_vector
    }
    return cube_dict
