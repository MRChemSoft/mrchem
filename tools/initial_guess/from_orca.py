#!/usr/bin/env python

import sys
import json

#
# CLI script to generate gto initial guess for MRChem
#
# Usage:
#
# python from_orca.py path/to/orca/output.json
#
# This generates the required mrchem.bas and mrchem.mop/moa/mob files
#
# Example orca 6 inputfile with orbital localization:
#

"""
!PBE cc-pVTZ Bohrs
%pal nprocs 4 end

%loc
  LocMet NewBoys
  Occ true
end

* XYZ 0 1
H      -1.4375  0.0000   1.0250
O       0.0000  0.0000  -0.1250
H       1.4375  0.0000   1.0250
*
"""

#
# After running with orca you will have a .gbw file with the canonical orbitals
# and a .loc file with the localized orbitals. You can run the orca_2json
# CLI-tool on these to convert to a .json file that can be given to this script
#

def consolidate_basis(basis):
    new_basis = {}

    for shell in basis:
        angmom = shell["Shell"]
        exp_hash = str(shell["Exponents"])

        if angmom not in new_basis:
            new_basis[angmom] = {}

        if exp_hash not in new_basis[angmom]:
            new_basis[angmom][exp_hash] = (shell["Exponents"],
                                           [shell["Coefficients"]])
        else:
            new_basis[angmom][exp_hash][1].append(shell["Coefficients"])

    return {angmom: list(shell_basis.values())
            for angmom, shell_basis in new_basis.items()}


def print_mrchem_bas_file(atoms, basis_sets, coord_scale):
    with open("mrchem.bas", "w") as f:
        f.write("Gaussian basis from ORCA gbw json file\n")
        f.write(f"{len(basis_sets):9}\n")

        for symbol, q, basis, atom_idxs in basis_sets:
            f.write(f"{q:9} {len(atom_idxs):4}")
            f.write(f" {len(basis):4}")
            for shells in basis.values():
                f.write(f" {len(shells):4}")
            f.write("\n")

            for i in atom_idxs:
                coord = atoms[i]["Coords"]
                f.write(f"{symbol} ")
                for x in coord:
                    f.write(f" {x * coord_scale:18.10f}")
                f.write("\n")

            for shells in basis.values():
                for exps, coeffs in shells:
                    n_prim = len(exps)
                    n_cont = len(coeffs)

                    f.write(f"{n_prim:9} {n_cont:4}\n")
                    for i in range(n_prim):
                        f.write(f"{exps[i]:15.7f}")
                        for j in range(n_cont):
                            f.write(f" {coeffs[j][i]:11.8f}")
                        f.write("\n")


def organize_basis(atoms):
    # Basis as key, atom index as value
    basis_sets = {}

    for atom in atoms:
        basis_hash = str((atom["ElementLabel"], atom["Basis"]))
        if basis_hash not in basis_sets:
            basis_sets[basis_hash] = (
                atom["ElementLabel"], atom["NuclearCharge"],
                consolidate_basis(atom["Basis"]), [atom["Idx"]]
            )
        else:
            basis_sets[basis_hash][3].append(atom["Idx"])

    return list(basis_sets.values())


def invert_perm(perm):
    invperm = [-1 for _ in perm]

    for i, j in enumerate(perm):
        invperm[j] = i

    return invperm


angmom_to_l = {
    's': 0,
    'p': 1,
    'd': 2,
    'f': 3,
    'g': 4,
    'h': 5,
}


def shell_permutation(angmom):
    if angmom == 'p':
        return [1, 2, 0]  # Hardcode xyz order for p orbitals

    l = angmom_to_l[angmom]

    m_to_i = {m: i for i, m in enumerate(range(-l, l+1))}

    perm = [m_to_i[0]]

    for m in range(l):
        perm.append(m_to_i[m + 1])
        perm.append(m_to_i[-(m + 1)])

    return invert_perm(perm)


def atom_basis_permutation(atom_basis):
    perm = []

    for angmom in atom_basis:
        shell_perm = shell_permutation(angmom)

        num_funcs = 0
        for _, coeffs in atom_basis[angmom]:
            num_funcs += len(coeffs)

        for _ in range(num_funcs):
            offset = len(perm)
            for i in shell_perm:
                perm.append(i + offset)

    return perm


def make_mo_permutation(atoms, basis_sets):
    atoms_perm = []

    for _, _, _, idxs in basis_sets:
        for i in idxs:
            atoms_perm.append(i)

    atom_basis_perms = [
        atom_basis_permutation(basis) for _, _, basis, _ in basis_sets
    ]

    atom_to_basis = [-1 for _ in atoms_perm]

    for i, (_, _, _, idxs) in enumerate(basis_sets):
        for j in idxs:
            atom_to_basis[j] = i

    atoms_offsets = []
    offset = 0
    for i in atom_to_basis:
        atoms_offsets.append(offset)
        offset += len(atom_basis_perms[i])

    perm = []

    for i in atoms_perm:
        for j in atom_basis_perms[atom_to_basis[i]]:
            perm.append(j + atoms_offsets[i])

    return perm


def permute(array, perm):
    return [array[i] for i in perm]


def read_occ_mo(mol, mo_perm):
    restricted = mol["HFTyp"] == "RHF"

    mo_buf = [0.0 for _ in mo_perm]

    moa = []

    i = 0

    orbitals = mol["MolecularOrbitals"]["MOs"]

    while orbitals[i]["Occupancy"] > 0:
        moa.extend(permute(orbitals[i]["MOCoefficients"], mo_perm))
        i += 1

    if restricted:
        return (moa, None)

    while orbitals[i]["Occupancy"] == 0:
        i += 1

    mob = []

    while orbitals[i]["Occupancy"] > 0:
        mob.extend(permute(orbitals[i]["MOCoefficients"], mo_perm))
        i += 1

    return (moa, mob)


def print_mrchem_mo_file(filename, mo, nao):
    with open(filename, "w") as f:
        f.write(f"{nao:12}\n")
        for x in mo:
            f.write(f"{x:20.15f}\n")


def print_mrchem_mo_files(mo, nao):
    moa, mob = mo

    if mob is None:
        print_mrchem_mo_file("mrchem.mop", moa, nao)
    else:
        print_mrchem_mo_file("mrchem.moa", moa, nao)
        print_mrchem_mo_file("mrchem.mob", mob, nao)


if __name__ == "__main__":
    with open(sys.argv[1]) as orca_file:
        orca_json = json.load(orca_file)

    mol = orca_json["Molecule"]

    atoms = mol["Atoms"]

    basis_sets = organize_basis(atoms)

    if mol["CoordinateUnits"] == "Bohrs":
        coord_scale = 1.0
    elif mol["CoordinateUnits"] == "Angs":
        coord_scale = 1.8897261339212517

    print_mrchem_bas_file(atoms, basis_sets, coord_scale)

    mo_perm = make_mo_permutation(atoms, basis_sets)

    mo = read_occ_mo(mol, mo_perm)

    print_mrchem_mo_files(mo, len(mo_perm))
