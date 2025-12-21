import numpy as np
import argparse
import basis_set_exchange as bse

from pyscf import gto, scf
from pyscf.scf.atom_ks import get_atm_nrks

def uncontract_dfunctions(basname='AHGBS-7', atomic_symbol='He'):
    """Fetch AHGBS-7 from basis set exchange and uncontract all functions if contractions are present."""
    # Fetch basis set in NWChem format from BSE.
    bas = bse.get_basis(name=basname, elements=[atomic_symbol], fmt='nwchem', header=True, uncontract_segmented=True)
    
    return bas

def get_basis_in_mrchem_format(atomic_symbol):
    bas = bse.convert_formatted_basis_str(uncontract_dfunctions(atomic_symbol=atomic_symbol), 'nwchem', 'dalton')

    Z = bse.lut.element_Z_from_sym(atomic_symbol)
    if Z > 20:
        nshells = 3
    elif Z > 2:
        nshells = 2
    elif Z > 0:
        nshells = 1
    else:
        raise ValueError('Invalid atomic number!')
    nfuncs = ' '.join(['1' for _ in range(nshells)])
        

    # Ensure correct formatting for MRChem
    new = []
    for line in bas.splitlines():
        if line.strip().startswith('!'):
            continue
        elif line.strip() == '':
            continue
        elif line.strip().startswith('a'):
            continue
        elif line.strip().startswith('H'):
            new.append(' '.join(line.split()[1:]))
        else:
            new.append(' '.join(list(map(lambda x: f'{float(x):.12e}', line.split()))))

    b = []
    b.append('Gaussian basis 3-21G\n')
    b.append('1\n')
    b.append(f'{Z}. 1 {nshells} {nfuncs}\n')
    b.append(f'{atomic_symbol} {0:.12f} {0:.12f} {0:.12f}\n')
    for line in new:
        b.append(line + '\n')

    return "".join(b)


def sad(atomic_symbol):
    """Compute density matrix for passed atom, and write density and basis set to file."""
    # Build molecule input
    basis = uncontract_dfunctions(atomic_symbol=atomic_symbol, basname='AHGBS-7')
    mol = gto.Mole(atom=f'{atomic_symbol} 0.0 0.0 0.0', basis=basis, symmetry=False, verbose=False)
    mol.spin = mol.tot_electrons() % 2
    mol.build()

    # Run DFT calculation for atom
    atm_scf_result = get_atm_nrks(mol, xc='slater', grid=(250, 434))

    # Compute AO density matrix with some thresholding
    scf_energy, _, mo_coeff, mo_occ = atm_scf_result[atomic_symbol]
    D = mo_coeff @ np.diag(mo_occ) @ mo_coeff.T
    zeros = np.abs(D) <= 1e-13
    D[zeros] = 0.0
    dim = D.shape[0]

    # Flatten density matrix and get basis set in mrchem format
    D_fmt = [f'{x:.12e}' for _ in D.tolist() for x in _]
    basis_mrc = get_basis_in_mrchem_format(atomic_symbol)

    # Write files
    with open(atomic_symbol + '.dens', 'w') as d, open(atomic_symbol + '.bas', 'w') as b:
        d.write(f'{dim}\n')
        for line in D_fmt:
            d.write(line + '\n')

        b.write(basis_mrc)
        
    return scf_energy
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate SAD guess files for MRChem')
    parser.add_argument('element', help='Atomic symbol for which you want to generate SAD files.')
    args = parser.parse_args()
    
    e = sad(args.element)
    print(f'{"LDA/AHGBS-7 energy for":15}: {args.element}')
    print(f'{"Total energy":15}: {e:.8f} a.u.')
    print(f'{"Density file":15}: {args.element}.dens')
    print(f'{"Basis set file":15}: {args.element}.bas')
