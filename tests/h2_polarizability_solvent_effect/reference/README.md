# Checking the polarizability with finite differences

This folder contains input files for checking the ZZ component of the static
polarizability by finite differences:

- `h2_m.inp`: electric field with -0.01 field strength.
- `h2_p.inp`: electric field with +0.01 field strength.
- `h2.inp`: polarizability.

After running the 3 MRChem calculations, run the `check-pol.py` script to
extract and compare the finite difference and analytical polarizabilities.

