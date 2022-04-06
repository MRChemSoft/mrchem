# How to update the input parser automatically with utility script

All constants used by MRChem are defined in `python/mrhcem/physical_constants.MRChemPhysConstants`.
They are in turn completely derived from `CODATA2018` provided by NIST, which is interfaced via
`qcelemental`.
If a constant cannot be derived from existing ones, then it must be explicitly defined.
The NIST source in ASCII format can be found here: https://physics.nist.gov/cuu/Constants/Table/allascii.txt

All constants defined in `MRChemPhysConstants` need to be included in `template.yml` with default values.
The utility script `python/mrchem/update_input_parser.py` will read the current template and replace the existing
`Constants` section with a new one generated automatically by accessing `MRChemPhysConstants`.
The script will by default also execute `parselglossy` to actually update the input parser, and update the
user reference by copying the new `python/mrchem/input_parser/docs/user_ref.rst` to `doc/users/user_ref.rst`.

To perform all three actions, run the script as follows:

``` bash
cd python/mrchem
$ python update_input_parser.py
```

For a help message, run:

``` bash
cd python/mrchem
$ python update_input_parser.py -h
```

# How to update the input parser manually

Make your edits to the template file. 
NOTE: Constants must be defined in `python/mrchem/physical_constants.MRChemPhysConstants`, 
and should be added manually to the template.

Run parselglossy to update the input parser:

``` bash
$ cd python
$ parselglossy generate --template template.yml --docfile user_ref.rst --doc-header="User input reference" --target="mrchem/input_parser"
```
Copy the user reference file to update the documentation:

``` bash
cp python/mrchem/input_parser/docs/user_ref.rst doc/users/user_ref.rst
```
