# How to update the input parser

Run:

``` bash
$ cd python
$ parselglossy generate --template template.yml --docfile user_ref.rst --doc-header="User input reference" --target="mrchem/input_parser"
```
Remember to also update the documentation:

``` bash
cp python/mrchem/input_parser/docs/user_ref.rst doc/users/user_ref.rst
```
