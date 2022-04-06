#!/usr/bin/env python

from pathlib import Path
from ruamel.yaml import YAML
import subprocess
import os
import shutil
import argparse

from physical_constants import MRChemPhysConstants

root = Path.cwd().parent
target = Path(__file__).parent.joinpath('input_parser')

yaml = YAML()

def update_constants():
    pc = MRChemPhysConstants()
    f_template = root.joinpath("template.yml")
    template = yaml.load(f_template)

    new = {
        "keywords": template["keywords"], 
        "sections": [section for section in template['sections'] if section['name'] != "Constants"]}

    # Build new constants section
    constants = {
        "name": "Constants", 
        "docstring": "Physical and mathematical constants used by MRChem", 
        "keywords": []
        }
    for name, unit, value, docstring in pc.data:
        constants["keywords"].append({
            "name": name,
            "default": float(value),
            "type": "float",
            "docstring": f"{docstring} (unit: {unit})"
        })

    new["sections"].append(constants)
    yaml.dump(new, f_template)


def run_parselglossy():
    os.chdir(root)
    cmd = [
        'parselglossy', 
        'generate', 
        '--template', 'template.yml', 
        '--docfile', 'user_ref.rst', 
        '--doc-header', '\"User input reference\"',
        '--target', target
         ]

    subprocess.call(cmd)

def update_doc():
    src = target.joinpath('docs/user_ref.rst')
    dst = root.parent.joinpath('doc/users/user_ref.rst')
    shutil.copyfile(src, dst)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sT', '--skip-template', action='store_true', help="Do not update constants section in template.yml")
    parser.add_argument('-sP', '--skip-parselglossy', action='store_true', help='Do not update input parser with parselglossy')
    parser.add_argument('-sD', '--skip-doc', action='store_true', help='Do not update the user reference file')
    args = parser.parse_args()

    if not args.skip_template:
        print(f'{"Updating template":20} ... ', end='')
        update_constants()
        print('done')

    if not args.skip_parselglossy:
        print(f'{"Running parselglossy":20} ... ', end='')
        run_parselglossy()
        print('done')

    if not args.skip_doc:
        print(f'{"Updating user ref":20} ... ', end='')
        update_doc()
        print('done')
