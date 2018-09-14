#!/usr/bin/env python
"""
Add or update the license information in the header of source files.

The script reads in the .gitattributes file, located in the project root
directory, to figure out which files need to be inspected and which license
header they need to have.
For example:

src/pedra/pedra_dlapack.F90 !licensefile
src/solver/*.hpp licensefile=.githooks/LICENSE-C++

The first line specifies that the file in src/pedra/pedra_dlapack.F90 should
not be touched, while the second line states that all .hpp files in src/solver
should get an header from the template in .githooks/LICENSE-C++
Location of files in .gitattributes are always specified with respect
to the project root directory.

The script reads in the appropriate license header template and prepares an
header with the correct year and authors information. These are encoded in the
variables YEAR and AUTHORS. The latter has to be modified by hand.
"""

import glob
import os
import re
import shutil
import sys
import tempfile
from datetime import date


def add_header(filepath, header, TEXT, YEAR, AUTHORS):
    """
    Add or update header in source file
    """
    tmpdir = tempfile.gettempdir()
    tmpfil = os.path.join(tmpdir, os.path.basename(filepath) + '.bak')
    shutil.copy2(filepath, tmpfil)
    with open(tmpfil, 'r') as tmp:
        inpt = tmp.readlines()
        output = []

        # Check if header is already present
        present = re.compile(TEXT)
        if list(filter(present.search, inpt)):
            # Check if year and authors in current file are up to date
            toupdate = re.compile(r'{0:s} (?!{1:d} {2:s}).*\n'.format(
                'Copyright \(C\)', YEAR, AUTHORS))
            if list(filter(toupdate.search, inpt)):
                print(('Updating header in {:s}'.format(filepath)))
                # Check to preserve '#!' at the top of the file
                if len(inpt) > 0 and inpt[0].startswith('#!'):
                    output.append('{}\n'.format(inpt[0]))
                    inpt = inpt[1:]
                regex = re.compile(r'Copyright \(C\).*\n')
                repl = r'Copyright (C) {0:d} {1:s}\n'.format(YEAR, AUTHORS)
                output.extend([re.sub(regex, repl, x) for x in inpt])
        else:
            print(('Adding header in {:s}'.format(filepath)))
            # Check to preserve '#!' at the top of the file
            if len(inpt) > 0 and inpt[0].startswith('#!'):
                output.append('{:s}\n'.format(inpt[0]))
                inpt = inpt[1:]
            output.append('{:s}\n\n'.format(header))
            for line in inpt:
                output.append(line)

        if output:
            try:
                f = open(filepath, 'w')
                f.writelines(output)
            except IOError as err:
                print(('Something went wrong trying to add header to {:s}: {:s}'.
                       format(filepath, err)))
            finally:
                f.close()
        os.remove(tmpfil)


def prepare_header(stub, YEAR, AUTHORS):
    """
    Update year and author information in license header template
    """
    with open(stub, 'r') as l:
        header = l.read().rstrip()
        # Insert correct YEAR and AUTHORS in stub
        rep = {'YEAR': str(YEAR), 'AUTHORS': AUTHORS}
        rep = dict((re.escape(k), v) for k, v in rep.items())
        pattern = re.compile("|".join(list(rep.keys())))
        header = pattern.sub(lambda m: rep[re.escape(m.group(0))], header)
    return header


def file_license(attributes):
    """
    Obtain dictionary { file : license } from .gitattributes
    """
    file_license = {}
    with open(attributes, 'r') as f:
        # Read in .gitattributes
        tmp = f.read()
        # Removing all comment lines and other attributes
        pattern = re.compile(r'(?m)^\#.*\n?|^((?!licensefile).)*$')
        gitattributes = re.sub(pattern, '', tmp).split()
        # Obtain list of files
        fil = [x for x in gitattributes if not 'licensefile' in x]
        # Remove licensefile= from strings
        lic = [
            re.sub(r'licensefile\=', '', x) for x in gitattributes
            if 'licensefile' in x
        ]
        # Create list of blacklisted files
        blacklist = [
            fname for key, value in list(dict(list(zip(fil, lic))).items())
            if value == '!licensefile' for fname in glob.glob(key)
        ]
        # Now create a dictionary with the files to be considered for
        # license header manipulation
        file_license = {
            key: value
            for k, value in list(dict(list(zip(fil, lic))).items())
            for key in glob.glob(k) if key not in blacklist
        }
    return file_license


def license_maintainer(TEXT, AUTHORS):
    """
    Maintain license header in source files
    """
    YEAR = date.today().year

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    headerize = file_license(os.path.join(project_root_dir, '.gitattributes'))

    for fname, license in list(headerize.items()):
        # Prepare header
        header = prepare_header(
            os.path.join(project_root_dir, license), YEAR, AUTHORS)
        add_header(fname, header, TEXT, YEAR, AUTHORS)


if __name__ == '__main__':
    TEXT = sys.argv[1]
    AUTHORS = sys.argv[2]
    license_maintainer(TEXT, AUTHORS)
