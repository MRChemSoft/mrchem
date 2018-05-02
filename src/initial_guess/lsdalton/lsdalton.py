import os
import re
from periodictable import PeriodicTableByZ as PTz
from periodictable import PeriodicTableByName as PTn

def makexyz(coords):
    nAtoms = len(coords)
    xyz = "    %s\n" % (nAtoms)

    for line in coords[:]:
        sp = line.split()
        elm = sp[0].lower()
        x, y, z = map(float, sp[1:])
        d = PTn[elm].Z
        xyz += "%s       % 10.10f    % 10.10f    % 10.10f\n" % (elm.title(), x, y, z)
    xyz +="\n"
    #print xyz
    return xyz

def makedal(coords, method, charge, mult):
    nel = getNElectrons(coords, charge)
    npaired = nel - (mult - 1)
    nalpha = npaired/2 + (mult - 1)
    nbeta = npaired/2

    lsdal = "**GENERAL\n"
    lsdal += ".NOGCBASIS\n"
    lsdal += "**WAVE FUNCTIONS\n"
    if method == "HF":
        lsdal += ".HF\n"
    else:
        lsdal += ".DFT\n"
        lsdal += method + "\n"
    lsdal += "*DENSOPT\n"
    lsdal += ".DIIS\n"
    lsdal += ".UNREST\n"
    lsdal += ".NALPHA\n"
    lsdal += "%s\n" % nalpha
    lsdal += ".NBETA\n"
    lsdal += "%s\n" % nbeta
    #lsdal += "**INTEGRAL\n"
    #lsdal += ".UNCONT\n"
    lsdal += "*END OF INPUT\n"
    #print lsdal
    return lsdal

def xyz2mol(coords, basis):
    nAtoms = len(coords)
    lsmol = "BASIS\n%s\nTitle 1\nTitle 2\n" % basis
    lsmol += "Atomtypes=%s\n" % nAtoms

    for line in coords[:]:
        sp = line.split()
        elm = sp[0].lower()
        x, y, z = map(float, sp[1:])
        d = PTn[elm].Z
        lsmol += "Charge=%s.0 Atoms=1\n" % d
        lsmol += "%s       % 10.10f    % 10.10f    % 10.10f\n" % (elm.title(), x, y, z)
    lsmol +="\n"
    #print lsmol
    return lsmol

def xyz2bas(coords, basis, basdir):
    baspath = os.path.join(basdir, basis)
    atomTypes = getNAtomTypes(coords)
    mrbas = "Gaussian basis %s\n" % (basis)
    mrbas += "        %s\n" % (atomTypes)

    match = ''
    N = 0
    middle = ''
    for line in coords[:]:
        name = line.split()[0]
        if (match != '' and name != match):
            top = getBas(baspath)[match].splitlines()[0]
            times = int(getBas(baspath)[match].splitlines()[0].split()[1])
            top = ("        %s    %s    %s"+times*"    1"+"\n") % \
                    (top.split()[0],N ,top.split()[1])
            bottom =''.join(str(num)+"\n" for num in
                    getBas(baspath)[match].splitlines()[2:])
            N = 0
            mrbas += top+middle+bottom
            middle = ''
        balls = line.split()
        middle += "%s       % 5.10f      % 5.10f      % 5.10f\n" % (balls[0], float(balls[1]),
               float(balls[2]), float(balls[3]))
        N += 1
        match = name
    times = int(getBas(baspath)[match].splitlines()[0].split()[1])
    top = getBas(baspath)[match].splitlines()[0]
    top = ("        %s    %s    %s"+times*"    1"+"\n") % \
            (top.split()[0],N ,top.split()[1])
    bottom =''.join(str(num)+"\n" for num in
            getBas(baspath)[match].splitlines()[2:])
    mrbas += top+middle+bottom
    #print mrbas
    return mrbas

def getNElectrons(coords, charge):
    nel = 0
    for line in coords[:]:
        sp = line.split()
        elm = sp[0].lower()
        x, y, z = map(float, sp[1:])
        nel += PTn[elm].Z
    return (nel - charge)

def getNAtomTypes(coords):
    atoms = set()
    for line in coords[:]:
        elm = line.split()[0]
        atoms.add(elm)
    return len(atoms)

def getBas(baspath):
    with open(baspath, "r") as f:
        lines = f.readlines()
    return genBas(lines)

def genBas(lines):
    molbas = dict()
    atombas = ''
    lqn = 0
    z = 1
    header = """        %s.    %s
    %s     %%4.10f    %%4.10f    %%4.10f
"""

    for l in lines:
        if re.match("A", l, re.I):
            if len(atombas) > 0:
                molbas[PTz[z].symbol] = (header + atombas[:-1]) % \
                    (z, lqn, PTz[z].symbol)
                z = int(l.split()[1])
                atombas = ''
                lqn = 0
        elif re.match(r'( +[0-9]+){3}', l):
            lqn += 1
            x = l.split()[0:2]
            atombas += "        " + x[0] + "    " + x[1] + "\n"
        elif re.match(r'[ \t]*$', l):
            pass
        elif not l[0] == "$":
            atombas += l
    molbas[PTz[z].symbol] = (header + atombas[:-1]) % \
        (z, lqn, PTz[z].symbol)
    return molbas

#def main():
    #run_initial_guess('coord.xyz', '3-21G')
#    xyz2bas(datadir+'/coord.xyz', '3-21G')
#    xyz2mol(datadir+'/coord.xyz', '3-21G')

#if __name__ == '__main__':
    #main()
# vim:ft=python
