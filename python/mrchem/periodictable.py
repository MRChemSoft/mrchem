#
# MRChem, a numerical real-space code for molecular electronic structure
# calculations within the self-consistent field (SCF) approximations of quantum
# chemistry (Hartree-Fock and Density Functional Theory).
# Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

from collections import OrderedDict, namedtuple


class Element(
    namedtuple(
        "Element",
        "radius covalent Z mass symbol bpt mpt density volume name debye a crystal cpera conf r_rms",
    )
):
    __slots__ = ()

    def __new__(cls, iterable):
        return super(cls, Element).__new__(
            cls,
            radius=float(iterable[0]),
            covalent=float(iterable[1]),
            Z=int(iterable[2]),
            mass=float(iterable[3]),
            symbol=iterable[4],
            bpt=float(iterable[5]),
            mpt=float(iterable[6]),
            density=float(iterable[7]),
            volume=float(iterable[8]),
            name=iterable[9],
            debye=float(iterable[10]),
            a=float(iterable[11]),
            crystal=iterable[12],
            cpera=float(iterable[13]),
            conf=iterable[14],
            r_rms=iterable[15], #RMS radius for smeared nucleus, from Visscher, L.; Dyall, K. At. Data Nucl. Data Tables 1997, 67, 207– 224
        )

    def __str__(self):
        return "{:s} ({:s}): {{{:s}}}, Z={:d}, m={:f}".format(
            self.name, self.symbol, self.conf, self.Z, self.mass
        )


# fmt: off
PeriodicTable = OrderedDict({
    'h':
    Element([
        '0.79', '0.32', '1', '1.00794', 'H', '20.268', '14.025', '0.0899',
        '14.4', 'Hydrogen', '110.0', '3.75', 'HEX', '1.731', '1s1', '2.6569547399E-05'
    ]),
    'he':
    Element([
        '0.49', '0.93', '2', '4.002602', 'He', '4.215', '0.95', '0.1787',
        '0.0', 'Helium', '-26.0', '3.57', 'HEX', '1.633', '1s2', '3.5849373401E-05'
    ]),
    'li':
    Element([
        '2.05', '1.23', '3', '6.941', 'Li', '1615', '453.7', '0.53', '13.10',
        'Lithium', '400.0', '3.49', 'BCC', '0.00', '1s2_2s1', '4.0992133976E-05'
    ]),
    'be':
    Element([
        '1.40', '0.90', '4', '9.012182', 'Be', '2745', '1560.0', '1.85', '5.0',
        'Beryllium', '1000.0', '2.29', 'HEX', '1.567', '1s2_2s2', '4.3632829651E-05 '
    ]),
    'b':
    Element([
        '1.17', '0.82', '5', '10.811', 'B', '4275', '2300.0', '2.34', '4.6',
        'Boron', '1250.0', '8.73', 'TET', '0.576', '1s2_2s2_2p1', '4.5906118608E-05'
    ]),
    'c':
    Element([
        '0.91', '0.77', '6', '12.011', 'C', '4470.0', '4100.0', '2.62', '4.58',
        'Carbon', '1860.0', '3.57', 'DIA', '0.00', '1s2_2s2_2p2', '4.6940079496E-05'
    ]),
    'n':
    Element([
        '0.75', '0.75', '7', '14.00674', 'N', '77.35', '63.14', '1.251',
        '17.3', 'Nitrogen', '-79.0', '4.039', 'HEX', '1.651', '1s2_2s2_2p3', '4.8847128967E-05'
    ]),
    'o':
    Element([
        '0.65', '0.73', '8', '15.9994', 'O', '90.18', '50.35', '1.429', '14.0',
        'Oxygen', '-46.0', '6.83', 'CUB', '0.00', '1s2_2s2_2p4', '5.0580178957E-05'
    ]),
    'f':
    Element([
        '0.57', '0.72', '9', '18.9984032', 'F', '84.95', '53.48', '1.696',
        '17.1', 'Fluorine', '0.0', '0.00', 'MCL', '0.00', '1s2_2s2_2p5', '5.2927138943E-05'
    ]),
    'ne':
    Element([
        '0.51', '0.71', '10', '20.1797', 'Ne', '27.096', '24.553', '0.901',
        '16.7', 'Neon', '63.0', '4.43', 'FCC', '0.00', '1s2_2s2_2p6', '5.3654104231E-05'
    ]),
    'na':
    Element([
        '2.23', '1.54', '11', '22.989768', 'Na', '1156', '371.0', '0.97',
        '23.7', 'Sodium', '150.0', '4.23', 'BCC', '0.00', '[Ne]3s1', '5.5699159416E-05'
    ]),
    'mg':
    Element([
        '1.72', '1.36', '12', '24.3050', 'Mg', '1363', '922', '1.74', '13.97',
        'Magnesium', '318.0', '3.21', 'HEX', '1.624', '[Ne]3s2', '5.6341070732E-05'
    ]),
    'al':
    Element([
        '1.82', '1.18', '13', '26.981539', 'Al', '2793', '933.25', '2.70',
        '10.0', 'Aluminum', '394.0', '4.05', 'FCC', '0.00', '[Ne]3s2_3p1', '5.8165765928E-05'
    ]),
    'si':
    Element([
        '1.46', '1.11', '14', '28.0855', 'Si', '3540.0', '1685', '2.33',
        '12.1', 'Silicon', '625.0', '5.43', 'DIA', '0.00', '[Ne]3s2_3p2', '5.8743802504E-05'
    ]),
    'p':
    Element([
        '1.23', '1.06', '15', '30.97362', 'P', '550.0', '317.30', '1.82',
        '17.0', 'Phosphorus', '0.0', '7.17', 'CUB', '0.00', '[Ne]3s2_3p3', '6.0399312923E-05'
    ]),
    's':
    Element([
        '1.09', '1.02', '16', '32.066', 'S', '717.75', '388.36', '2.07',
        '15.5', 'Sulfur', '0.0', '10.47', 'ORC', '0.00', '[Ne]3s2_3p4', '6.0927308666E-05'
    ]),
    'cl':
    Element([
        '0.97', '0.99', '17', '35.4527', 'Cl', '239.1', '172.16', '3.17',
        '22.7', 'Chlorine', '0.0', '6.24', 'ORC', '0.00', '[Ne]3s2_3p5', '6.2448101115E-05'
    ]),
    'ar':
    Element([
        '0.88', '0.98', '18', '39.948', 'Ar', '87.30', '83.81', '1.784',
        '28.5', 'Argon', '85.0', '5.26', 'FCC', '0.00', '[Ne]3s2_3p6', '6.4800211825E-05'
    ]),
    'k':
    Element([
        '2.77', '2.03', '19', '39.0983', 'K', '1032', '336.35', '0.86',
        '45.46', 'Potassium', '100.0', '5.23', 'BCC', '0.00', '[Ar]4s1', '6.4346167051E-05'
    ]),
    'ca':
    Element([
        '2.23', '1.91', '20', '40.078', 'Ca', '1757', '1112', '1.55', '29.9',
        'Calcium', '230.0', '5.58', 'FCC', '0.00', '[Ar]4s2', '6.4800211825E-05'
    ]),
    'sc':
    Element([
        '2.09', '1.62', '21', '44.955910', 'Sc', '3104', '1812', '3.0', '15.0',
        'Scandium', '-359.0', '3.31', 'HEX', '1.594', '[Ar]3d1_4s2', '6.6963627201E-05'
    ]),
    'ti':
    Element([
        '2.00', '1.45', '22', '47.88', 'Ti', '3562', '1943', '4.50', '10.64',
        'Titanium', '380.0', '2.95', 'HEX', '1.588', '[Ar]3d2_4s2', '6.8185577480E-05'
    ]),
    'v':
    Element([
        '1.92', '1.34', '23', '50.9415', 'V', '3682', '2175', '5.8', '8.78',
        'Vanadium', '390.0', '3.02', 'BCC', '0.00', '[Ar]3d3_4s2', '6.9357616830E-05'
    ]),
    'cr':
    Element([
        '1.85', '1.18', '24', '51.9961', 'Cr', '2945', '2130.0', '7.19',
        '7.23', 'Chromium', '460.0', '2.88', 'BCC', '0.00', '[Ar]3d5_4s1', '6.9738057221E-05'
    ]),
    'mn':
    Element([
        '1.79', '1.17', '25', '54.93085', 'Mn', '2335', '1517', '7.43', '1.39',
        'Manganese', '400.0', '8.89', 'CUB', '0.00', '[Ar]3d5_4s2', '7.0850896638E-05'
    ]),
    'fe':
    Element([
        '1.72', '1.17', '26', '55.847', 'Fe', '3135', '1809', '7.86', '7.1',
        'Iron', '460.0', '2.87', 'BCC', '0.00', '[Ar]3d6_4s2', '7.1212829817E-05'
    ]),
    'co':
    Element([
        '1.67', '1.16', '27', '58.93320', 'Co', '3201', '1768', '8.90', '6.7',
        'Cobalt', '385.0', '2.51', 'HEX', '0.00', '[Ar]3d7_4s2', '7.2273420879E-05'
    ]),
    'ni':
    Element([
        '1.62', '1.15', '28', '58.69', 'Ni', '3187', '1726', '8.90', '6.59',
        'Nickel', '375.0', '3.52', 'FCC', '0.00', '[Ar]3d8_4s2', '7.1923970253E-05'
    ]),
    'cu':
    Element([
        '1.57', '1.17', '29', '63.546', 'Cu', '2836', '1357.6', '8.96', '7.1',
        'Copper', '315.0', '3.61', 'FCC', '0.00', '[Ar]3d10_4s1', '7.3633018675E-05'
    ]),
    'zn':
    Element([
        '1.53', '1.25', '30', '65.39', 'Zn', '1180.0', '692.73', '7.14', '9.2',
        'Zinc', '234.0', '2.66', 'HEX', '0.00', '[Ar]3d10_4s2', '7.3963875193E-05'
    ]),
    'ga':
    Element([
        '1.81', '1.26', '31', '69.723', 'Ga', '2478', '302.90', '5.91', '11.8',
        'Gallium', '240.0', '4.51', 'ORC', '0.00', '[Ar]3d10_4s2_4p1', '7.5568424848E-05'
    ]),
    'ge':
    Element([
        '1.52', '1.22', '32', '72.61', 'Ge', '3107', '1210.4', '5.32', '13.6',
        'Germanium', '360.0', '5.66', 'DIA', '0.00', '[Ar]3d10_4s2_4p2', '7.7097216161E-05'
    ]),
    'as':
    Element([
        '1.33', '1.20', '33', '74.92159', 'As', '876', '1081', '5.72', '13.1',
        'Arsenic', '285.0', '4.13', 'RHL', '54.16', '[Ar]3d10_4s2_4p3', '7.7394645153E-05'
    ]),
    'se':
    Element([
        '1.22', '1.16', '34', '78.96', 'Se', '958', '494', '4.80', '16.45',
        'Selenium', '-150.0', '4.36', 'HEX', '0.00', '[Ar]3d10_4s2_4p4', '7.8843427408E-05'
    ]),
    'br':
    Element([
        '1.12', '1.14', '35', '79.904', 'Br', '332.25', '265.90', '3.12',
        '23.5', 'Bromine', '0.0', '6.67', 'ORC', '0.00', '[Ar]3d10_4s2_4p5', '7.8558604038E-05'
    ]),
    'kr':
    Element([
        '1.03', '1.12', '36', '83.80', 'Kr', '119.80', '115.78', '3.74',
        '38.9', 'Krypton', '-73.0', '5.72', 'FCC', '0.00', '[Ar]3d10_4s2_4p6', '7.9959560033E-05'
    ]),
    'rb':
    Element([
        '2.98', '2.16', '37', '85.4678', 'Rb', '961', '312.64', '1.53', '55.9',
        'Rubidium', '-56.0', '5.59', 'BCC', '0.00', '[Kr]5s1', '8.0233033713E-05'
    ]),
    'sr':
    Element([
        '2.45', '1.91', '38', '87.62', 'Sr', '1650.0', '1041', '2.6', '33.7',
        'Strontium', '-147.0', '6.08', 'FCC', '0.00', '[Kr]5s2', '8.1040799081E-05'
    ]),
    'y':
    Element([
        '2.27', '1.62', '39', '88.90585', 'Y', '3611', '1799', '4.5', '19.8',
        'Yttrium', '-256.0', '3.65', 'HEX', '1.571', '[Kr]4d1_5s2', '8.1305968993E-05'
    ]),
    'zr':
    Element([
        '2.16', '1.45', '40', '91.224', 'Zr', '4682', '2125', '6.49', '14.1',
        'Zirconium', '250.0', '3.23', 'HEX', '1.593', '[Kr]4d2_5s2', '8.1569159980E-05'
    ]),
    'nb':
    Element([
        '2.09', '1.34', '41', '92.90638', 'Nb', '5017', '2740.0', '8.55',
        '10.87', 'Niobium', '275.0', '3.30', 'BCC', '0.00', '[Kr]4d4_5s1', '8.2347219223E-05'
    ]),
    'mo':
    Element([
        '2.01', '1.30', '42', '95.94', 'Mo', '4912', '2890.0', '10.2', '9.4',
        'Molybdenum', '380.0', '3.15', 'BCC', '0.00', '[Kr]4d5_5s1', '8.3607614434E-05'
    ]),
    'tc':
    Element([
        '1.95', '1.27', '43', '-98', 'Tc', '4538', '2473', '11.5', '8.5',
        'Technetium', '0.0', '2.74', 'HEX', '1.604', '[Kr]4d5_5s2', '8.3607614434E-05'
    ]),
    'ru':
    Element([
        '1.89', '1.25', '44', '101.07', 'Ru', '4423', '2523', '12.2', '8.3',
        'Ruthenium', '-382.0', '2.70', 'HEX', '1.584', '[Kr]4d7_5s1', '8.4585397905E-05'
    ]),
    'rh':
    Element([
        '1.83', '1.25', '45', '102.90550', 'Rh', '3970.0', '2236', '12.4',
        '8.3', 'Rhodium', '-350.0', '3.80', 'FCC', '0.00', '[Kr]4d8_5s1', '8.4825835954E-05'
    ]),
    'pd':
    Element([
        '1.79', '1.28', '46', '106.42', 'Pd', '3237', '1825', '12.0', '8.9',
        'Palladium', '275.0', '3.89', 'FCC', '0.00', '[Kr]4d10_5s0', '8.5537941156E-05'
    ]),
    'ag':
    Element([
        '1.75', '1.34', '47', '107.8682', 'Ag', '2436', '1234', '10.5', '10.3',
        'Silver', '215.0', '4.09', 'FCC', '0.00', '[Kr]4d10_5s1', '8.5772320442E-05'
    ]),
    'cd':
    Element([
        '1.71', '1.48', '48', '112.411', 'Cd', '1040.0', '594.18', '8.65',
        '13.1', 'Cadmium', '120.0', '2.98', 'HEX', '1.886', '[Kr]4d10_5s2', '8.7373430179E-05'
    ]),
    'in':
    Element([
        '2.00', '1.44', '49', '114.82', 'In', '2346', '429.76', '7.31', '15.7',
        'Indium', '129.0', '4.59', 'TET', '1.076', '[Kr]4d10_5s2_5p1', '8.7596760865E-05'
    ]),
    'sn':
    Element([
        '1.72', '1.41', '50', '118.710', 'Sn', '2876', '505.06', '7.30',
        '16.3', 'Tin', '170.0', '5.82', 'TET', '0.546', '[Kr]4d10_5s2_5p2', '8.8694413774E-05'
    ]),
    'sb':
    Element([
        '1.53', '1.40', '51', '121.75', 'Sb', '1860.0', '904', '6.68', '18.23',
        'Antimony', '200.0', '4.51', 'RHL', '57.10', '[Kr]4d10_5s2_5p3', '8.8910267995E-05'
    ]),
    'te':
    Element([
        '1.42', '1.36', '52', '127.60', 'Te', '1261', '722.65', '6.24', '20.5',
        'Tellurium', '-139.0', '4.45', 'HEX', '1.33', '[Kr]4d10_5s2_5p4', '9.0801452955E-05'
    ]),
    'i':
    Element([
        '1.32', '1.33', '53', '126.90447', 'I', '458.4', '386.7', '4.92',
        '25.74', 'Iodine', '0.0', '7.27', 'ORC', '0.00', '[Kr]4d10_5s2_5p5', '9.0181040290E-05'
    ]),
    'xe':
    Element([
        '1.24', '1.31', '54', '131.29', 'Xe', '165.03', '161.36', '5.89',
        '37.3', 'Xenon', '-55.0', '6.20', 'FCC', '0.00', '[Kr]4d10_5s2_5p6', '9.1209776425E-05'
    ]),
    'cs':
    Element([
        '3.34', '2.35', '55', '132.90543', 'Cs', '944', '301.55', '1.87',
        '71.07', 'Cesium', '-40.0', '6.05', 'BCC', '0.00', '[Xe]6s1', '9.1412392742E-05'
    ]),
    'ba':
    Element([
        '2.78', '1.98', '56', '137.327', 'Ba', '2171', '1002', '3.5', '39.24',
        'Barium', '-110.0', '5.02', 'BCC', '0.00', '[Xe]6s2', '9.2410525664E-05'
    ]),
    'la':
    Element([
        '2.74', '1.69', '57', '138.9055', 'La', '3730.0', '1193', '6.7',
        '20.73', 'Lanthanum', '132.0', '3.75', 'HEX', '1.619', '[Xe]5d1_6s2', '9.2607247118E-05'
    ]),
    'hf':
    Element([
        '2.16', '1.44', '72', '178.49', 'Hf', '4876', '2500.0', '13.1', '13.6',
        'Hafnium', '0.0', '3.20', 'HEX', '1.582', '[Xe]4f14_5d2_6s2', '9.9970978172E-05'
    ]),
    'ta':
    Element([
        '2.09', '1.34', '73', '180.9479', 'Ta', '5731', '3287', '16.6',
        '10.90', 'Tantalum', '225.0', '3.31', 'BCC', '0.00', '[Xe]4f14_5d3_6s2', '1.0013585755E-04'
    ]),
    'w':
    Element([
        '2.02', '1.30', '74', '183.85', 'W', '5828', '3680.0', '19.3', '9.53',
        'Tungsten', '310.0', '3.16', 'BCC', '0.00', '[Xe]4f14_5d4_6s2', '1.0062688070E-04'
    ]),
    're':
    Element([
        '1.97', '1.28', '75', '186.207', 'Re', '5869', '3453', '21.0', '8.85',
        'Rhenium', '416.0', '2.76', 'HEX', '1.615', '[Xe]4f14_5d5_6s2', '1.0111259523E-04'
    ]),
    'os':
    Element([
        '1.92', '1.26', '76', '190.2', 'Os', '5285', '3300.0', '22.4', '8.49',
        'Osmium', '-400.0', '2.74', 'HEX', '1.579', '[Xe]4f14_5d6_6s2', '1.0191070333E-04'
    ]),
    'ir':
    Element([
        '1.87', '1.27', '77', '192.22', 'Ir', '4701', '2716', '22.5', '8.54',
        'Iridium', '430.0', '3.84', 'FCC', '0.00', '[Xe]4f14_5d7_6s2', '1.0206865731E-04'
    ]),
    'pt':
    Element([
        '1.83', '1.30', '78', '195.08', 'Pt', '4100.0', '2045', '21.4', '9.10',
        'Platinum', '230.0', '3.92', 'FCC', '0.00', '[Xe]4f14_5d10_6s0', '1.0238293593E-04'
    ]),
    'au':
    Element([
        '1.79', '1.34', '79', '196.96654', 'Au', '3130.0', '1337.58', '19.3',
        '10.2', 'Gold', '170.0', '4.08', 'FCC', '0.00', '[Xe]4f14_5d10_6s1', '1.0269507292E-04'
    ]),
    'hg':
    Element([
        '1.76', '1.49', '80', '200.59', 'Hg', '630.0', '234.28', '13.53',
        '14.82', 'Mercury', '100.0', '2.99', 'RHL', '70.75',
        '[Xe]4f14_5d10_6s2', '1.0346628039E-04'
    ]),
    'tl':
    Element([
        '2.08', '1.48', '81', '204.3833', 'Tl', '1746', '577', '11.85', '17.2',
        'Thallium', '96.0', '3.46', 'HEX', '1.599', '[Xe]4f14_5d10_6s2_6p1', '1.0392291259E-04'
    ]),
    'pb':
    Element([
        '1.81', '1.47', '82', '207.2', 'Pb', '2023', '600.6', '11.4', '18.17',
        'Lead', '88.0', '4.95', 'FCC', '0.00', '[Xe]4f14_5d10_6s2_6p2', '1.0437511130E-04'
    ]),
    'bi':
    Element([
        '1.63', '1.46', '83', '208.98037', 'Bi', '1837', '544.52', '9.8',
        '21.3', 'Bismuth', '120.0', '4.75', 'RHL', '57.23',
        '[Xe]4f14_5d10_6s2_6p3', '1.0452487744E-04'
    ]),
    'po':
    Element([
        '1.53', '1.46', '84', '-209', 'Po', '1235', '527', '9.4', '22.23',
        'Polonium', '0.0', '3.35', 'SC', '0.00', '[Xe]4f14_5d10_6s2_6p4', '1.0452487744E-04'
    ]),
    'at':
    Element([
        '1.43', '1.45', '85', '210.0', 'At', '610.0', '575', '0.0', '0.0',
        'Astatine', '0.0', '0.00', '', '0.00', '[Xe]4f14_5d10_6s2_6p5', '1.0467416660E-04'
    ]),
    'rn':
    Element([
        '1.34', '1.43', '86', '-222', 'Rn', '211', '202', '9.91', '50.5',
        'Radon', '0.0', '0.00', 'FCC', '0.00', '[Xe]4f14_5d10_6s2_6p6', '1.0642976299E-04'
    ]),
    'fr':
    Element([
        '3.50', '2.50', '87', '-223', 'Fr', '950.0', '300.0', '0.0', '0.0',
        'Francium', '0.0', '0.00', 'BCC', '0.00', '[Rn]7s1', '1.0657317899E-04'
    ]),
    'ra':
    Element([
        '3.00', '2.40', '88', '226.025', 'Ra', '1809', '973', '5', '45.20',
        'Radium', '0.0', '0.00', '', '0.00', '[Rn]7s2', '1.0700087100E-04'
    ]),
    'ac':
    Element([
        '3.20', '2.20', '89', '227.028', 'Ac', '3473', '1323', '10.07',
        '22.54', 'Actinium', '0.0', '5.31', 'FCC', '0.00', '[Rn]6d1_7s2', '1.0714259349E-04'
    ]),
    'rf':
    Element([
        '0.0', '0.0', '104', '-257.0', 'Rf', '0.0', '0.0', '0.0', '0.0',
        'Rutherfordium', '0.0', '0.00', '', '0.00', '4-5s', '1.1173204420E-04'
    ]),
    'db':
    Element([
        '0.0', '0.0', '105', '-262.0', 'Db', '0.0', '0.0', '0.0', '0.0',
        'Dubnium', '0.0', '0.00', '', '0.00', '40s', '1.1186082063E-04'
    ]),
    'sg':
    Element([
        '0.0', '0.0', '106', '-263.0', 'Sg', '0.0', '0.0', '0.0', '0.0',
        'Seaborgium', '0.0', '0.00', '', '0.00', '0.9s', '1.1198926979E-04'
    ]),
    'bh':
    Element([
        '0.0', '0.0', '107', '-262.0', 'Bh', '0.0', '0.0', '0.0', '0.0',
        'Bohrium', '0.0', '0.00', '', '0.00', '2ms', '1.1186082063E-04'
    ]),
    'hs':
    Element([
        '0.0', '0.0', '108', '-264.0', 'Hs', '0.0', '0.0', '0.0', '0.0',
        'Hassium', '0.0', '0.00', '', '0.00', '', '1.1224519460E-04'
    ]),
    'mt':
    Element([
        '0.0', '0.0', '109', '-266', 'Mt', '0.0', '0.0', '0.0', '0.0',
        'Meitnerium', '0.0', '0.00', '', '0.00', '5ms', '1.1237267433E-04'
    ]),
    '110':
    Element([
        '0.0', '0.0', '110', '-269', '110', '0.0', '0.0', '0.0', '0.0',
        '(recent_disc.)', '0.0', '0.00', '', '0.00', '', '-1.0'
    ]),
    '111':
    Element([
        '0.0', '0.0', '111', '-272', '111', '0.0', '0.0', '0.0', '0.0',
        '(recent_disc.)', '0.0', '0.00', '', '0.00', '4/1000s', '-1.0'
    ]),
    '112':
    Element([
        '0.0', '0.0', '112', '-277', '112', '0.0', '0.0', '0.0', '0.0',
        '(recent_disc.)', '0.0', '0.00', '', '0.00', '280\265s', '-1.0'
    ]),
    'ce':
    Element([
        '2.70', '1.65', '58', '140.115', 'Ce', '3699', '1071', '6.78', '20.67',
        'Cerium', '-139.0', '5.16', 'FCC', '0.00', '[Xe]4f2_5d0_6s2', '9.2803027311E-05'
    ]),
    'pr':
    Element([
        '2.67', '1.65', '59', '140.90765', 'Pr', '3785', '1204', '6.77',
        '20.8', 'Praseodymium', '-152.0', '3.67', 'HEX', '1.614',
        '[Xe]4f3_5d0_6s2', '9.2997877424E-05'
    ]),
    'nd':
    Element([
        '2.64', '1.64', '60', '144.24', 'Nd', '3341', '1289', '7.00', '20.6',
        'Neodymium', '-157.0', '3.66', 'HEX', '1.614', '[Xe]4f4_5d0_6s2', '9.3576955934E-05'
    ]),
    'pm':
    Element([
        '2.62', '1.63', '61', '-145', 'Pm', '3785', '1204', '6.475', '22.39',
        'Promethium', '0.0', '0.00', '', '0.00', '[Xe]4f5_5d0_6s2', '9.3768193375E-05'
    ]),
    'sm':
    Element([
        '2.59', '1.62', '62', '150.36', 'Sm', '2064', '1345', '7.54', '19.95',
        'Samarium', '166.0', '9.00', 'RHL', '23.22', '[Xe]4f6_5d0_6s2', '9.5082839751E-05'
    ]),
    'eu':
    Element([
        '2.56', '1.85', '63', '151.965', 'Eu', '1870.0', '1090.0', '5.26',
        '28.9', 'Europium', '-107.0', '4.61', 'BCC', '0.00', '[Xe]4f7_5d0_6s2', '9.5267329183E-05'
    ]),
    'gd':
    Element([
        '2.54', '1.61', '64', '157.25', 'Gd', '3539', '1585', '7.89', '19.9',
        'Gadolinium', '-176.0', '3.64', 'HEX', '1.588', '[Xe]4f7_5d1_6s2', '9.6177915369E-05'
    ]),
    'tb':
    Element([
        '2.51', '1.59', '65', '158.92534', 'Tb', '3496', '1630.0', '8.27',
        '19.2', 'Terbium', '-188.0', '3.60', 'HEX', '1.581', '[Xe]4f9_5d0_6s2', '9.6357719009E-05'
    ]),
    'dy':
    Element([
        '2.49', '1.59', '66', '162.50', 'Dy', '2835', '1682', '8.54', '19.0',
        'Dysprosium', '-186.0', '3.59', 'HEX', '1.573', '[Xe]4f10_5d0_6s2', '9.6892647152E-05'
    ]),
    'ho':
    Element([
        '2.47', '1.58', '67', '164.93032', 'Ho', '2968', '1743', '8.80',
        '18.7', 'Holmium', '-191.0', '3.58', 'HEX', '1.570', '[Xe]4f11_5d0_6s2', '9.6892647152E-05'
    ]),
    'er':
    Element([
        '2.45', '1.57', '68', '167.26', 'Er', '3136', '1795', '9.05', '18.4',
        'Erbium', '-195.0', '3.56', 'HEX', '1.570', '[Xe]4f12_5d0_6s2', '9.7943009317E-05'
    ]),
    'tm':
    Element([
        '2.42', '1.56', '69', '168.93421', 'Tm', '2220.0', '1818', '9.33',
        '18.1', 'Thulium', '-200.0', '3.54', 'HEX', '1.570', '[Xe]4f13_5d0_6s2', '9.8115626740E-05'
    ]),
    'yb':
    Element([
        '2.40', '1.74', '70', '173.04', 'Yb', '1467', '1097', '6.98', '24.79',
        'Ytterbium', '-118.0', '5.49', 'FCC', '0.00', '[Xe]4f14_5d0_6s2', '9.8968651305E-05'
    ]),
    'lu':
    Element([
        '2.25', '1.56', '71', '174.967', 'Lu', '3668', '1936', '9.84', '17.78',
        'Lutetium', '-207.0', '3.51', 'HEX', '1.585', '[Xe]4f14_5d1_6s2', '9.9137288835E-05'
    ]),
    'th':
    Element([
        '3.16', '1.65', '90', '232.0381', 'Th', '5061', '2028', '11.7', '19.9',
        'Thorium', '100.0', '5.08', 'FCC', '0.00', '[Rn]6d2_7s2', '1.0784503195E-04'
    ]),
    'pa':
    Element([
        '3.14', '0.0', '91', '231.03588', 'Pa', '0.0', '0.0', '15.4', '15.0',
        'Protactinium', '0.0', '3.92', 'TET', '0.825', '[Rn]5f2_6d1_7s2', '1.0770535752E-04'
    ]),
    'u':
    Element([
        '3.11', '1.42', '92', '238.0289', 'U', '4407', '1405', '18.90',
        '12.59', 'Uranium', '-210.0', '2.85', 'ORC', '0.00', '[Rn]5f3_6d1_7s2', '1.0867476102E-04'
    ]),
    'np':
    Element([
        '3.08', '0.0', '93', '237.048', 'Np', '0.0', '910.0', '20.4', '11.62',
        'Neptunium', '-188.0', '4.72', 'ORC', '0.00', '[Rn]5f4_6d1_7s2', '1.0853744903E-04'
    ]),
    'pu':
    Element([
        '3.05', '0.0', '94', '-244', 'Pu', '3503', '913', '19.8', '12.32',
        'Plutonium', '-150.0', '0.00', 'MCL', '0.00', '[Rn]5f6_6d0_7s2', '1.0949065967E-04'
    ]),
    'am':
    Element([
        '3.02', '0.0', '95', '-243', 'Am', '2880.0', '1268', '13.6', '17.86',
        'Americium', '0.0', '0.00', '', '0.00', '[Rn]5f7_6d0_7s2', '1.0935561268E-04'
    ]),
    'cm':
    Element([
        '2.99', '0.0', '96', '-247', 'Cm', '0.0', '1340.0', '13.511', '18.28',
        'Curium', '0.0', '0.00', '', '0.00', '[Rn]5f7_6d1_7s2', '1.0989359973E-04'
    ]),
    'bk':
    Element([
        '2.97', '0.0', '97', '-247', 'Bk', '0.0', '0.0', '0.0', '0.0',
        'Berkelium', '0.0', '0.00', '', '0.00', '[Rn]5f8_6d1_7s2', '1.0989359973E-04'
    ]),
    'cf':
    Element([
        '2.95', '0.0', '98', '-251', 'Cf', '0.0', '900.0', '0.0', '0.0',
        'Californium', '0.0', '0.00', '', '0.00', '[Rn]5f10_6d0_7s2', '1.1042580946E-04'
    ]),
    'es':
    Element([
        '2.92', '0.0', '99', '-252', 'Es', '0.0', '0.0', '0.0', '0.0',
        'Einsteinium', '0.0', '0.00', '', '0.00', '[Rn]5f11_6d0_7s2', '1.1055797721E-04'
    ]),
    'fm':
    Element([
        '2.90', '0.0', '100', '-257', 'Fm', '0.0', '0.0', '0.0', '0.0',
        'Fermium', '0.0', '0.00', '', '0.00', '[Rn]5f12_6d0_7s2', '1.1121362374E-04'
    ]),
    'md':
    Element([
        '2.87', '0.0', '101', '-258', 'Md', '0.0', '0.0', '0.0', '0.0',
        'Mendelevium', '0.0', '0.00', '', '0.00', '[Rn]5f13_6d0_7s2', '1.1134373034E-04'
    ]),
    'no':
    Element([
        '2.85', '0.0', '102', '-259', 'No', '0.0', '0.0', '0.0', '0.0',
        'Nobelium', '0.0', '0.00', '', '0.00', '[Rn]5f14_6d0_7s2', '1.1147350119E-04'
    ]),
    'lr':
    Element([
        '2.82', '0.0', '103', '-260', 'Lr', '0.0', '0.0', '0.0', '0.0',
        'Lawrencium', '0.0', '0.00', '', '0.00', '[Rn]5f14_6d1_7s2', '1.1186082063E-04'
    ]),
    'none':
    Element([
        '0.0', '0.0', '0', '0', 'None', '0.0', '0.0', '0.0', '0.0', 'None',
        '0.0', '0.00', '', '0.00', '-', '-1.0'
    ]),
    'x':
    Element([
        '0.0', '0.0', '0', '0', 'X', '0.0', '0.0', '0.0', '0.0', 'X', '0.0',
        '0.00', '', '0.00', '-', '-1.0'
    ]),
    'q':
    Element([
        '0.0', '0.0', '0', '0', 'Q', '0.0', '0.0', '0.0', '0.0', 'Q', '0.0',
        '0.00', '', '0.00', '-', '-1.0'
    ])
})
# fmt: on

PeriodicTableByName = PeriodicTable

PeriodicTableByZ = OrderedDict(
    {
        32: PeriodicTableByName["ge"],
        64: PeriodicTableByName["gd"],
        31: PeriodicTableByName["ga"],
        57: PeriodicTableByName["la"],
        3: PeriodicTableByName["li"],
        81: PeriodicTableByName["tl"],
        69: PeriodicTableByName["tm"],
        103: PeriodicTableByName["lr"],
        90: PeriodicTableByName["th"],
        22: PeriodicTableByName["ti"],
        52: PeriodicTableByName["te"],
        65: PeriodicTableByName["tb"],
        43: PeriodicTableByName["tc"],
        73: PeriodicTableByName["ta"],
        70: PeriodicTableByName["yb"],
        0: PeriodicTableByName["none"],
        66: PeriodicTableByName["dy"],
        54: PeriodicTableByName["xe"],
        1: PeriodicTableByName["h"],
        15: PeriodicTableByName["p"],
        0: PeriodicTableByName["x"],
        30: PeriodicTableByName["zn"],
        111: PeriodicTableByName["111"],
        110: PeriodicTableByName["110"],
        112: PeriodicTableByName["112"],
        63: PeriodicTableByName["eu"],
        40: PeriodicTableByName["zr"],
        68: PeriodicTableByName["er"],
        44: PeriodicTableByName["ru"],
        75: PeriodicTableByName["re"],
        104: PeriodicTableByName["rf"],
        88: PeriodicTableByName["ra"],
        37: PeriodicTableByName["rb"],
        86: PeriodicTableByName["rn"],
        45: PeriodicTableByName["rh"],
        4: PeriodicTableByName["be"],
        56: PeriodicTableByName["ba"],
        83: PeriodicTableByName["bi"],
        97: PeriodicTableByName["bk"],
        35: PeriodicTableByName["br"],
        6: PeriodicTableByName["c"],
        19: PeriodicTableByName["k"],
        8: PeriodicTableByName["o"],
        16: PeriodicTableByName["s"],
        74: PeriodicTableByName["w"],
        76: PeriodicTableByName["os"],
        27: PeriodicTableByName["co"],
        96: PeriodicTableByName["cm"],
        17: PeriodicTableByName["cl"],
        20: PeriodicTableByName["ca"],
        91: PeriodicTableByName["pa"],
        98: PeriodicTableByName["cf"],
        58: PeriodicTableByName["ce"],
        48: PeriodicTableByName["cd"],
        55: PeriodicTableByName["cs"],
        24: PeriodicTableByName["cr"],
        29: PeriodicTableByName["cu"],
        59: PeriodicTableByName["pr"],
        78: PeriodicTableByName["pt"],
        94: PeriodicTableByName["pu"],
        82: PeriodicTableByName["pb"],
        71: PeriodicTableByName["lu"],
        46: PeriodicTableByName["pd"],
        84: PeriodicTableByName["po"],
        61: PeriodicTableByName["pm"],
        108: PeriodicTableByName["hs"],
        67: PeriodicTableByName["ho"],
        105: PeriodicTableByName["db"],
        72: PeriodicTableByName["hf"],
        80: PeriodicTableByName["hg"],
        2: PeriodicTableByName["he"],
        101: PeriodicTableByName["md"],
        12: PeriodicTableByName["mg"],
        5: PeriodicTableByName["b"],
        9: PeriodicTableByName["f"],
        42: PeriodicTableByName["mo"],
        25: PeriodicTableByName["mn"],
        7: PeriodicTableByName["n"],
        109: PeriodicTableByName["mt"],
        23: PeriodicTableByName["v"],
        89: PeriodicTableByName["ac"],
        47: PeriodicTableByName["ag"],
        77: PeriodicTableByName["ir"],
        95: PeriodicTableByName["am"],
        13: PeriodicTableByName["al"],
        33: PeriodicTableByName["as"],
        18: PeriodicTableByName["ar"],
        79: PeriodicTableByName["au"],
        85: PeriodicTableByName["at"],
        49: PeriodicTableByName["in"],
        28: PeriodicTableByName["ni"],
        102: PeriodicTableByName["no"],
        11: PeriodicTableByName["na"],
        41: PeriodicTableByName["nb"],
        60: PeriodicTableByName["nd"],
        10: PeriodicTableByName["ne"],
        99: PeriodicTableByName["es"],
        93: PeriodicTableByName["np"],
        107: PeriodicTableByName["bh"],
        87: PeriodicTableByName["fr"],
        21: PeriodicTableByName["sc"],
        26: PeriodicTableByName["fe"],
        100: PeriodicTableByName["fm"],
        53: PeriodicTableByName["i"],
        38: PeriodicTableByName["sr"],
        106: PeriodicTableByName["sg"],
        0: PeriodicTableByName["q"],
        36: PeriodicTableByName["kr"],
        14: PeriodicTableByName["si"],
        92: PeriodicTableByName["u"],
        50: PeriodicTableByName["sn"],
        62: PeriodicTableByName["sm"],
        39: PeriodicTableByName["y"],
        51: PeriodicTableByName["sb"],
        34: PeriodicTableByName["se"],
    }
)


def main():
    print("PeriodicTableByZ = {")
    for k, v in PeriodicTableByName.items():
        print("{:d} : PeriodicTableByName['{:s}'],".format(v.Z, k))
    print("}")


if __name__ == "__main__":
    main()
