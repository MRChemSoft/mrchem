import inspect
import json
import shlex
import subprocess
import sys
from functools import reduce
from math import log10
from optparse import OptionParser
from pathlib import Path

frame = inspect.stack()[-1]
module = inspect.getmodule(frame[0])
caller_file = module.__file__
caller_dir = Path(caller_file).resolve().parent

SUM_OCCUPIED = ("output", "properties", "orbital_energies", "sum_occupied")
E_KIN = ("output", "properties", "scf_energy", "E_kin")
E_EN = ("output", "properties", "scf_energy", "E_en")
E_EE = ("output", "properties", "scf_energy", "E_ee")
E_X = ("output", "properties", "scf_energy", "E_x")
E_XC = ("output", "properties", "scf_energy", "E_xc")
E_EEXT = ("output", "properties", "scf_energy", "E_eext")
E_NEXT = ("output", "properties", "scf_energy", "E_next")
E_EL = ("output", "properties", "scf_energy", "E_el")
ER_TOT = ("output", "properties", "scf_energy", "Er_tot")
ER_EL = ("output", "properties", "scf_energy", "Er_el")
ER_NUC = ("output", "properties", "scf_energy", "Er_nuc")


def DIPOLE_MOMENT(index):
    return ("output", "properties", "dipole_moment", f"dip-{index}", "vector")


def DIPOLE_MOMENT_EL(index):
    return ("output", "properties", "dipole_moment", f"dip-{index}",
            "vector_el")


def DIPOLE_MOMENT_NUC(index):
    return ("output", "properties", "dipole_moment", f"dip-{index}",
            "vector_nuc")


def QUADRUPOLE_MOMENT(index):
    return ("output", "properties", "quadrupole_moment", f"quad-{index}",
            "tensor")


def QUADRUPOLE_MOMENT_EL(index):
    return ("output", "properties", "quadrupole_moment", f"quad-{index}",
            "tensor_el")


def QUADRUPOLE_MOMENT_NUC(index):
    return ("output", "properties", "quadrupole_moment", f"quad-{index}",
            "tensor_nuc")


def POLARIZABILITY(frequency):
    return ("output", "properties", "polarizability", f"pol-{frequency:.6f}",
            "tensor")

def MAGNETIZABILITY(frequency):
    return ("output", "properties", "magnetizability", f"mag-{frequency:.6f}",
            "tensor")

def NMR_SHIELDING(atom):
    return ("output", "properties", "nmr_shielding", f"nmr-{atom}", "tensor")


def GEOMETRIC_DERIVATIVE(index, comp):
    return ("output", "properties", "geometric_derivative", f"geom-{index}", comp)

def HIRSHFELD_CHARGES(index, comp):
    return ("output", "properties", "hirshfeld_charges", f"hirshfeld-{index}", comp)

def run(options, *, input_file, filters=None, extra_args=None):
    launcher = "mrchem"
    launcher_full_path = Path(options.binary_dir).joinpath(launcher).resolve()

    input_file = caller_dir / input_file

    command = []
    command.append(str(launcher_full_path))
    command.append(str(input_file))
    command.append(f"--executable={options.binary_dir}/mrchem.x")
    if extra_args is not None:
        command += extra_args

    if options.launch_agent is not None:
        command.append(f"--launcher={options.launch_agent}")

    inp_no_suffix = input_file.stem
    output_prefix = inp_no_suffix

    sys.stdout.write(
        f"\nrunning {' '.join(command)}\ntest with input files {input_file} and args {extra_args}\n"
    )

    if options.skip_run:
        sys.stdout.write("(skipped run with -s|--skip-run)\n")
    else:
        child = subprocess.run(
            command,
            cwd=options.work_dir,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        # might still be the returncode is zero, but something happened
        if child.returncode != 0 or child.stderr:
            print(f"\nstdout\n{child.stdout}")
            print(f"\nstderr\n{child.stderr}")
            return 137

    if filters is None:
        sys.stdout.write("\nno filters passed: skipped automated verification\n")
        return 0

    computed = Path(options.work_dir) / f"{inp_no_suffix}.json"
    expected = caller_dir / f"reference/{inp_no_suffix}.json"
    with computed.open("r") as o_json, expected.open("r") as r_json:
        out = json.load(o_json)
        ref = json.load(r_json)

    success = True
    for what, threshold in filters.items():
        cptd = nested_get(out, what)
        xptd = nested_get(ref, what)
        where = location_in_dict(address=what)
        if options.no_verification:
            sys.stdout.write("\nskipped verification")
        else:
            passed, message = compare_values(cptd,
                                             xptd,
                                             where,
                                             rtol=threshold.rtol,
                                             atol=threshold.atol)
            sys.stdout.write(f"\n{message}")
            success = success and passed
    sys.stdout.write("\n")

    return 0 if success else 137


class Tolerance:
    def __init__(self, *, rtol=1e-5, atol=1e-8):
        self.rtol = rtol
        self.atol = atol


def rel_tolerance(rtol):
    """Mimics the behavior of the ``rel_tolerance`` filter in ``runtest``"""
    return Tolerance(rtol=rtol, atol=rtol)


def abs_tolerance(atol):
    """Mimics the behavior of the ``abs_tolerance`` filter in ``runtest``"""
    return Tolerance(rtol=0.0, atol=atol)


def is_close(computed, expected, *, rtol=1e-5, atol=1e-8, equal_phase=False):
    """
    .. code-block:: python
        absolute(computed - expected) <= (atol + rtol * absolute(expected))
    """
    def _is_close(c, x, r, a):
        return abs(c - x) <= (a + r * abs(x))

    if isinstance(computed, list):
        ok = [_is_close(c, x, rtol, atol) for c, x in zip(computed, expected)]
        allclose = all(ok)
    else:
        ok = _is_close(computed, expected, rtol, atol)
        allclose = ok

    if not allclose and equal_phase:
        # redo comparison negating the computed value
        if isinstance(computed, list):
            nok = [
                _is_close(-c, x, rtol, atol)
                for c, x in zip(computed, expected)
            ]
            allclose = all(nok)
        else:
            nok = _is_close(-c, x, rtol, atol)
            allclose = nok

    return allclose


def nested_get(d, ks):
    """Get value from a nested dictionary.

    Parameters
    ----------
    d : JSONDict
    ks : Tuple[str]

    Returns
    -------
    v : Optional[Any]

    Notes
    -----
    Adapted from: https://stackoverflow.com/a/40675868/2528668
    """
    def _func(x, k):
        return x.get(k, None) if isinstance(x, dict) else None

    return reduce(_func, ks, d)


def location_in_dict(*, address, dict_name="JSON"):
    """Convert tuple of keys of a ``JSONDict`` to its representation in code.

    For example, given ``("a", "b", "c")`` returns the string ``JSON['a']['b']['c']``.

    Parameters
    ----------
    address : Tuple[str]
    dict_name : str

    Returns
    -------
    where : str
    """
    return reduce(lambda x, y: x + f"['{y}']", address, dict_name)


def script_cli():
    parser = OptionParser(
        description="MRChem test runner. Heavily inspired by runtest.")

    parser.add_option(
        "--binary-dir",
        "-b",
        action="store",
        default=caller_dir,
        help="directory containing the binary/runscript [default: %default]",
    )
    parser.add_option(
        "--work-dir",
        "-w",
        action="store",
        default=caller_dir,
        help="working directory [default: %default]",
    )
    parser.add_option(
        "--launch-agent",
        "-l",
        action="store",
        default=None,
        help=
        'prepend a launch agent command (e.g. "mpirun -np 8" or "valgrind --leak-check=yes") [default: %default]',
    )
    parser.add_option(
        "--skip-run",
        "-s",
        action="store_true",
        default=False,
        help="skip actual calculation(s) [default: %default]",
    )
    parser.add_option(
        "--no-verification",
        "-n",
        action="store_true",
        default=False,
        help="run calculation(s) but do not verify results [default: %default]",
    )

    (options, args) = parser.parse_args(args=sys.argv[1:])

    return options


def compare_values(
    computed,
    expected,
    label,
    *,
    atol=1.0e-8,
    rtol=1.0e-5,
    equal_phase=False,
):
    """Returns True if two (lists of) floats are equal within a tolerance.

    Parameters
    ----------
    computed : float or float array-like
        Input value to compare against `expected`.
    expected : float or float array-like
        Reference value against which `computed` is compared.
    label : str, optional
        Label for passed and error messages.
    rtol : float, optional
        Relative tolerance (see formula below).
    equal_phase : bool, optional
        Compare computed *or its opposite* as equal.

    Returns
    -------
    allclose : bool
        Returns True if `expected` and `computed` are equal within tolerance; False otherwise.
    message : str
        Return passed or error message.

    Notes
    -----
    * Adapted from https://github.com/MolSSI/QCElemental/blob/master/qcelemental/testing.py
    * For scalar float-comparable types and for arbitrary-dimension, np.ndarray-castable, uniform-type,
      float-comparable types. For mixed types, use :py:func:`compare_recursive`.
    * Closeness measured as:
    .. code-block:: python
        absolute(computed - expected) <= (atol + rtol * absolute(expected))
    """
    pass_message = f"\t{label:.<120}PASSED"

    digits1 = abs(int(log10(atol))) + 2
    digits_str = f"to atol={atol}"
    if rtol > 1.0e-12:
        digits_str += f", rtol={rtol}"

    # check that type of computed and expected match, shortcircuit if not
    types_match = isinstance(computed, type(expected))
    if not types_match:
        return False, f"\tType of computed value ({type(computed).__name__}) does not match expected type ({type(expected).__name__})."

    expected_is_list = isinstance(expected, list)
    computed_is_list = isinstance(computed, list)

    # filter None from expected and computed lists
    if expected_is_list:
        expected = [_ for _ in expected if _ is not None]
    if computed_is_list:
        computed = [_ for _ in computed if _ is not None]

    ok = is_close(computed,
                  expected,
                  rtol=rtol,
                  atol=atol,
                  equal_phase=equal_phase)

    if ok:
        message = pass_message
    else:
        if not expected_is_list:
            xptd_str = f"{float(expected):.{digits1}f}"
        else:
            xptd_str = str(expected)
            xptd_str = "\n".join("    " + ln for ln in xptd_str.splitlines())

        if not computed_is_list:
            cptd_str = f"{float(computed):.{digits1}f}"
        else:
            cptd_str = str(computed)
            cptd_str = "\n".join("    " + ln for ln in cptd_str.splitlines())

        if not expected_is_list:
            diff = computed - expected
            diff_str = f"{float(diff):.{digits1}f}"
            message = f"\t{label}: computed value ({cptd_str}) does not match ({xptd_str}) {digits_str} by difference ({diff_str})."
        else:
            if len(computed) == len(expected):
                diff = [c - x for c, x in zip(computed, expected)]
                diff_str = str(diff)
                diff_str = "\n".join("    " + ln
                                     for ln in diff_str.splitlines())
                message = f"\t{label}: computed value does not match {digits_str}.\n  Expected:\n{xptd_str}\n  Observed:\n{cptd_str}\n  Difference:\n{diff_str}\n"
            else:
                message = f"\tDimension of computed value ({len(computed)}) does not match expected dimension ({len(expected)})."

    return ok, message
