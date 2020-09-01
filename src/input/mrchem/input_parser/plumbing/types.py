# -*- coding: utf-8 -*-
#
# parselglossy -- Generic input parsing library, speaking in tongues
# Copyright (C) 2020 Roberto Di Remigio, Radovan Bast, and contributors.
#
# This file is part of parselglossy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# For information on the complete list of contributors to the
# parselglossy library, see: <http://parselglossy.readthedocs.io/>
#

import re
from typing import Any, Callable, List, Optional, Tuple, Union

from typing import Dict  # noqa: F401; noqa: F401

RunCallable = Callable[[Any, str], Tuple[str, Any]]

ScalarTypes = Union[bool, str, int, float, complex]

allowed_scalar_types = ["str", "int", "float", "complex", "bool"]


ListTypes = Union[List[bool], List[str], List[int], List[float], List[complex]]

allowed_list_types = [f"List[{t}]" for t in allowed_scalar_types]

AllowedTypes = Union[ScalarTypes, ListTypes]

allowed_types = allowed_scalar_types + allowed_list_types


def type_matches(value: AllowedTypes, expected_type: str) -> Optional[bool]:
    """Checks whether a value is of the expected type.

    Parameters
    ----------
    value : AllowedTypes
        Value whose type needs to be checked
    expected_type : str

    Notes
    -----
    Allowed types T are: `str`, `int`, `float`, `complex`, `bool`,
    as well as `List[T]`.

    Returns
    -------
    True if value has the type expected_type, otherwise False.

    Raises
    ------
    ValueError
        If expected_type is not among the allowed types.

    Notes
    -----
    Complex numbers are a tad more fastidious, as they *might be* read in as
    strings. To avoid false negatives, we have a nasty hack, for which RDR will
    forever burn in hell. For complex and List[complex] we will re-run type
    checking after type
    casting if and only if:
      - We expected complex and List[complex].
      - We got a str or List[str].
      - The string are in the right format for the type casting operator.
    """

    # first verify whether expected_type is allowed
    if expected_type not in allowed_types:
        raise ValueError(f"could not recognize expected_type: {expected_type}")

    expected_complex = expected_type == "complex" or expected_type == "List[complex]"

    expected_type_is_list = re.search(r"^List\[(\w+)\]$", expected_type)

    matches = False
    if expected_type_is_list is not None:
        matches = _type_check_list(
            value, expected_type_is_list.group(1)
        )  # type: ignore
    else:
        matches = _type_check_scalar(value, expected_type)  # type: ignore

    recheck_complex = (
        expected_complex
        and _scalar_type_is_str(value)
        and _str_is_complex_number(value)
    )
    if not matches and recheck_complex:
        # Run the type cast, this might fail with ValueError and return None.
        value = _complex_helper(value)  # type: ignore
        if value is not None:
            matches = type_matches(value, expected_type)  # type:ignore
        else:
            matches = False

    return matches


def _complex_helper(
    value: Union[complex, List[complex]]
) -> Optional[Union[complex, List[complex]]]:
    try:
        if type(value).__name__ == "list":
            return type_fixers["List[complex]"](value)
        else:
            return type_fixers["complex"](value)
    except ValueError:
        return None


def _scalar_type_is_str(value: AllowedTypes) -> bool:
    return any(
        [_type_check_scalar(value, "str"), _type_check_list(value, "str")]
    )  # type: ignore


def _str_is_complex_number(value: AllowedTypes) -> bool:
    cmplx = re.compile(r"^.*\d+j$")
    if type(value).__name__ == "list":
        return all((cmplx.match(x) is not None for x in value))
    else:
        return cmplx.match(value) is not None


def _type_check_scalar(value: ScalarTypes, expected_type: str) -> bool:
    return type(value).__name__ == expected_type


def _type_check_list(value: ListTypes, expected_type: str) -> bool:
    if type(value).__name__ == "list":
        type_checks = all((_type_check_scalar(x, expected_type) for x in value))
    else:
        type_checks = False

    return type_checks


type_fixers = {
    "bool": bool,
    "complex": complex,
    "float": float,
    "int": int,
    "str": str,
    "List[bool]": lambda x: list(map(bool, x)),
    "List[complex]": lambda x: list(map(complex, x)),
    "List[float]": lambda x: list(map(float, x)),
    "List[int]": lambda x: list(map(int, x)),
    "List[str]": lambda x: list(map(str, x)),
}  # type: Dict[str, Callable[[Any], Any]]
"""Dict[str, Callable[[Any], Any]]: dictionary holding functions for type fixation."""
