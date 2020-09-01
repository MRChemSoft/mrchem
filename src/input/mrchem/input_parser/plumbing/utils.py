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

"""Common utilities."""

import json
from functools import reduce
from pathlib import Path
from shutil import copy
from typing import Any, Dict, List, Optional, Tuple, Union

JSONDict = Dict[str, Any]


class ComplexEncoder(json.JSONEncoder):
    """JSON encoder for complex numbers."""

    def default(self, obj):
        if isinstance(obj, complex):
            return {"__complex__": [obj.real, obj.imag]}
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)


def as_complex(dct):
    """JSON decoder for complex numbers."""
    if "__complex__" in dct:
        return complex(dct["__complex__"][0], dct["__complex__"][1])
    return dct


def location_in_dict(*, address: Tuple, dict_name: str = "user") -> str:
    """Convert tuple of keys of a ``JSONDict`` to its representation in code.

    For example, given ``("a", "b", "c")`` returns the string ``user['a']['b']['c']``.

    Parameters
    ----------
    address : Tuple[str]
    dict_name : str

    Returns
    -------
    where : str
    """
    return reduce(lambda x, y: x + f"['{y}']", address, "user")


def path_resolver(f: Union[str, Path], *, touch: bool = True) -> Path:
    """Resolve a path.

    Parameters
    ----------
    f : Union[str, Path]
        File whose path needs to be resolved.
    touch : bool
        Create file is not already existent. Default: True

    Returns
    -------
    path : Path
        File as a ``Path`` object.

    Notes
    -----
    The file will be created if not already existent.
    """

    path = Path(f) if isinstance(f, str) else f

    if not path.exists() and touch:
        path.touch()

    return path.resolve()


def default_outfile(*, fname: Union[str, Path], suffix: str) -> str:
    """Default name for output file.

    Parameters
    ----------
    fname : Union[str, Path]
        Name to use as stencil.
    suffix : str
        Suffix to append.

    Returns
    -------
    The name of the output file.
    """
    fname = Path(fname) if isinstance(fname, str) else fname

    base = fname.name

    return base.rsplit(".", 1)[0] + suffix


def copier(src: Path, dest: Path) -> None:
    """Copy file, ensuring it can be overwritten."""
    # copy file
    _ = copy(src, dest, follow_symlinks=True)
    # Ensure we can overwrite, by doing a chmod uga+rw
    Path(_).chmod(0o666)


def flatten_list(nested_list: List[Any]) -> List[Any]:
    """Flattens a nested list into a flat list.

    Parameters
    ----------
    nested_list : List[Any]
         Nested list

    Returns
    -------
    List[Any]
         Flattened list
    """
    if nested_list == []:
        return nested_list
    if isinstance(nested_list[0], list):
        return flatten_list(nested_list[0]) + flatten_list(nested_list[1:])
    return nested_list[:1] + flatten_list(nested_list[1:])


def dict_to_list(d: JSONDict) -> List[Any]:
    """Converts a nested dictionary to a nested list.

    Parameters
    ----------
    d : JSONDict
         Nested dictionary

    Returns
    -------
    list: List[Any]
         Flattened list
    """
    nested_list = []
    for k, v in d.items():
        if isinstance(v, dict):
            nested_list.append([k, dict_to_list(v)])
        else:
            nested_list.append([k, v])
    return nested_list


def nested_get(d: JSONDict, *ks: str) -> Optional[Any]:
    """Get value from a nested dictionary.

    Parameters
    ----------
    d : JSONDict
    ks : str

    Returns
    -------
    v : Optional[Any]

    Notes
    -----
    Adapted from: https://stackoverflow.com/a/40675868/2528668
    """

    def _func(x: JSONDict, k: str) -> Optional[JSONDict]:
        return x.get(k, None) if isinstance(x, dict) else None

    return reduce(_func, ks, d)  # type: ignore


def nested_set(d: JSONDict, ks: Tuple[Any, ...], v: Any) -> None:
    """Set value in nested dictionary.

    Parameters
    ----------
    d : JSONDict
    ks : Tuple[str]
    v : Any

    Notes
    -----
    Adapted from: https://stackoverflow.com/a/13688108/2528668
    """
    for k in ks[:-1]:
        d = d.setdefault(k, {})
    d[ks[-1]] = v
