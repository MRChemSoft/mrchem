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

"""Tools to extract views of dictionaries."""

from typing import Any, Callable, List, Optional, Type

from .exceptions import ParselglossyError
from .utils import JSONDict


def view_by_type(d: JSONDict) -> JSONDict:
    """Partial application of :func:`view_by` for types.

    Parameters
    ----------
    d: JSONDict

    Returns
    -------
    outgoing: JSONDict
       A dictionary with a view by types.
    """
    return view_by("type", d, missing=ParselglossyError)


def view_by_default(d: JSONDict) -> JSONDict:
    """Partial application of :func:`view_by` for defaults.

    Parameters
    ----------
    d: JSONDict

    Returns
    -------
    outgoing: JSONDict
       A dictionary with a view by defaults.
    """
    return view_by("default", d)


def view_by_default_keywords(keywords: List[JSONDict]) -> JSONDict:
    """View by defaults only for lists of keywords.

    Parameters
    ----------
    keywords: List[JSONDict]

    Returns
    -------
    outgoing: JSONDict
       A dictionary with a view by defaults for the keywords.
    """
    return {v["name"]: v["default"] if "default" in v else None for v in keywords}


def view_by_docstring(d: JSONDict) -> JSONDict:
    """Partial application of :func:`view_by` for docstrings.

    Parameters
    ----------
    d: JSONDict

    Returns
    -------
    outgoing: JSONDict
       A dictionary with a view by docstrings.
    """

    def docstring_not_empty(x: Any, y: str) -> bool:
        """Check that a docstring is not empty."""
        return x[y].strip() != ""

    def docstring_rstrip(x: str) -> str:
        """Apply rstrip to a docstring"""
        return x.rstrip()

    return view_by(
        "docstring",
        d,
        predicate=docstring_not_empty,
        missing=ParselglossyError,
        transformer=docstring_rstrip,
    )


def view_by_predicates(d: JSONDict) -> JSONDict:
    """Partial application of :func:`view_by` for predicates.

    Parameters
    ----------
    d: JSONDict

    Returns
    -------
    outgoing: JSONDict
       A dictionary with a view by predicates.
    """
    return view_by("predicates", d)


def view_by(
    what: str,
    d: JSONDict,
    *,
    predicate: Callable[[Any, str], bool] = lambda x, y: True,
    missing: Optional[Type[Exception]] = None,
    transformer: Callable[[Any], Any] = lambda x: x,
) -> JSONDict:
    """Recursive decimation of a template into an input.

    Parameters
    ----------
    what : str
        What view to extract from the dictionary. Any of ``type``, ``default``,
        ``docstring`` or ``predicates`` is allowed.
    d: JSONDict
    predicate: Callable
       A predicate accepting two arguments.
    missing: Optional[Exception]
    transformer: Callable

    Returns
    -------
    outgoing: JSONDict
       A dictionary with the desired view.

    Raises
    ------
    exc:`ValueError` if `what` is not among the allowed views.
    """

    allowed = ["type", "default", "docstring", "predicates"]
    if what not in allowed:
        raise ValueError(
            f"Requested view {what:s} not among possible views ({allowed})"
        )

    # Determine whether we have only sections and only keywords
    has_keywords = True if "keywords" in d else False
    has_sections = True if "sections" in d else False

    view = {}  # type: JSONDict
    if has_keywords:
        view = {
            v["name"]: transformer(v[what])
            if all([what in v, predicate(v, what)])
            else missing
            for v in d["keywords"]
        }

    if has_sections:
        for section in d["sections"]:
            view[section["name"]] = view_by(
                what,
                section,
                predicate=predicate,
                missing=missing,
                transformer=transformer,
            )

    return view
