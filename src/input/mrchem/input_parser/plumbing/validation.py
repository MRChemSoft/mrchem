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

"""Validation facilities."""

import json
from pathlib import Path
from typing import Union

from .exceptions import ParselglossyError, collate_errors
from .utils import ComplexEncoder, JSONDict, path_resolver
from .validation_plumbing import (
    _rec_check_predicates,
    _rec_fix_defaults,
    _rec_merge_ours,
)
from .views import view_by_default, view_by_predicates, view_by_type


def validate_from_dicts(
    *, ir: JSONDict, template: JSONDict, fr_file: Union[str, Path] = None
) -> JSONDict:
    """Validate intermediate representation into final representation.

    Parameters
    ----------
    ir : JSONDict
        Intermediate representation of the input file.
    template : JSONDict
        A _validated_ template.
    fr_file : Union[str, Path]
         File to write final representation to (JSON format).
         None by default, which means file is not written out.

    Returns
    -------
    fr : JSONDict
        The validated input.

    Raises
    ------
    :exc:`ParselglossyError`
    """
    stencil = view_by_default(template)
    types = view_by_type(template)
    predicates = view_by_predicates(template)

    fr = merge_ours(theirs=stencil, ours=ir)
    fr = fix_defaults(fr, types=types)
    check_predicates(fr, predicates=predicates)

    if fr_file is not None:
        fr_file = path_resolver(fr_file)
        with fr_file.open("w") as out:
            json.dump(fr, out, cls=ComplexEncoder, indent=4)

    return fr


def merge_ours(*, theirs: JSONDict, ours: JSONDict) -> JSONDict:
    """Recursively merge two `dict`-s with "ours" strategy.

    Parameters
    ----------
    theirs : JSONDict
    ours : JSONDict

    Returns
    -------
    outgoing : JSONDict

    Raises
    ------
    :exc:`ParselglossyError`

    Notes
    -----
    The ``theirs`` dictionary is supposed to be the view by defaults of the
    validation specification, whereas ``ours`` is the dictionary from user
    input. The recursive merge action will generate a complete, but not
    validated, input dictionary by using default values where these are not
    overridden by user input, hence the naming "ours" for the merge strategy.

    This is porcelain over the recursive function :func:`_rec_merge_ours`.
    """
    outgoing, errors = _rec_merge_ours(theirs=theirs, ours=ours)

    if errors:
        msg = collate_errors(when="merging", errors=errors)
        raise ParselglossyError(msg)

    return outgoing


def fix_defaults(incoming: JSONDict, *, types: JSONDict) -> JSONDict:
    """Fix default values and perform type checking.

    Parameters
    ----------
    incoming: JSONDict
        The input ``dict``. This is supposed to be the one obtained by merging
        user and template ``dict``-s.
    types: JSONDict
        Types of all keywords in the input. Generated from :func:`view_by_types`.

    Returns
    -------
    outgoing: JSONDict
        A dictionary with all default values fixed.

    Raises
    ------
    :exc:`ParselglossyError`

    Notes
    -----
    Since we allow callables to appear as defaults, we need to run them
    to determine the actual default values.

    This operation must be done with some care, to avoid false negatives or
    ambiguous type checks. For example:

    * If the type is ``str`` and the default a callable, the type will
      match, but the default will make no sense.
    * If the type is numerical, *e.g.* ``int``, the type will not match.

    However, by design the callables **must** refer to some other field in the
    input tree, hence they **must** contain the reserved token "user". This
    allows us to disambiguate a callable as default from a value as default.

    The final strategy adopted is then:

    1. Perform type checking with :func:`type_matches`. If successful, we
    coerce the type.
    2. If types did not match, we further check whether the value is a string,
    containing the reserved token "value". This means the default value is
    actually a callable. We run the callable, which internally coerces the type
    of the result to the expected one.
    3. If even this check was unsuccessful, types really were unmatched. We
    report the error and move on.

    This is porcelain over recursive function :func:`_rec_fix_defaults`.
    """

    outgoing, errors = _rec_fix_defaults(incoming, types=types)

    if errors:
        msg = collate_errors(when="fixing defaults", errors=errors)
        raise ParselglossyError(msg)

    return outgoing


def check_predicates(incoming: JSONDict, *, predicates: JSONDict) -> None:
    """Run predicates on input tree with fixed defaults.

    Parameters
    ----------
    incoming : JSONDict
        The input ``dict``. This is supposed to be the result of :func:`fix_defaults`.
    predicates : JSONDict
        A view-by-predicates of the template ``dict``.

    Raises
    ------
    :exc:`ParselglossyError`

    Notes
    -----
    This is porcelain over recursive function :func:`_rec_check_predicates`.
    """

    errors = _rec_check_predicates(incoming, predicates=predicates)

    if errors:
        msg = collate_errors(when="checking predicates", errors=errors)
        raise ParselglossyError(msg)
