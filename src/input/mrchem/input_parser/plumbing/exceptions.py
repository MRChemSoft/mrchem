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

"""Error-handling facilities."""

from collections import namedtuple
from typing import List

from .utils import location_in_dict


class ParselglossyError(Exception):
    """Exception raised when parsing fails."""

    pass


class Error(namedtuple("Error", ["address", "message"])):
    """Detailed error reporting for dictionaries.

    Attributes
    ----------
    address : Tuple
         The keys needed to access the offending element in the `dict`.
    message : str
         The error message.

    Notes
    -----
    As we need to support Python 3.6, we cannot use the defaults field of
    ``namedtuple`` or ``dataclass``.
    """

    __slots__ = ()

    def __new__(cls, address=(), message=""):
        return super(Error, cls).__new__(cls, address, message)

    def __eq__(self, other):
        return self.address == other.address and self.message == other.message

    def __repr__(self):
        msg = f"{self.message:s}"
        if self.address != ():
            msg = f"At {location_in_dict(address=self.address):s}:\n  {msg:s}"
        return "- " + msg


def collate_errors(*, when: str, errors: List[Error]) -> str:
    """Collate a list of error into an informative message.

    Parameters
    ----------
    when: str
        When the error occurred.
    errors: List[Error]
        List of errors.

    Returns
    -------
    msg: str
        An error message with details about where in the dictionary the error
        arose and what the error is. For example::

            Error(s) occurred when fixing defaults:
            - At 'user['scf']['another_number']'
                KeyError 'min_num_iterations' in closure
                'user['scf']['min_num_iterations'] / 2'
    """
    plural = "s" if len(errors) > 1 else ""
    preamble = f"\nError{plural:s} occurred when {when:s}:"
    msgs = [preamble] + [f"{e}" for e in errors]

    return "\n".join(msgs)
