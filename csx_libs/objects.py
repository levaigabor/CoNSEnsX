#!/usr/bin/python
# -*- coding: utf-8 -*-

import math
import os

class RDC_Record(object):
    """Class for storing RDC data"""

    def __init__(self, resnum1, atom1, resnum2, atom2, RDC_value):
        self.RDC_type  = (str(int(resnum1) - int(resnum2))
                          + '_' + atom1 + '_' + atom2)
        self.resnum1   = int(resnum1)
        self.atom1     = atom1
        self.resnum2   = int(resnum2)
        self.atom2     = atom2
        self.value     = float(RDC_value)


class S2_Record(object):
    """Class for storing S2 data"""

    def __init__(self, resnum, S2_type, S2_value):
        self.resnum = int(resnum)
        self.type   = S2_type
        self.value  = float(S2_value)


class Vec_3D(object):
    """Vector class for S2 calculation"""

    def __init__(self, v):
        self.v = v

    def __sub__(self, other):
        v = self.v
        u = other.v
        return [u[i] - v[i] for i in range(len(u))]

    def __iadd__(self, other):
        v = self.v
        u = other.v
        return Vec_3D([u[i] + v[i] for i in range(len(u))])

    def __idiv__(self, divisor):
        v = self.v
        return Vec_3D([v[i] / divisor for i in range(len(v))])

    def magnitude(self):
        v = self.v
        return math.sqrt(sum(v[i] * v[i] for i in range(len(v))))

    def normalize(self):
        vmag = self.magnitude()
        v = self.v
        return Vec_3D([v[i] / vmag  for i in range(len(v))])


# Define a context manager to suppress stdout and stderr.
class suppress_output(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited.
    more info at -> http://stackoverflow.com/questions/11130156/
    suppress-stdout-stderr-print-from-python-functions

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])
