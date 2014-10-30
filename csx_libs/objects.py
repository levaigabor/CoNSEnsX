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


class JCoup_Record(object):
    """Class for storing J-Coupling data"""
    def __init__(self, resnum, jcoup_type, JCoup_value):
        self.resnum = int(resnum)
        self.type   = jcoup_type
        self.value  = float(JCoup_value)


class ChemShift_Record(object):
    """Class for storing chemical shift data"""
    def __init__(self, resnum, res_name, atom_name, ChemShift_value):
        self.resnum    = int(resnum)
        self.res_name  = res_name
        self.atom_name = atom_name
        self.value     = float(ChemShift_value)


class Restraint_Record(object):
    """Class for storing restraint data"""
    def __init__(self, curr_distID, seq_ID1, seq_ID2, atom_ID1, atom_ID2, dist_max):
        self.curr_distID = int(curr_distID)
        self.seq_ID1     = int(seq_ID1)
        self.seq_ID2     = int(seq_ID2)
        self.atom_ID1    = str(atom_ID1)
        self.atom_ID2    = str(atom_ID2)
        self.dist_max    = float(dist_max)


class Vec_3D(object):
    """Vector class for calculations"""
    def __init__(self, v):
        self.v = v

    def __getitem__(self, index):
        return self.v[index]

    def __len__(self):
        return len(self.v)

    def __sub__(self, other):
        v = self.v
        u = other.v
        return Vec_3D([u[i] - v[i] for i in range(len(u))])

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

    def cross(one, other):
        c = [one.v[1] * other.v[2] - one.v[2] * other.v[1],
             one.v[2] * other.v[0] - one.v[0] * other.v[2],
             one.v[0] * other.v[1] - one.v[1] * other.v[0]]

        return Vec_3D(c)

    def dihedAngle(one, other):
        return math.degrees(math.acos((one.v[0] * other.v[0] +
                                       one.v[1] * other.v[1] +
                                       one.v[2] * other.v[2]) /
                                      (one.magnitude() * other.magnitude())))


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
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])
