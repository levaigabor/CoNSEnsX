#!/usr/bin/python
# -*- coding: utf-8 -*-

import math

class RDC_Record(object):
    """Class for storing RDC data"""

    def __init__(self, resnum1, atom1, resnum2, atom2, RDC_value):
        self.RDC_type  = (str(int(resnum1) - int(resnum2))
                          + '_' + atom1 + '_' + atom2)
        self.resnum1   = int(resnum1)
        self.atom1     = atom1
        self.resnum2   = int(resnum2)
        self.atom2     = atom2
        self.RDC_value = float(RDC_value)


class Vec_3D(object):

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
