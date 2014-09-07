#!/usr/bin/python
# -*- coding: utf-8 -*-



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
