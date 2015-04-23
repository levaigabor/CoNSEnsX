import math
import os
import copy

class RDC_Record(object):
    """Class for storing RDC data"""
    def __init__(self, resnum1, atom1, resnum2, atom2, RDC_value):
        self.RDC_type = (str(int(resnum1) - int(resnum2))
                         + '_' + atom1 + '_' + atom2)
        self.resnum   = int(resnum1)
        self.atom     = atom1
        self.resnum2  = int(resnum2)
        self.atom2    = atom2
        self.value    = float(RDC_value)


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


class PDB_model(object):
    """Class for storing PDB model data"""
    model_list = []

    def __init__(self, atomgroup,):
        self.atomgroup = atomgroup
        PDB_model.model_list.append(self.atomgroup)


class Restraint_Record(object):
    """Class for storing NOE restraint data"""
    all_restraints   = []

    def __init__(self, curr_distID, seq_ID1, seq_ID2, seq_name1, seq_name2,
                 atom_ID1, atom_ID2, dist_max):

        self.curr_distID = int(curr_distID)
        self.seq_ID1     = int(seq_ID1)
        self.seq_ID2     = int(seq_ID2)
        self.seq_name1   = str(seq_name1)
        self.seq_name2   = str(seq_name2)
        self.dist_max    = float(dist_max)
        resol = {
                "MET" : {
                            "ME"  : ["HE1",  "HE2",  "HE3"],
                            "MD"  : ["HD11", "HD12", "HD13"],
                            "MG"  : ["HG11", "HG12", "HG13"],
                            "MD1" : ["HD11", "HD12", "HD13"],
                            "MD2" : ["HD21", "HD22", "HD23"],
                            "MG1" : ["HG11", "HG12", "HG13"],
                            "MG2" : ["HG21", "HG22", "HG23"]
                        },
                "ILE" : {
                            "MD"  : ["HD11", "HD12", "HD13"],
                            "MG"  : ["HG21", "HG22", "HG23"],
                            "MB"  : ["HB1",  "HB2",  "HB3"],
                            "MD1" : ["HD11", "HD12", "HD13"],
                            "MG2" : ["HG21", "HG22", "HG23"]
                        },
                "LYS" : {
                            "MB"  : ["HB1",  "HB2",  "HB3"],
                            "MD1" : ["HD11", "HD12", "HD13"],
                            "MD2" : ["HD21", "HD22", "HD23"],
                            "MD"  : ["HD11", "HD12", "HD13"],
                            "MG2" : ["HG21", "HG22", "HG23"]
                        },
                "ALA" : {
                            "MD"  : ["HD11", "HD12", "HD13"],
                            "MB"  : ["HB1",  "HB2",  "HB3"],
                            "MG1" : ["HG11", "HG12", "HG13"],
                            "MG2" : ["HG21", "HG22", "HG23"],
                            "MD1" : ["HD11", "HD12", "HD13"],
                            "MD2" : ["HD21", "HD22", "HD23"]
                        },
                "LEU" : {
                            "MD1" : ["HD11", "HD12", "HD13"],
                            "MD2" : ["HD21", "HD22", "HD23"]
                        },
                "VAL" : {
                            "MG1" : ["HG11", "HG12", "HG13"],
                            "MG2" : ["HG21", "HG22", "HG23"]
                        },
                "GLN" : {
                            "MG"  : ["HG11", "HG12", "HG13"],
                            "MG1" : ["HG11", "HG12", "HG13"],
                            "MG2" : ["HG21", "HG22", "HG23"]
                        },
                "HIS" : {
                            "MG"  : ["HG11", "HG12", "HG13"],
                            "MG1" : ["HG11", "HG12", "HG13"]
                        },
                "PHE" : {
                            "MG2" : ["HG21", "HG22", "HG23"]
                        },
                "THR" : {
                            "MG"  : ["HG21", "HG22", "HG23"]
                        }
        }

        if atom_ID1.startswith('M'):
            atom_list1 = resol[seq_name1][atom_ID1]
        else:
            atom_list1 = [atom_ID1]

        if atom_ID2.startswith('M'):
            atom_list2 = resol[seq_name2][atom_ID2]
        else:
            atom_list2 = [atom_ID2]

        for atom1 in atom_list1:
            for atom2 in atom_list2:
                self.atom_ID1 = atom1
                self.atom_ID2 = atom2
                me = copy.deepcopy(self)
                Restraint_Record.all_restraints.append(me)

    @staticmethod
    def getRestraintCount():
        return Restraint_Record.all_restraints[-1].curr_distID

    @staticmethod
    def getPRIDE_restraints():
        PRIDE_restraints = {}
        prev_id = -1
        seq1_ok, seq2_ok, distance_ok = False, False, False

        for restraint in Restraint_Record.all_restraints:
            curr_id  = restraint.curr_distID
            atom_ID1 = restraint.atom_ID1
            atom_ID2 = restraint.atom_ID2
            seq_ID1  = restraint.seq_ID1
            seq_ID2  = restraint.seq_ID2

            if prev_id != curr_id:
                prev_id = curr_id

                if seq1_ok and seq2_ok and distance_ok:
                    if seq_dist in PRIDE_restraints:
                        PRIDE_restraints[seq_dist] += 1
                    else:
                        PRIDE_restraints[seq_dist] = 1

                seq1_ok  = atom_ID1 in ["H", "HA"] or atom_ID1.startswith("HB")
                seq2_ok  = atom_ID2 in ["H", "HA"] or atom_ID2.startswith("HB")
                seq_dist = abs(seq_ID1 - seq_ID2)
                distance_ok = seq_dist > 2

                id_distance = seq_dist

            else:
                seq1_ok  &= atom_ID1 in ["H", "HA"] or atom_ID1.startswith("HB")
                seq2_ok  &= atom_ID2 in ["H", "HA"] or atom_ID2.startswith("HB")
                seq_dist = abs(seq_ID1 - seq_ID2)
                distance_ok &= seq_dist > 2 and id_distance == seq_dist

        return PRIDE_restraints


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

    @classmethod
    def cross(cls, one, other):
        c = [one.v[1] * other.v[2] - one.v[2] * other.v[1],
             one.v[2] * other.v[0] - one.v[0] * other.v[2],
             one.v[0] * other.v[1] - one.v[1] * other.v[0]]

        return cls(c)

    @staticmethod
    def dihedAngle(one, other):
        return math.degrees(math.acos((one.v[0] * other.v[0] +
                                       one.v[1] * other.v[1] +
                                       one.v[2] * other.v[2]) /
                                      (one.magnitude() * other.magnitude())))
