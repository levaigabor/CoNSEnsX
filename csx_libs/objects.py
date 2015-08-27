import math
import copy


class RDC_modell_corr(object):

    """Class for per model RDC correlation"""
    RDC_lists = []

    def __init__(self, RDC_list_permodel):
        RDC_modell_corr.RDC_lists.append(RDC_list_permodel)

    @staticmethod
    def get_best_model(RDC_list_num):
        best_value = -1
        best_num   = -1

        for num, value in enumerate(RDC_modell_corr.RDC_lists[RDC_list_num]):
            if value > best_value:
                best_value = value
                best_num   = num

        return best_num, best_value


class ThirdParty(object):

    """Class to store 3rd party software information"""
    pales    = ""
    shiftx   = ""
    prideDB  = ""
    prideNMR = ""

    @staticmethod
    def get_thirdparty(config_file):
        cfg = open(config_file)
        for line in cfg:
            if line.startswith("#"):
                continue
            if line.startswith("pales"):
                ThirdParty.pales = line.split("'")[1]
            elif line.startswith("shiftx"):
                ThirdParty.shiftx = line.split("'")[1]
            elif line.startswith("prideDB"):
                ThirdParty.prideDB = line.split("'")[1]
            elif line.startswith("prideNMR"):
                ThirdParty.prideNMR = line.split("'")[1]


class CSV_buffer(object):

    """Class which stores data for values.CSV"""
    working_dir = ""
    max_resnum = -1
    min_resnum = 100000
    csv_data = []

    def __init__(self, name, calced, experimental):
        self.name = name
        self.calced = calced
        self.exp = {}

        for i in experimental:
            self.exp[i.resnum] = i.value

        CSV_buffer.csv_data.append(self)

    @staticmethod
    def writeCSV():
        filename = CSV_buffer.working_dir + "values.csv"
        # filename = "values.csv"
        output_csv = open(filename, 'w')
        output_csv.write(',')
        for data in CSV_buffer.csv_data:
            output_csv.write(data.name + " EXP, " + data.name + " CALC,")
        output_csv.write("\n")
        for resnum in range(CSV_buffer.min_resnum, CSV_buffer.max_resnum):
            output_csv.write(str(resnum) + ',')
            for data in CSV_buffer.csv_data:
                try:
                    output_csv.write("{0:.2f}".format(data.exp[resnum]) + ',' +
                                     "{0:.2f}".format(data.calced[resnum])
                                     + ',')
                except KeyError:
                    output_csv.write('')

            output_csv.write("\n")


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
        self.calced = None


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
    is_fitted  = False

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
            "MET": {
                "ME":  ["HE1",  "HE2",  "HE3"],
                "MD":  ["HD11", "HD12", "HD13"],
                "MG":  ["HG11", "HG12", "HG13"],
                "MD1": ["HD11", "HD12", "HD13"],
                "MD2": ["HD21", "HD22", "HD23"],
                "MG1": ["HG11", "HG12", "HG13"],
                "MG2": ["HG21", "HG22", "HG23"]
            },
            "ILE": {
                "MD":  ["HD11", "HD12", "HD13"],
                "MG":  ["HG21", "HG22", "HG23"],
                "MB":  ["HB1",  "HB2",  "HB3"],
                "MD1": ["HD11", "HD12", "HD13"],
                "MG2": ["HG21", "HG22", "HG23"]
            },
            "LYS": {
                "MB":  ["HB1",  "HB2",  "HB3"],
                "MD1": ["HD11", "HD12", "HD13"],
                "MD2": ["HD21", "HD22", "HD23"],
                "MD":  ["HD11", "HD12", "HD13"],
                "MG2": ["HG21", "HG22", "HG23"]
            },
            "ALA": {
                "MD":  ["HD11", "HD12", "HD13"],
                "MB":  ["HB1",  "HB2",  "HB3"],
                "MG1": ["HG11", "HG12", "HG13"],
                "MG2": ["HG21", "HG22", "HG23"],
                "MD1": ["HD11", "HD12", "HD13"],
                "MD2": ["HD21", "HD22", "HD23"]
            },
            "LEU": {
                "MD1": ["HD11", "HD12", "HD13"],
                "MD2": ["HD21", "HD22", "HD23"]
            },
            "VAL": {
                "MG1": ["HG11", "HG12", "HG13"],
                "MG2": ["HG21", "HG22", "HG23"]
            },
            "GLN": {
                "MG":  ["HG11", "HG12", "HG13"],
                "MG1": ["HG11", "HG12", "HG13"],
                "MG2": ["HG21", "HG22", "HG23"]
            },
            "HIS": {
                "MG":  ["HG11", "HG12", "HG13"],
                "MG1": ["HG11", "HG12", "HG13"]
            },
            "PHE": {
                "MG2": ["HG21", "HG22", "HG23"]
            },
            "THR": {
                "MG": ["HG21", "HG22", "HG23"]
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
        seq_dist = -1

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
