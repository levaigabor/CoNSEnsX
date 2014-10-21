import nmrpystar
import time
import math
import sys

ts = time.time()


class Restraint_Record(object):
    """Class for storing restraint data"""
    def __init__(self, curr_distID, seq_ID1, seq_ID2, atom_ID1, atom_ID2, dist_val):

        self.curr_distID = int(curr_distID)
        self.seq_ID1     = int(seq_ID1)
        self.seq_ID2     = int(seq_ID2)
        self.atom_ID1    = str(atom_ID1)
        self.atom_ID2    = str(atom_ID2)
        self.dist_val    = float(dist_val)


def parseSTR(STR_file):
    star_file = open(STR_file)
    myString = ""

    for line in star_file:
        myString += line

    parsed = nmrpystar.parse(myString)

    if parsed.status != 'success':
        print('Error during STR parsing: ', parsed)
        raise SystemExit
    else:
        return parsed


sys.argv[1]

my_parsed  = parseSTR(sys.argv[1])
saveShifts = my_parsed.value

saveShiftName = "CNS/XPLOR_distance_constraints_"
NOEkey        = "Gen_dist_constraint_list.Constraint_type"

for frame in saveShifts.saves.keys():
    if frame.startswith(saveShiftName):
        if NOEkey in saveShifts.saves[frame].datums.keys():
            if saveShifts.saves[frame].datums[NOEkey] == "NOE":
                loopShifts = saveShifts.saves[frame].loops[-1]
                break

prev_distID = 0

# dist_dict = {}
# dist_acc  = []

restraints = []

for ix in range(len(loopShifts.rows)):   # fetch values from STR file
    row = loopShifts.getRowAsDict(ix)

    curr_distID = row["Gen_dist_constraint.ID"]
    seq_ID1     = row["Gen_dist_constraint.Seq_ID_1"]
    seq_ID2     = row["Gen_dist_constraint.Seq_ID_2"]
    atom_ID1    = row["Gen_dist_constraint.Atom_ID_1"]
    atom_ID2    = row["Gen_dist_constraint.Atom_ID_2"]
    dist_val    = row["Gen_dist_constraint.Distance_val"]

    restraints.append(Restraint_Record(curr_distID, seq_ID1, seq_ID2,
                                       atom_ID1, atom_ID2, dist_val))


#     if curr_distID == prev_distID:
#         dist_acc.append(dist_val)

#     else:
#         if dist_acc:

#             avg_dist = 0.0

#             for distance in dist_acc:
#                 avg_dist += math.pow(float(distance), -6)

#             dist_dict[prev_distID] = math.pow(avg_dist, -1/6)


#         prev_distID = curr_distID
#         dist_acc    = [dist_val]

# print dist_dict


te = time.time()
print "runtime:", te - ts
