#!/usr/bin/python

import re
import glob
import os
import subprocess


PDBs = []
STRs = []
log = open('consensx.log', 'w')

# collect PDB files into PDBs
for file in glob.glob("*_fit.pdb"):
    PDBs.append(file)

# collect SRT files into STRs
for file in glob.glob("STRs/*.str"):
    STRs.append(file)

print '{0:>2} {1:>8}'.format(len(PDBs), "fitted PDB file(s) found")
print '{0:>2} {1:>8}'.format(len(STRs), "STR file(s) found")


# ==========================================================================
# consensx_s2rdc.pl -b <STR_file> -f <PDB_file>
#             [-l <RDC_LC_MODEL: bic|pf1>] [-R] [-O<txt|html>] [-o<outfile>]
#             [-T<linethickness_in_graphs>]
#       -R causes to do SVD for back-calculating RDC  data
# ==========================================================================

for pop in PDBs:
    for str_data in STRs:
        print("PROCESSING: " + pop + " WITH " + str_data)

        # set output file name from pop and str names
        output_file = (pop.split('.')[0] + "_" +
                                     str_data.split('/')[1].split('.')[0])
        # calling ConsensX
        subprocess.call(["consensx",
                        "-b", str_data,    # restraints (.srt)
                        "-f", pop,         # population (.pdb)
                        "-R",              # back-calculating with SVD
                        "-O html",         # set plain text output
                        "-o", output_file, # set output file name
                        "-T", "0",         # do not generate graphs
                        ], stdout=log, stderr=log)

# clean temp folder
to_delete = []

for file in glob.glob("temp/*.pdb"):
    to_delete.append(file)

for file in glob.glob("temp/*_list"):
    to_delete.append(file)

for file in glob.glob("temp/*dummy*"):
    to_delete.append(file)

for item in to_delete:
    os.remove(item)

# process HTMLs and generate CSV output
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


class Csv_line:
    def __init__(self, pop, exp):
        self.pop  = pop
        self.exp  = exp
        self.s2   = {}
        self.rdcs = []

    def write_line(self):
        return_line = (self.pop + ';' + self.exp + ';' +
                       self.s2["corr"] + ';' + self.s2["q_val"] + ';')

        my_rdcs = []

        for rdc in self.rdcs:
            my_rdcs.append(rdc["bond"])

        if "0_N_H" in my_rdcs:
            for rdc in self.rdcs:
                if rdc["bond"] == "0_N_H":
                    return_line += rdc["corr"] + ';' + rdc["q_val"] + ';'
        else:
            return_line += ';;'

        if "0_CA_HA" in my_rdcs:
            for rdc in self.rdcs:
                if rdc["bond"] == "0_CA_HA":
                    return_line += rdc["corr"] + ';' + rdc["q_val"] + ';'
        else:
            return_line += ';;'

        if "0_CA_C" in my_rdcs:
            for rdc in self.rdcs:
                if rdc["bond"] == "0_CA_C":
                    return_line += rdc["corr"] + ';' + rdc["q_val"] + ';'
        else:
            return_line += ';;'

        if "0_C_N" in my_rdcs:
            for rdc in self.rdcs:
                if rdc["bond"] == "0_C_N":
                    return_line += rdc["corr"] + ';' + rdc["q_val"] + ';'
        else:
            return_line += ';;'

        if "0_C_H" in my_rdcs:
            for rdc in self.rdcs:
                if rdc["bond"] == "0_C_H":
                    return_line += rdc["corr"] + ';' + rdc["q_val"] + ';'
        else:
            return_line += ';;'

        return return_line + "\n"


HTMLs = glob.glob("*.html")
HTMLs.sort(key = natural_keys)
Csv_line_buffer = []

for html in HTMLs:
    pop = html.split('_')[1]
    exp = html.split('_')[-1].split('.')[0]

    html_file = open(html)
    my_line   = Csv_line(pop, exp)

    for line in html_file:

        if "2</sup> N" in line:

            s2_corr   = line.split()[3]
            s2_q_val  = line.split()[5]

            my_line.s2 = {"corr" : s2_corr, "q_val" : s2_q_val}

        elif "RDC " in line:

            bond_type = line.split()[2].split('<')[0]
            corr_val  = line.split()[3]
            q_value   = line.split()[5]

            my_line.rdcs.append({
                                "bond"  : bond_type,
                                "corr"  : corr_val,
                                "q_val" : q_value
                                })

    Csv_line_buffer.append(my_line)

csv_file = open("results.csv", 'w')

header   = ("POP;EXP;S2(0_N_H) CORR;Q_VAL;"
            "RDC(0_N_H) CORR;Q_VAL;RDC(0_CA_HA) CORR;Q_VAL;"
            "RDC(0_CA_C) CORR;Q_VAL;RDC(0_C_N) CORR;Q_VAL;"
            "RDC(0_C_H) CORR;Q_VAL\n")

csv_file.write(header)

for line in Csv_line_buffer:
    csv_file.write(line.write_line())
