#!/usr/bin/python
# -*- coding: utf-8 -*-

### Creates a PALES dummy file based on the PDB structure supplied.


import argparse
import subprocess
import os


pales   = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"
rdc_lc_model = "bic"

shortcodes = {
    'ALA':'A',  'ASP':'D',  'ASN':'N',  'ARG':'R',  'CYS':'C',  'GLY':'G',
    'GLU':'E',  'GLN':'Q',  'HIS':'H',  'ILE':'I',  'LEU':'L',  'LYS':'K',
    'MET':'M',  'PHE':'F',  'PRO':'P',  'SER':'S',  'THR':'T',  'TRP':'W',
    'TYR':'Y',  'VAL':'V'
}


#--------------------  Setting up and parsing arguments   --------------------#
usage  = '%(prog)s <pdb_file> pales_dummy_file'
parser = argparse.ArgumentParser(add_help = True,
                                 usage    = usage)
parser.add_argument("pdb_file", help="input PDB file")
args = parser.parse_args()


#-----------------------  Open file and read PDB data  -----------------------#
class AA(object):
    """Class for storing amino acid data"""
    def __init__(self, resname, resnum):
        self.resname = resname
        self.resnum  = resnum

try:
    input_pdb = open(args.pdb_file)     # open input PDB file or quit
except IOError:
    print("Input file " + args.pdb_file + " was not found")

seg = []                                # list for storing sequence data

for line in input_pdb:
    if line.startswith("ATOM") and line.split()[2] == "CA":
        resname = line.split()[3]       # get residue name
        resnum  = line.split()[5]       # get residue number
        seg.append(AA(resname, resnum))


#---------------------------  Write sequence data  ---------------------------#
short_seg = ""

for aa in seg:
    short_seg += shortcodes[aa.resname]

seg_lines = []
my_line = "DATA SEQUENCE "
char_counter = 0
row_counter  = 0

for char in short_seg:
    if char_counter == 10:              # write aa output in 10 wide blocks
        my_line += " "
        char_counter = 0
        row_counter += 1

        if row_counter == 5:            # write 5 block per line
            seg_lines.append(my_line + "\n")
            char_counter = 0
            row_counter  = 0
            my_line = "DATA SEQUENCE "

    my_line += char
    char_counter += 1

seg_lines.append(my_line + "\n")        # write last line of aa output

pales_dummy = open('pales_dummy.txt', 'w')

for line in seg_lines:
    pales_dummy.write(line)


#---------------------------  Write dummy dipoles  ---------------------------#
pales_dummy.write("""
VARS RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D DD W
FORMAT %5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f \n
""")

for i, aa in enumerate(seg):
    if aa.resname == "PRO":   # skip PRO named aa-s
        continue

    # print aligned dummy dipole output
    pales_dummy.write(" %5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
        str(i+1)+'A', aa.resname, 'H',
        str(i+1)+'A', aa.resname, 'N',
        0.000, 1.000,  1.00))
    pales_dummy.write(" %5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
        str(i+1)+'A', aa.resname, 'N',
        str(i+1)+'A', aa.resname, 'C',
        0.000, 1.000,  1.00))
    pales_dummy.write(" %5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
        str(i+1)+'A', aa.resname, 'C',
        str(i+1)+'A', aa.resname, 'CA',
        0.000, 1.000,  1.00))

pales_dummy.close()



# {`$PALES -inD $palesdummy -pdb $_ -$rdc_lc_model >> $rdcpalesout 2> /dev/null`;}
# {`$PALES -inD $palesdummy -pdb $_ -$rdc_lc_model -bestFit >> $rdcpalesout 2> /dev/null`;}

outfile = open("pales.out", 'a')
DEVNULL = open(os.devnull, 'w')

# http://stackoverflow.com/questions/16450788/python-running-subprocess-in-parallel

# processes = []
# for file in files_output:
#     f = os.tmpfile()
#     p = subprocess.Popen(['md5sum',file],stdout=f)
#     processes.append((p, f))

# for p, f in processes:
#     p.wait()
#     f.seek(0)
#     logfile.write(f.read())
#     f.close()



subprocess.call([pales,
                "-inD", "pales_dummy.txt",
                "-pdb", args.pdb_file,
                '-'+rdc_lc_model,
                "-bestFit"],
                stdout = outfile,
                stderr = DEVNULL)
