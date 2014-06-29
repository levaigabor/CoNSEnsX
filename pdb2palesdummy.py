#!/usr/bin/python
# -*- coding: utf-8 -*-

###################################################
# pdb2palesdummy.py <pdb_file> pales_dummy_file   #
#                                                 #
# Creates a PALES dummy file based on the PDB     #
# structure supplied.                             #
#                                                 #
# Version 0.1                         15.11.2004. #
###################################################

import argparse


shortcodes = {
    'ALA':'A',
    'ASP':'D',
    'ASN':'N',
    'ARG':'R',
    'CYS':'C',
    'GLY':'G',
    'GLU':'E',
    'GLN':'Q',
    'HIS':'H',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V'
}

#--------------------  Setting up and parsing arguments   --------------------#
usage  = '%(prog)s <pdb_file> pales_dummy_file'
parser = argparse.ArgumentParser(add_help = True,
                                 usage    = usage)
parser.add_argument("pdb_file", help="input PDB file")
args = parser.parse_args()


#-----------------------  Open file and read PDB data  -----------------------#
try:
    input_pdb = open(args.pdb_file)
except IOError:
    print("Input file " + args.pdb_file + " was not found")

seg = []    # list for storing sequence data

for line in input_pdb:
    if line.startswith("ATOM") and line.split()[2] == "CA":
        resname = line.split()[3]
        resnum  = line.split()[5]
        seg.append([resname, resnum])


#---------------------------  Write sequence data  ---------------------------#
short_seg = ""

for _ in seg:
    short_seg += shortcodes[_[0]]

seg_lines = []
my_line = "DATA SEQUENCE "
char_counter = 0
row_counter  = 0

for char in short_seg:
    if char_counter == 10:
        my_line += " "
        char_counter = 0
        row_counter += 1

        if row_counter == 5:
            seg_lines.append(my_line)
            char_counter = 0
            row_counter  = 0
            my_line = "DATA SEQUENCE "


    my_line += char
    char_counter += 1
seg_lines.append(my_line)

for line in seg_lines:
    print(line)


#---------------------------  Write dummy dipoles  ---------------------------#
print("""
VARS RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D DD W
FORMAT %5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f
""")


    # 1A     MET       H     1A     MET       N      0.000      1.000  1.00
    # 1A     MET       N     1A     MET       C      0.000      1.000  1.00

for i, _ in enumerate(seg):
    if _[0] == "PRO":
        continue

    print " %5s  %6s  %6s  %5s  %6s  %6s      0.000      1.000  1.00" % (
        str(i+1)+'A', _[0], 'H', str(i+1)+'A', _[0], 'N')
    print " %5s  %6s  %6s  %5s  %6s  %6s      0.000      1.000  1.00" % (
        str(i+1)+'A', _[0], 'N', str(i+1)+'A', _[0], 'C')
    print " %5s  %6s  %6s  %5s  %6s  %6s      0.000      1.000  1.00" % (
        str(i+1)+'A', _[0], 'C', str(i+1)+'A', _[0], 'CA')

