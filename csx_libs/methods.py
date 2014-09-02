#!/usr/bin/python
# -*- coding: utf-8 -*-


import subprocess
import os
import re


pales   = "/home/daniel/Dokumente/önlab/gz_pack/pales/linux/pales"
rdc_lc_model = "bic"

shortcodes = {
    'ALA':'A',  'ASP':'D',  'ASN':'N',  'ARG':'R',  'CYS':'C',  'GLY':'G',
    'GLU':'E',  'GLN':'Q',  'HIS':'H',  'ILE':'I',  'LEU':'L',  'LYS':'K',
    'MET':'M',  'PHE':'F',  'PRO':'P',  'SER':'S',  'THR':'T',  'TRP':'W',
    'TYR':'Y',  'VAL':'V'
}


class AA(object):
        """Class for storing amino acid data"""
        def __init__(self, resname, resnum):
            self.resname = resname
            self.resnum  = resnum


def callPalesOn(pdb_files):
    try:
        os.remove("pales.out")                  # remove output file if present
    except OSError:
        pass

    for pdb_file in pdb_files:
        #-------------------  Open file and read PDB data  -------------------#
        try:
            input_pdb = open(pdb_file)          # open input PDB file
        except IOError:
            print("Input file " + pdb_file + " was not found")
            raise SystemExit                    # exit if input file not found

        seg = []                                # list storing sequence data

        for line in input_pdb:
            if line.startswith("ATOM") and line.split()[2] == "CA":
                resname = line.split()[3]       # get residue name
                resnum  = line.split()[5]       # get residue number
                seg.append(AA(resname, resnum)) # append new objecjt to list


        #-----------------------  Write sequence data  -----------------------#
        short_seg = ""

        for aa in seg:
            short_seg += shortcodes[aa.resname]

        my_line      = "DATA SEQUENCE "
        char_counter = 0
        row_counter  = 0
        pales_dummy  = open('pales_dummy.txt', 'w')

        for char in short_seg:
            if char_counter == 10:          # write aa output in 10 wide blocks
                my_line      += " "
                char_counter = 0
                row_counter  += 1

                if row_counter == 5:        # write 5 block per line
                    pales_dummy.write(my_line + "\n")
                    char_counter = 0
                    row_counter  = 0
                    my_line = "DATA SEQUENCE "

            my_line      += char
            char_counter += 1

        pales_dummy.write(my_line + "\n")    # write last line of aa output


        #-----------------------  Write dummy dipoles  -----------------------#
        pales_dummy.write(
        "\nVARS RESID_I RESNAME_I ATOMNAME_I " +
        "RESID_J RESNAME_J ATOMNAME_J D DD W\n" +
        "FORMAT %5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f \n\n"
        )

        for i, aa in enumerate(seg):
            if aa.resname == "PRO":         # skip PRO named aa-s
                continue

            # print aligned dummy dipole output
            pales_dummy.write(
                "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
                str(i+1)+'A', aa.resname, 'H',
                str(i+1)+'A', aa.resname, 'N',
                0.000, 1.000,  1.00))
            pales_dummy.write(
                "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
                str(i+1)+'A', aa.resname, 'N',
                str(i+1)+'A', aa.resname, 'C',
                0.000, 1.000,  1.00))
            pales_dummy.write(
                "%5s  %6s  %6s  %5s  %6s  %6s  %9.3f  %9.3f  %.2f\n" % (
                str(i+1)+'A', aa.resname, 'C',
                str(i+1)+'A', aa.resname, 'CA',
                0.000, 1.000,  1.00))

        pales_dummy.close()

        outfile = open("pales.out", 'a')    # open output file with append mode
        DEVNULL = open(os.devnull, 'w')     # open systems /dev/null

        print("calculating: " + pdb_file)

        subprocess.call([pales,
                        "-inD", "pales_dummy.txt",  # pales dummy file
                        "-pdb", pdb_file,           # pdb file
                        '-' + rdc_lc_model,         # rdc lc model
                        "-bestFit"],
                        stdout = outfile,
                        stderr = DEVNULL)


def avgPalesRDCs(pales_out):
    pales_out       = open(pales_out)
    n_of_structures = 0
    averageRDC      = {}
    npair           = 0


    for line in pales_out:
        if re.match("REMARK \d+ couplings", line):
            n_of_structures += 1 # n_of_structures to divide by

        elif re.match("\s+ \d+", line):
            resnum1 = int(line.split()[0])
            resnum2 = int(line.split()[3])
            atom1   = line.split()[2]
            atom2   = line.split()[5]
            D       = float(line.split()[6])
            RDCtype = str(resnum2 - resnum1) + "_" + atom1 + "_" + atom2

            if RDCtype in averageRDC.keys():
                if resnum1 in averageRDC[RDCtype].keys():
                    averageRDC[RDCtype][resnum1] += D
                else:
                    averageRDC[RDCtype][resnum1] = D
            else:
                averageRDC[RDCtype] = {}
                averageRDC[RDCtype][resnum1] = D



    for RDCtype in averageRDC.keys():
        for resnum in averageRDC[RDCtype].keys():
            averageRDC[RDCtype][resnum] /= n_of_structures


    print averageRDC

