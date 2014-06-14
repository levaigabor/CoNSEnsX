#!/usr/bin/python3

import argparse

###############################################################################
# usage: cons_test.py [-h] [-l LC_MODEL] [-R] [-O OUTPUT_FORMAT]
#                   [-o OUTPUT_FILE_NAME] [-T LINE_THICKNESS]
#                   STR_file PDB_file
#
# CoNSEnsX: assessing the compliance of varios NMR data with a protein
#           structural ensemble
# The program invokes SHIFTX for chemical shift and PALES for RDC calculation
# See their respective documentation for more.
###############################################################################


parser = argparse.ArgumentParser()
parser.add_argument("STR_file", help="restraints file")
parser.add_argument("PDB_file", help="PDB file (fitted)")

parser.add_argument("-l", "--lc_model", choices=["bic", "pf1"], default="bic",
                    help="rdc lc model: <bic|pf1>")
parser.add_argument("-R", action='store_true', default=False,
                    help="causes to do SVD for back-calculating RDC data")
parser.add_argument("-O", "--output_format", choices=["txt", "html"],
                    default="html|txt", help="output format")
parser.add_argument("-o", "--output_file_name", default="consensx",
                    help="output file name")
parser.add_argument("-T", "--line_thickness",   help="line thickness in raphs")
args = parser.parse_args()

print(args.lc_model)
print(args.R)
print(args.output_format)
print(args.output_file_name)
print(args.line_thickness)
print(args.STR_file)
print(args.PDB_file)
