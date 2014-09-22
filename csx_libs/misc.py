#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse

#-----------------------  Help and credit information   ----------------------#
help_text = """
usage: cons_test.py -b <STR_file> -f <PDB_file>
             [-l <RDC_LC_MODEL: <bic|pf1>] [-R] [-O <txt|html>] [-o <outfile>]
             [-T <linethickness_in_graphs>]

  -h, --help            show this help message and exit
  -c                    show credits

  -b STR_file           restraints file
  -f PDB_file           PDB file (fitted)

  -l {bic,pf1}          rdc lc model: <bic|pf1>
  -R                    causes to do SVD for back-calculating RDC data
  -r                    X-PLOR restraint file
  -O {txt,html}         output format
  -o OUTPUT_FILE_NAME   output file name
  -T LINE_THICKNESS     line thickness in raphs


CoNSEnsX: assessing the compliance of varios NMR data with a protein
          structural ensemble
The program invokes SHIFTX for chemical shift and PALES for RDC calculation
See their respective documentation for more.
"""

cred = """
* Description of the CoNSEnsX concept and method can be found in:
-----------------------------------------------------------------
  CoNSEnsX: assessing the accuracy of NMR-derived protein structural ensembles
  Ángyán et al, 2009, submitted

* SHIFTX is described in:
------------------------
   Rapid and accurate calculation of protein 1H. 13C amd 15N
   Neal et al. (2003) J. Biomol. NMR 26:215.

* PALES is described in:
------------------------
  Prediction of sterically induced alignment in a dilute liquid
  crystalline phase: aid to protein strcuture determination by NMR
  Zweckstetter & Bax (2000) J. Am. Chem. Soc. 122:3791
"""

def createParser():
    """Create parser object for CoNSEnsX"""
    consensx_usage = '%(prog)s -b STR_file -f PDB_file [options]'

    parser = argparse.ArgumentParser(add_help = False,
                                    usage    = consensx_usage)

    parser.add_argument("-b", "--STR_file", help = "restraints file")
    parser.add_argument("-f", "--PDB_file", help = "PDB file (fitted)")
    parser.add_argument("-r", "--XPLOR_file", default="",
                       help = "X-PLOR restraint file")

    parser.add_argument("-h", "--help", help="show help and exit",
                       action='store_true')
    parser.add_argument('-c', action='store_true',
                       help="show credits")

    parser.add_argument("-l", "--lc_model", choices=["bic", "pf1"],
                       default="bic", help="rdc lc model: <bic|pf1>")
    parser.add_argument("-R", action='store_true', default=False,
                       help="causes to do SVD for back-calculating RDC data")
    parser.add_argument("-O", "--output_format", choices=["txt", "html"],
                       default="html|txt", help="output format")
    parser.add_argument("-o", "--output_file_name", default="consensx",
                       help="output file name")
    parser.add_argument("-T", "--line_thickness", help="raphs line thickness")

    return parser       # return created parser object
