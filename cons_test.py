#!/usr/bin/python3

import argparse

###############################################################################
# usage: cons_test.py [-h] STR_file PDB_file
#        [-l {bic,pf1}] [-R] [-O {txt,html}]
#        [-o OUTPUT_FILE_NAME] [-T LINE_THICKNESS]
#
# CoNSEnsX: assessing the compliance of varios NMR data with a protein
#           structural ensemble
# The program invokes SHIFTX for chemical shift and PALES for RDC calculation
# See their respective documentation for more.
###############################################################################


class RawUsageHelpFormatter(argparse.HelpFormatter):
    def _format_usage(self, usage, actions, groups, prefix):
        # use actions in the order that they are define
        # no line wrapping
        if prefix is None:
            prefix = 'usage: '
        if usage is not None:
            usage = usage % dict(prog=self._prog)
        elif usage is None and not actions:
            usage = '%(prog)s' % dict(prog=self._prog)
        elif usage is None:
            prog = '%(prog)s' % dict(prog=self._prog)
            format = self._format_actions_usage
            action_usage = format(actions, groups)
            usage = ' '.join([s for s in [prog, action_usage] if s])
        return '%s%s\n\n' % (prefix, usage)


parser = argparse.ArgumentParser(formatter_class=RawUsageHelpFormatter)
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
