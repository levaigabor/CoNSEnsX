#!/usr/bin/python

# TODO -> format (align) output

import sys
import os


####        hadling parameters      ####
if len(sys.argv) < 2:
    sys.exit("Usage: top2strs.py <topology file> <output prefix>")

if sys.argv[1] == "-h":
    print("Split RDC data from topology file to seperate str files\n" +
          "Usage: top2strs.py <topology file> <output prefix>")
    sys.exit(0)

try:
    top_file        = open(sys.argv[1], 'r')
    output_prefix   = sys.argv[2]
except IndexError:
    sys.exit("Usage: top2strs.py <topology file> <output prefix>")
except IOError:
    sys.exit("Topology file was not found")


####        initialize globals      ####
rdc_sets   = []
rdc_buffer = []
in_block   = False


####        processing RDC data     ####
for line in top_file:
    line = line.strip()
    if line.endswith("[ orientation_restraints ]"):

        in_block = True
        continue

    if not line:
        in_block = False

        if len(rdc_buffer) > 0:
            rdc_sets.append(rdc_buffer)
            rdc_buffer = []

        continue

    if in_block and not line.startswith(";"):
        data = line.split()
        rdc_buffer.append([data[10], data[11][0:-1], data[13], data[7]])

print(str(len(rdc_sets)) + " RDC sets were found")


####        set output folder       ####
if not os.path.exists("STRs"):
    os.makedirs("STRs")


####        writing RDC data        ####
for i, rdc_set in enumerate(rdc_sets):

    output_name = "STRs/" + output_prefix + str(i + 1) + ".str"
    output_file = open(output_name, 'w')
    str_master  = open("master.str", 'r')

    for line in str_master:
        if line.strip() == "#  BITE ME PYTHON  #":
            for i in rdc_set:
                if i[1] == 'N' and i[2] == 'H':
                    short = "DNH"
                else:
                    short = i[1] + i[2]

                rdc_line = ("      " + short + " " + i[0] + " XXX " + i[1] + " "
                            + i[0] + " XXX " + i[2] + " " + i[3] + "  ? ? . . ?\n")

                output_file.write(rdc_line)
        else:
            output_file.write(line)

    output_file.close()
    str_master.close()

