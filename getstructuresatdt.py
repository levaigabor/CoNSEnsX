#!/usr/bin/python

import os

program_suffix  = "_455s_mumo"
molname         = "1d3z"
np              = 64
group           = 1             # backbone in the default index file
dt              = 100           # dt for trjconv

# os.environ["GMX_DISRE_ENSEMBLE_SIZE"] = str(np)
# os.environ["LD_RUN_PATH"] = "/home/szpari/GROMACS/lib" + program_suffix
# os.environ["LD_LIBRARY_PATH"] = "/home/szpari/GROMACS/lib" + program_suffix

# env_vars =  "setenv GMX_DISRE_ENSEMBLE_SIZE " + str(np) + ";\n" +\
#             "setenv LD_RUN_PATH /home/szpari/GROMACS/lib" + program_suffix + "\n" +\
#             "setenv LD_LIBRARY_PATH /home/szpari/GROMACS/lib" + program_suffix + "\n\n"

# Invoking trjconv
for i in range(np):
    current_file = open("getpdb.sh", 'w')

    trjconv_command = "trjconv" + program_suffix +\
                    " -f " + molname + "_8rep" + str(i) + ".trr" +\
                    " -s " + molname +"_full"  + str(i) + ".tpr" +\
                    " -o " + molname + str(i) + "_dt" + str(dt) + ".pdb" +\
                    " -dt " + str(dt) + " << HERE\n"

    current_file.write("#!/bin/csh\n\n")
    # current_file.write(env_vars)
    current_file.write(trjconv_command)
    current_file.write(str(group))
    current_file.write("\nHERE\n")
    current_file.close()

    os.system("chmod +x getpdb.sh")
    os.system("./getpdb.sh")

os.remove("getpdb.sh")
