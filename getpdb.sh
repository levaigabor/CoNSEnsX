#!/bin/csh

setenv GMX_DISRE_ENSEMBLE_SIZE 64;
setenv LD_RUN_PATH /home/szpari/GROMACS/lib_455mumo1
setenv LD_LIBRARY_PATH /home/szpari/GROMACS/lib_455mumo1

trjconv_455mumo1 -f 1d3z_8rep63.trr -s 1d3z_full63.tpr -o 1d3z63_dt100.pdb -dt 100 << HERE
1
HERE
