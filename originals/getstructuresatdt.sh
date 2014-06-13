#!/bin/csh

set program_suffix = _455mumo1
set molname=1d3z
set np=64;

setenv GMX_DISRE_ENSEMBLE_SIZE $np;
setenv LD_RUN_PATH /home/szpari/GROMACS/lib$program_suffix
setenv LD_LIBRARY_PATH /home/szpari/GROMACS/lib$program_suffix


# backbone in the default index file
set group=1

# dt for trjconv
set dt=100


# Invoking trjconv

set i=0
while ($i < $np)
echo "trjconv$program_suffix -f $molname\_8rep$i\.trr -s $molname\_full$i\.tpr -o $molname$i\_dt$dt\.pdb -dt $dt << HERE " > getpdb.sh
echo "$group\n" >> getpdb.sh
echo "HERE\n" >> getpdb.sh
chmod +x getpdb.sh
./getpdb.sh
@ i = $i + 1
end

rm getpdb.sh
