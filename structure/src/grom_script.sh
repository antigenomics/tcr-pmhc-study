#!/bin/bash

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
    	"$@" 2>$LOG_PATH
        echo "error with: $@"
        echo
        mv em$NAME.edr groups$NAME.dat trash
        exit 0
    fi
    return 0
}

function pytest {
    SEQS=$("$@")
    local status=$?
    if [ $status -ne 0 ]; then
    	echo
    	"$@" 2>$GROMACS_PATH/$LOG_PATH
        echo "error with: $@"
        echo
        mv groups$NAME.dat trash
        exit 0
    fi
    return 0
}

#split='====================================================n'
NAME=$3
LOG_PATH=fails/$NAME.log
XPM_PATH=$2
GROMACS_PATH=$1

cd $GROMACS_PATH 

mkdir fails
test gmx pdb2gmx -f $NAME.pdb -o $NAME.gro -water spce -missing -ff oplsaa
test gmx editconf -f $NAME.gro -o $NAME.gro -c -d 1.0 -bt cubic

cd -

pytest python enemat.py $*

cd -

test gmx grompp -f params/minim.mdp -c $NAME.gro -p topol.top -n index.ndx -o em$NAME.tpr
test gmx mdrun -v -deffnm em$NAME -nb cpu
test gmx enemat -groups groups$NAME.dat -nlevels 10000 -emat $NAME.xpm -f em$NAME.edr
echo $SEQS >> total*.xpm
mv em$NAME.edr groups$NAME.dat trash

cd -

mkdir $XPM_PATH
mv ${GROMACS_PATH}/*.xpm $XPM_PATH
echo 'XPM files can be found in '${XPM_PATH}

rm -f ${GROMACS_PATH}/*

