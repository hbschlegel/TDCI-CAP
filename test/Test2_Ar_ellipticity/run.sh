#!/bin/bash


gausspath=$1
tdcipath=$2

source ${gausspath}/bsd/gdv.profile

rm -rf rundir
mkdir rundir
mkdir rundir/g rundir/t_eps1 rundir/t_eps05

cp g/input.gdv rundir/g/input.gdv
cp t_eps1/input rundir/t_eps1/input
cp t_eps05/input rundir/t_eps05/input

cd rundir/g
${gausspath}/gdv < "input.gdv" >"output.log" 2>&1

cd ../t_eps1
ln -s ../g/MatrixElements.faf .
$tdcipath

cd ../t_eps05
ln -s ../g/MatrixElements.faf .
$tdcipath




