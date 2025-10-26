#!/bin/bash


gausspath=$1
tdcipath=$2

source ${gausspath}/bsd/gdv.profile

rm -rf rundir
mkdir rundir
mkdir rundir/g rundir/t

cp g/input.gdv rundir/g/input.gdv
cp t/input rundir/t

cd rundir/g
${gausspath}/gdv < "input.gdv" >"output.log" 2>&1

cd ../t
ln -s ../g/MatrixElements.faf .
$tdcipath




