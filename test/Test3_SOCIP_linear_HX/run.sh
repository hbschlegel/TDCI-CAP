#!/bin/bash

gdvpath=$1
tdcipath=$2

source ${gdvpath}/bsd/gdv.profile

rm -rf rundir
mkdir rundir
mkdir rundir/g_HCl rundir/t_HCl_soc0 rundir/t_HCl_soc1 rundir/g_HF rundir/t_HF_soc0 rundir/t_HF_soc1

cp g_HCl/input.gdv rundir/g_HCl/input.gdv
cp g_HF/input.gdv rundir/g_HF/input.gdv
cp t_HCl_soc0/input rundir/t_HCl_soc0/input
cp t_HCl_soc1/input rundir/t_HCl_soc1/input
cp t_HF_soc0/input rundir/t_HF_soc0/input
cp t_HF_soc1/input rundir/t_HF_soc1/input

cd rundir/g_HCl
${gdvpath}/gdv < "input.gdv" >"output.log" 2>&1

cd ../t_HCl_soc0
ln -s ../g_HCl/MatrixElements.faf .
$tdcipath

cd ../t_HCl_soc1
ln -s ../g_HCl/MatrixElements.faf .
$tdcipath

cd ../g_HF
${gdvpath}/gdv < "input.gdv" >"output.log" 2>&1

cd ../t_HF_soc0
ln -s ../g_HF/MatrixElements.faf .
$tdcipath

cd ../t_HF_soc1
ln -s ../g_HF/MatrixElements.faf .
$tdcipath


