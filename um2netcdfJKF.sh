#!/bin/bash

# script to convert UM output to netcdf.  written specifically and kludgily for 
# the 2 year short runs December 2006 - November 2008. -Jennifer K Fletcher

### needs to be specified each time:
RUNNAME=uamwa
FILETYPE=daily # files in daily vs monthly chunks

### shouldn't have to be modified
OUTDIR=/short/k10/jkf565/${RUNNAME}
mkdir ${OUTDIR}/ncfiles

module use ~access/modules
module load pythonlib/ScientificPython
module load pythonlib/cdat-lite

# set up to convert name according to UM model output naming convention
if [ "$FILETYPE" = "daily" ]
 then 
 filenames=$(ls ${OUTDIR}/${RUNNAME}a.p*)
 for f in $filenames; do
   foo=$(basename "$f")
   dattype1="$(echo ${foo[@]:8:1})"
   case $dattype1 in 
     a)dattype="met" ;;
     c)dattype="instmet" ;;
     d)dattype="1Dvar" ;;
     e)dattype="cld" ;;
     i)dattype="ISCCPhist" ;;
     *)dattype=$dattype1 ;;
   esac
   dec="200"
   yr="$(echo ${foo[@]:10:1})"
   mon1="$(echo ${foo[@]:11:1})"
   case $mon1 in
     1)mon="01" ;;
     2)mon="02" ;;
     3)mon="03" ;;
     4)mon="04" ;;
     5)mon="05" ;;
     6)mon="06" ;;
     7)mon="07" ;;
     8)mon="08" ;;
     9)mon="09" ;;
     a)mon="10" ;;
     b)mon="11" ;;
     c)mon="12" ;;
     *)mon=$mon1 ;;
   esac
   dy1="$(echo ${foo[@]:12:1})"
   case $dy1 in
     1)dy="01" ;;
     2)dy="02" ;;
     3)dy="03" ;;
     4)dy="04" ;;
     5)dy="05" ;;
     6)dy="06" ;;
     7)dy="07" ;;
     8)dy="08" ;;
     9)dy="09" ;;
     a)dy="10" ;;
     b)dy="11" ;;
     c)dy="12" ;;
     d)dy="13" ;;
     e)dy="14" ;;
     f)dy="15" ;;
     g)dy="16" ;;
     h)dy="17" ;;
     i)dy="18" ;;
     j)dy="19" ;;
     k)dy="20" ;;
     l)dy="21" ;;
     m)dy="22" ;;
     n)dy="23" ;;
     o)dy="24" ;;
     p)dy="25" ;;
     q)dy="26" ;;
     r)dy="27" ;;
     s)dy="28" ;;
     t)dy="29" ;;
     u)dy="30" ;;
     v)dy="31" ;;
     *)dy=$dy1 ;;
  esac
  if [ -e "${OUTDIR}/ncfiles/${RUNNAME}.${dattype}${dec}${yr}${mon}${dy}.nc" ]
  then
    continue
  fi
  if [ ${dattype} = "ISCCPhist" ]
  then
    /short/k10/jkf565/extractfield $f 2337 ${OUTDIR}/ncfiles/${RUNNAME}${dattype}${dec}${yr}${mon}${dy}.nc
  else
    um2netcdf.py -i $f -o  ${OUTDIR}/ncfiles/${RUNNAME}.${dattype}${dec}${yr}${mon}${dy}.nc --nomask
  fi
 done
fi
  

### -- daily average ISCCP histograms

# /short/k10/jkf565/extractfield ${OUTDIR}/${RUNNAME}a.pik6dec 2337 ${OUTDIR}/ncfiles/${RUNNAME}.200612.ISCCPhist.nc
