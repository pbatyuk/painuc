#!/bin/bash
#
# set -x   # comment out to display script lines
#
#=========================================================================
# function declaration
#
COMPILE_ONE()
{
  # 1st argument $1  is the sign ("positive" or "negative")
  # 2nd argument $2  is the energy ("106" or "68")
  # 3rd argument $3  is the destination directory
  LOGFILE=`pwd`/Makefile_allinclude_$1_$2.log
  EXEFILE=digitization_$1_allinclude
  DESTFILE=measure_${1:0:3}_$2
  echo "Making executable for $1 ; see log in ${LOGFILE}"
  make -f Makefile_allinclude_$1 clean
  make -f Makefile_allinclude_$1 >& ${LOGFILE}
  if test "X`grep warning\: ${LOGFILE}`" != "X"; then
    echo "... Warnings from make : see ${LOGFILE} for more details"
  fi
  if test "X`grep make ${LOGFILE} | grep Error`" != "X"; then
    echo "... Errors from make : ../${EXEFILE} executable was not created"
    echo "      see ${LOGFILE} for more details"
  else
    echo "... Creating executable `pwd`/${DESTFILE}"
    mv ../${EXEFILE} `pwd`/${DESTFILE}
 #   echo "... Saving previous executable in $3/${DESTFILE}"
 #   TODAY=`ls -l $3/${DESTFILE} | mawk '{ print( $6 "-" $7) }' | tr : - | tr - _ `
    if [ -e $3/${DESTFILE} ] ; then
      mv $3/${DESTFILE} $3/${DESTFILE}_${TODAY}
    fi
    echo "... Copying new executable into $3"
    cp `pwd`/${DESTFILE} $3/${DESTFILE}
  fi
}


TOPDIR=`pwd`/..
cd ${TOPDIR}/src
echo "======================================================================"
echo "Script Build_all.sh for running make of all 4 versions "
echo "of the digitization_*_allinclude_* program "
echo "WARNING: this script is running in `pwd` "
echo "NOTE: you may need to enter the sudo password to save old executable"
echo "      and copy new executable in ${TOPDIR}/ directory"
echo "======================================================================"
#
echo " "
echo "Running make for positive pions @ 106 MeV ..."
export CXXFLAGS=
SEGNO=positive
ENERGY=106
# make -f Makefile_allinclude_${SEGNO} >& $LOGFILE
COMPILE_ONE ${SEGNO} ${ENERGY} ${TOPDIR}
#
#
echo " "
echo "Running make for positive pions @ 68 MeV ..."
export CXXFLAGS=-DPI70
SEGNO=positive
ENERGY=68
COMPILE_ONE ${SEGNO} ${ENERGY} ${TOPDIR}
#
#
echo " "
echo "Running make for negative pions @ 106 MeV ..."
export CXXFLAGS=
SEGNO=negative
ENERGY=106
COMPILE_ONE ${SEGNO} ${ENERGY} ${TOPDIR}
#
#
echo " "
echo "Running make for negative pions @ 68 MeV ..."
export CXXFLAGS=-DPI70
SEGNO=negative
ENERGY=68
COMPILE_ONE ${SEGNO} ${ENERGY} ${TOPDIR}
#
#
echo "================================"
echo "Exiting from script Build_all.sh"
echo "================================"
echo " "
#
exit
