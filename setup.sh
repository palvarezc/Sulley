#!/bin/bash
SCRIPT=Sulley
SCRIPTPATH=${PWD}

echo ""
echo "==============================================================================="
echo "    Setting up paths and environment for "$SCRIPT"."
echo "==============================================================================="
echo ""
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCRIPTPATH/lib
if [ ! -d $SCRIPTPATH/lib ]
then
  mkdir $SCRIPTPATH/lib
fi
if [ ! -d $SCRIPTPATH/obj ]
then
  mkdir $SCRIPTPATH/obj
fi
if [ ! -d $SCRIPTPATH/bin ]
then
  mkdir $SCRIPTPATH/bin
fi

export SULLEYROOT=${PWD}


# export ROOTSYS=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_79/ROOT/6.04.02/x86_64-slc6-gcc48-opt/root
# source $ROOTSYS/bin/thisroot.sh


