#example of submit jobs 44 and 45:
#qsub -cwd -q hepshort.q -t 44-45:1 runOnBatch.sh


. ~lhcb/grouplogin/lhcb_login.sh
# ROOTSYS=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_79/ROOT/6.04.02/x86_64-slc6-gcc48-opt
. SetupProject.sh Urania v5r0
ROOTSYS=/vols/build/lhcb/ROOT_6.06.02_patch_histfactory_old
source $ROOTSYS/bin/thisroot.sh



LINE="$(sed $SGE_TASK_ID'q;d' submitLines.dat)"
OUTPUTDIR="$(echo $LINE | cut -d " " -f 12)"

echo $OUTPUTDIR

if [ ! -d $OUTPUTDIR ]
then
   mkdir $OUTPUTDIR
fi

$LINE

