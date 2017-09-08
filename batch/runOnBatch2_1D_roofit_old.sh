# Fitter arguments
# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
# Reload control samples from tuples: 0=YES, 1=NO
#"mode: 1 (1D fit Mvis), 2 (2D fit Mvis x misPT), 3 (2D fit simultaneous with Kemu), 4 (binned), 5 (2D fit Mvis x misPT, RooPolyTimesX fit function),
# 6: (2D fit Mvis x misPT, RooPolyTimesX fit function, simultaneous with kemu)"
# 7: (2D fit with HistFactory Mvis x misPT hist)"
# 8: (1D fit with HistFactory Mvis)"
RELOAD_CTRL_SAMPLES=0
Y_SIG=170
Y_PART_RECO=54
Y_COMB=54
Y_JPSILEAK=0
# Y_SIG=199
# Y_PART_RECO=92
# Y_COMB=440
# Y_JPSILEAK=11
TRIGGER=0
N_TOYS=50
BDT_VAR_NAME="BDTNewu4bR"
BDT_CUT=0.7647 #0.1428  #-0.0187   0.297
OUTPUTDIR="/vols/lhcb/palvare1/RK_analysis/Fit_Toys/OldRK"
CONSTRAINED=1
MODE=1
WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4880
MAX_B_MASS=6200
N_BATCH_JOBS=20
N_KEMU=420   #420 for ETOS, 542 for HTOS, 299 for TIS                     #OLD:: 563 for ETOS, 553 for HTOS, 425 for TIS
PPERP_CUT=600

if [ ! -d $OUTPUTDIR ]
then
   mkdir -p $OUTPUTDIR
fi


#create subdirectories for the job

FRACDIR="output"

for i in `seq 1 $N_BATCH_JOBS`;
do
   OUTPUTSUBDIR=$OUTPUTDIR"/"$FRACDIR"_"$i
   if [ ! -d $OUTPUTSUBDIR ]
   then
      mkdir $OUTPUTSUBDIR
   fi
done

#create the file with the lines that run the main

CURRENTDIR=${PWD}"/"
FILESUBMITLINES=$OUTPUTDIR"/submitLines.dat"

if [ -a $FILESUBMITLINES ]
then
   rm $FILESUBMITLINES
fi


for i in `seq 1 $N_BATCH_JOBS`;
do
   echo -e $CURRENTDIR/bin/toystudy $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $TRIGGER $N_TOYS $BDT_VAR_NAME  $BDT_CUT $CONSTRAINED $MODE $FRACDIR"_"$i $WANT_HOP_CUT $MIN_B_MASS $MAX_B_MASS $N_KEMU $PPERP_CUT >> $FILESUBMITLINES
done

#copy the script that reads one line

runOnBatchScript=$OUTPUTDIR"/runOnBatch.sh"
cp batch/runOnBatch.sh $runOnBatchScript 

#prepare hadd  


if [ $TRIGGER = 0 ]
then
   TRIGSTR="L0ETOSOnly_d"
fi

if [ $TRIGGER = 1 ]
then
   TRIGSTR="L0HTOSOnly_d"
fi

if [ $TRIGGER = 2 ]
then
   TRIGSTR="L0TISOnly_d"
fi

if [ $TRIGGER = -1 ]
then
   TRIGSTR="B_plus_ETA"
fi


NDIMS=1
if [ $MODE = 2 ] || [ $MODE = 3 ]
then
   NDIMS=2
fi
OUTPUTTREE="toystudyHistBremCatTsallisBkg_resultsMode"$MODE""$TRIGSTR".root"

HADDSTRING="hadd -k -f "$OUTPUTDIR$OUTPUTTREE

for i in `seq 1 $N_BATCH_JOBS`;
do
   HADDSTRING=$HADDSTRING" "$OUTPUTDIR$FRACDIR"_"$i"/"$OUTPUTTREE
done

#prepare cleanup script

CLEANUPSCRIPT=$OUTPUTDIR"/cleanup.sh"
OUTPUTTREEINSIDE="params_"$TRIGSTR
OUTPUTTABLEDAT=$OUTPUTDIR"/toyresults.dat"
OUTPUTTABLETEX=$OUTPUTDIR"/toyresults.tex"

if [ -a $CLEANUPSCRIPT ]
then
   rm $CLEANUPSCRIPT
fi

echo -e $HADDSTRING >> $CLEANUPSCRIPT
echo -e "mkdir iodir" >> $CLEANUPSCRIPT
echo -e "mv runOnBatch.sh.* iodir/." >> $CLEANUPSCRIPT
echo -e $CURRENTDIR/setup.sh
echo -e $CURRENTDIR/bin/maintables $OUTPUTDIR$OUTPUTTREE $OUTPUTTREEINSIDE $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $OUTPUTTABLEDAT 0 >> $CLEANUPSCRIPT
echo -e $CURRENTDIR/bin/maintables $OUTPUTDIR$OUTPUTTREE $OUTPUTTREEINSIDE $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $OUTPUTTABLETEX 1 >> $CLEANUPSCRIPT
echo -e pdflatex $OUTPUTTABLETEX >> $CLEANUPSCRIPT


#send to the batch

if [ $MODE = 1 ]
then
   # QUEUE="hepshort.q"
   QUEUE="hep.q"
fi

if [ $MODE = 2 ]
then
   QUEUE="hepmedium.q"
fi

if [ $MODE = 3 ]
then
   QUEUE="hepmedium.q"
fi

if [ $MODE = 4 ]
then
   # QUEUE="hepshort.q"
   QUEUE="hep.q"
fi

if [ $MODE = 5 ]
then
   # QUEUE="hepmedium.q"
   QUEUE="hep.q"
fi

if [ $MODE = 6 ]
then
   QUEUE="hepmedium.q"
fi

if [ $MODE = 7 ]
then
   QUEUE="hep.q"
fi

if [ $MODE = 8 ]
then
   QUEUE="hep.q"
fi


QSUBSTRING="qsub -wd "$OUTPUTDIR" -l h_rt=48:00:00 -M paula.alvarez@cern.ch -q "$QUEUE" -t 1-"$N_BATCH_JOBS":1 "$runOnBatchScript

echo  $QSUBSTRING
read -p "Do you want to submit to the batch? " -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
   echo
   $QSUBSTRING
fi

