# Fitter arguments
# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
# Reload control samples from tuples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0
Y_SIG=19900
Y_PART_RECO=9200
Y_COMB=46000
Y_JPSILEAK=1600
TRIGGER=0
BDT_VAR_NAME="BDTNewu4bR"
BDT_CUT=0.1428  #-0.0187   0.297
OUTPUTDIR="toy_result_12_02/toy_result_ETOS_BDTNewu4bR_Syst_TMinusOneSigma/"
CONSTRAINED=0
WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4880
MAX_B_MASS=6200
PARAM_TO_CHANGE="T"
N_SIGMA=-1
N_TOYS=100
N_BATCH_JOBS=25








MODE=2

if [ ! -d $OUTPUTDIR ]
then
   mkdir $OUTPUTDIR
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
   echo -e $CURRENTDIR/bin/systematicstudy $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $TRIGGER $N_TOYS $BDT_VAR_NAME  $BDT_CUT $CONSTRAINED $PARAM_TO_CHANGE $N_SIGMA $FRACDIR"_"$i $WANT_HOP_CUT $MIN_B_MASS $MAX_B_MASS >> $FILESUBMITLINES
done

#copy the script that reads one line

runOnBatchScript=$OUTPUTDIR"/runOnBatch.sh"
cp runOnBatch.sh $runOnBatchScript 

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

OUTPUTTREE="toystudyHistBremCatTsallisBkg_results"$MODE"D"$TRIGSTR".root"

HADDSTRING="hadd -k -f "$CURRENTDIR$OUTPUTDIR$OUTPUTTREE

for i in `seq 1 $N_BATCH_JOBS`;
do
   HADDSTRING=$HADDSTRING" "$CURRENTDIR$OUTPUTDIR$FRACDIR"_"$i"/"$OUTPUTTREE
done

#prepare cleanup script

CLEANUPSCRIPT=$OUTPUTDIR"/cleanup.sh"
OUTPUTTREEINSIDE="params_"$TRIGSTR
OUTPUTTABLEDAT=$CURRENTDIR$OUTPUTDIR"/toyresults.dat"
OUTPUTTABLETEX=$CURRENTDIR$OUTPUTDIR"/toyresults.tex"

if [ -a $CLEANUPSCRIPT ]
then
   rm $CLEANUPSCRIPT
fi

echo -e $HADDSTRING >> $CLEANUPSCRIPT
echo -e "mkdir iodir" >> $CLEANUPSCRIPT
echo -e "mv runOnBatch.sh.* iodir/." >> $CLEANUPSCRIPT
echo -e $CURRENTDIR/setup.sh
echo -e $CURRENTDIR/bin/maintables $CURRENTDIR$OUTPUTDIR$OUTPUTTREE $OUTPUTTREEINSIDE $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $OUTPUTTABLEDAT 0 >> $CLEANUPSCRIPT
echo -e $CURRENTDIR/bin/maintables $CURRENTDIR$OUTPUTDIR$OUTPUTTREE $OUTPUTTREEINSIDE $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $OUTPUTTABLETEX 1 >> $CLEANUPSCRIPT
echo -e pdflatex $OUTPUTTABLETEX >> $CLEANUPSCRIPT


#send to the batch

if [ $MODE = 1 ]
then
   QUEUE="hepshort.q"
fi

if [ $MODE = 2 ]
then
   QUEUE="hepmedium.q"
fi

if [ $MODE = 3 ]
then
   QUEUE="hepmedium.q"
fi

QSUBSTRING="qsub -wd "$CURRENTDIR"/"$OUTPUTDIR" -q "$QUEUE" -t 1-"$N_BATCH_JOBS":1 "$runOnBatchScript

echo  $QSUBSTRING
read -p "Do you want to submit to the batch? " -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
   echo
   $QSUBSTRING
fi

