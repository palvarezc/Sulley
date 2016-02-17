# Fitter arguments
# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
# Reload control samples from tuples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0
Y_SIG=32
Y_PART_RECO=16
Y_COMB=12
Y_JPSILEAK=0
TRIGGER=2
N_TOYS=400
BDT_VAR_NAME="BDTNewu4bR"
BDT_CUT=0.932  #-0.0187   0.297
OUTPUTDIR="toy_result_12_02/toy_result_TIS_BDTNewu4bRB/"
CONSTRAINED=0
NDIMS=2
WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4880
MAX_B_MASS=6200
N_BATCH_JOBS=25

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
   echo -e $CURRENTDIR/bin/toystudy $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $TRIGGER $N_TOYS $BDT_VAR_NAME  $BDT_CUT $CONSTRAINED $NDIMS $FRACDIR"_"$i $WANT_HOP_CUT $MIN_B_MASS $MAX_B_MASS >> $FILESUBMITLINES
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

OUTPUTTREE="toystudyHistBremCatTsallisBkg_results"$NDIMS"D"$TRIGSTR".root"

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

if [ $NDIMS = 1 ]
then
   QUEUE="hepshort.q"
fi

if [ $NDIMS = 2 ]
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

