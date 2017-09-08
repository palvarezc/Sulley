#!/bin/sh                                                                                                                            
VERSION=V9_fixPR
  
# Reload control samples from tuples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0

# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
TRIGGER=-1

#Yields(trig0, trig1, trig2)
Y_SIG_1D=(188 49 54)
Y_PART_RECO_1D=(86 5 27)
Y_COMB_1D=(180 130 43)
Y_JPSILEAK_1D=(8 5 2)
# For 2D working point
Y_SIG_2D=(199 55 51)
Y_PART_RECO_2D=(92 6 26)
Y_COMB_2D=(440 300 101)
Y_JPSILEAK_2D=(16 6 3)

N_TOYS=50
BDT_VAR_NAME="BDTNewu4bR"
# BDT_CUTS=(0.7647 0.7647) #0.1428  #-0.0187   0.297
BDT_CUTS_2D=(0.142857 0.0756303 0.579832)
BDT_CUTS_1D=(0.764706 0.865546 0.932773)

OUTPUT=/vols/lhcb/palvare1/RK_analysis/Fit_Toys/$VERSION
CONSTRAINED=1
QUEUE="hep.q"

#"mode: 
# 1 (1D fit Mvis), 2 (2D fit Mvis x misPT), 3 (2D fit simultaneous with Kemu)
# 4 (binned), 5 (2D fit Mvis x misPT, RooPolyTimesX fit function),
# 6: (2D fit Mvis x misPT, RooPolyTimesX fit function, simultaneous with kemu)"
# 7: (2D fit with HistFactory Mvis x misPT hist)"
# 8: (1D fit with HistFactory Mvis)"
MODE_CODES=(1 5)
MODE_NAMES=(1D_RooFit 2D_RooFit) 

WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4880
MAX_B_MASS=6200
N_BATCH_JOBS=20
N_KEMU=(420 542 299)   #420 for ETOS, 542 for HTOS, 299 for TIS                     #OLD:: 563 for ETOS, 553 for HTOS, 425 for TIS
PPERP_CUT=600

if [ -d $OUTPUT ]
then
echo "Output folder:"$OUTPUT" already exist, please delete it or update OUTPUT"
return 1
fi

mkdir -p $OUTPUT

cp batch/processall.sh $OUTPUT/processall.sh
sed -i "s/VERSION/$VERSION/g" $OUTPUT/processall.sh


read -p "Do you want to submit to the batch? " -n 1 -r
if [[ !($REPLY =~ ^[Yy]$) ]]
then
    echo
    return 1
fi
echo


    
# for RUN in $(echo 1 2);do                                                                                                              
for RUN in $(echo 1);do                                                                                                              
    for TRIGGER in $(echo 0 1 2);do 
        # 1D modes
	for MODEINDX in $(echo 0 1);do

	    MODE=${MODE_CODES[MODEINDX]} 
	    OUTPUTDIR=$OUTPUT/Run$RUN/Trig$TRIGGER/${MODE_NAMES[MODEINDX]}/
	    let "N_KEMU = N_KEMU[TRIGGER]*RUN"

	    mkdir -p $OUTPUTDIR
	    
	    if [ $MODEINDX = 0 ] || [ $MODEINDX = 2 ]
		then
		let "Y_SIG = Y_SIG_1D[TRIGGER]*RUN"
		let "Y_PART_RECO = Y_PART_RECO_1D[TRIGGER]*RUN"
		let "Y_COMB = Y_COMB_1D[TRIGGER]*RUN"
		let "Y_JPSILEAK = Y_JPSILEAK_1D[TRIGGER]*RUN"
		
		BDT_CUT=${BDT_CUTS_1D[$TRIGGER]}	    
		fi
	    
	    if [ $MODEINDX = 1 ] || [ $MODEINDX = 3 ]
		then
		let "Y_SIG = Y_SIG_2D[TRIGGER]*RUN"
		let "Y_PART_RECO = Y_PART_RECO_2D[TRIGGER]*RUN"
		let "Y_COMB = Y_COMB_2D[TRIGGER]*RUN"
		let "Y_JPSILEAK = Y_JPSILEAK_2D[TRIGGER]*RUN"

		BDT_CUT=${BDT_CUTS_2D[$TRIGGER]}	    
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
		TRIGSTR="passTrigCat0"
	    fi
	    
	    if [ $TRIGGER = 1 ]
	    then
		TRIGSTR="passTrigCat1"
	    fi

	    if [ $TRIGGER = 2 ]
	    then
		TRIGSTR="passTrigCat2"
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

	    QSUBSTRING="qsub -wd "$OUTPUTDIR" -l h_rt=48:00:00 -M paula.alvarez@cern.ch -q "$QUEUE" -t 1-"$N_BATCH_JOBS":1 "$runOnBatchScript

	    echo "Submitting..."
	    echo $QSUBSTRING
	    $QSUBSTRING
	    
 	done	
    done
done



