#!/bin/sh                                                                                                                            
VERSION=UBDTs
  
# Reload control samples from tuples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0

# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
TRIGGER=-1

PART_RECO_FRACS=(0.46 0.46 0.46)

#Yields(trig0, trig1, trig2)
Y_SIG_2D_RUN1=(221 225 226)
Y_COMB_2D_RUN1=(400 338 325)
Y_JPSILEAK_2D_RUN1=(8 10 9)
BDT_VAR_NAMES_2D_RUN1=("BDTu1alpha0Run1R" "BDTu1alpha4Run1R" "BDTu1alpha8Run1R")
BDT_CUTS_2D_RUN1=(0.23 0.32 0.31)


Y_SIG_2D_RUN2=(347 358 344)
Y_COMB_2D_RUN2=(910 770 387)
Y_JPSILEAK_2D_RUN2=(17 17 17)
BDT_VAR_NAMES_2D_RUN2=("BDTu1alpha4Run2R" "BDTu1alpha32Run2R" "BDTu1alpha64Run2R")
BDT_CUTS_2D_RUN2=(0.15 0.18 0.66)


OUTPUT=/vols/lhcb/palvare1/RK_analysis/Fit_Toys/NewBDT/$VERSION
CONSTRAINED=0
QUEUE="hep.q"

#"mode: 
# 1 (1D fit Mvis), 2 (2D fit Mvis x misPT), 3 (2D fit simultaneous with Kemu)
# 4 (binned), 5 (2D fit Mvis x misPT, RooPolyTimesX fit function),
# 6: (2D fit Mvis x misPT, RooPolyTimesX fit function, simultaneous with kemu)"
# 7: (2D fit with HistFactory Mvis x misPT hist)"
# 8: (1D fit with HistFactory Mvis)"
MODE_CODES=(1 5)
MODE_NAMES=(1D_RooFit 2D_RooFit) 

N_TOYS=50
N_BATCH_JOBS=20
WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4880
MAX_B_MASS=6200
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

MODEINDX=1
TRIGGER=0
for RUN in $(echo 1 2);do
    for UBDT in $(echo 0 1 2);do 
        # 1D modes

	    MODE=${MODE_CODES[MODEINDX]} 
	    OUTPUTDIR=$OUTPUT/Run$RUN/BDT$UBDT/${MODE_NAMES[MODEINDX]}/
	    let "N_KEMU = N_KEMU[UBDT]"

	    mkdir -p $OUTPUTDIR
	    
	    
	    if [ $RUN = 1 ] 
		then
		let "Y_SIG = Y_SIG_2D_RUN1[UBDT]"
		Y_PART_RECO=$(echo ${PART_RECO_FRACS[UBDT]}"*"${Y_SIG_2D_RUN1[UBDT]} | bc)
		let "Y_COMB = Y_COMB_2D_RUN1[UBDT]"
		let "Y_JPSILEAK = Y_JPSILEAK_2D_RUN1[UBDT]"

		BDT_VAR_NAME=${BDT_VAR_NAMES_2D_RUN1[$UBDT]}	    
		BDT_CUT=${BDT_CUTS_2D_RUN1[$UBDT]}	    
		fi

	    if [ $RUN = 2 ] 
		then
		let "Y_SIG = Y_SIG_2D_RUN2[UBDT]"
		Y_PART_RECO=$(echo ${PART_RECO_FRACS[UBDT]}"*"${Y_SIG_2D_RUN2[UBDT]} | bc)
		let "Y_COMB = Y_COMB_2D_RUN2[UBDT]"
		let "Y_JPSILEAK = Y_JPSILEAK_2D_RUN2[UBDT]"

		BDT_VAR_NAME=${BDT_VAR_NAMES_2D_RUN2[$UBDT]}	    
		BDT_CUT=${BDT_CUTS_2D_RUN2[$UBDT]}	    
		fi

		
	    Y_PART_RECO=`printf '%.0f\n' $Y_PART_RECO` #makes the yield integer
	    # echo "Y part reco is "$Y_PART_RECO
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
		echo -e $CURRENTDIR/bin/toystudy $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $TRIGGER $N_TOYS $BDT_VAR_NAME  $BDT_CUT $CONSTRAINED $MODE $FRACDIR"_"$i $WANT_HOP_CUT $MIN_B_MASS $MAX_B_MASS  $RUN $N_KEMU $PPERP_CUT >> $FILESUBMITLINES
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



