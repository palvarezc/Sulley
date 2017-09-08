#!/bin/sh                                                                                                                         

MODE_NAMES=(1D_RooFit 2D_RooFit)
TOYMCDIR=${PWD}
HTMLDIR=/home/hep/palvare1/public_html/RK_fit/NewBDT/VERSION

for RUN in $(echo 1 2);do                                                                                   
    for TRIGGER in $(echo 0 1 2);do
        for MODEINDX in $(echo 0 1);do
	    cd $TOYMCDIR/Run$RUN/Trig$TRIGGER/${MODE_NAMES[MODEINDX]}
	    source cleanup.sh
	    if [ ! -d $HTMLDIR/Run$RUN/Trig$TRIGGER/${MODE_NAMES[MODEINDX]} ];
		then
		mkdir -p $HTMLDIR/Run$RUN/Trig$TRIGGER/${MODE_NAMES[MODEINDX]}
		fi
	    \cp *pdf $HTMLDIR/Run$RUN/Trig$TRIGGER/${MODE_NAMES[MODEINDX]}
	    cd -
	done
    done 
done
