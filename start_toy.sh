#!/bin/bash

# Fitter arguments
# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
# Reload control samples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0
Y_SIG=56
Y_PART_RECO=52
Y_COMB=12
TRIGGER=2
N_TOYS=1000
BDT_CUT=0.369

if [[ "$RELOAD_CTRL_SAMPLES" == '1' ]]; 
then 
reply="No"; 
else 
reply="Yes"; 
fi

echo "================================================================================"
echo "   Starting "${N_TOYS}" toys with the following configuration: "
echo "   "
echo "       # signal events      : "${Y_SIG}
echo "       # part. reco. events : "${Y_PART_RECO}
echo "       # comb. bkg. events  : "${Y_COMB}
echo "       Trigger category (0=ETOS, 1=HTOS, 2=TIS) : "${TRIGGER}" "
echo "       BDT cut : "${BDT_CUT}
echo "       Reloading control samples? "$reply 
echo "================================================================================"

./bin/toystudy2DHistBremCatTsallisBkg $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $TRIGGER $N_TOYS $BDT_CUT
