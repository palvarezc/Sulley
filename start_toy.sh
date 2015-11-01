#!/bin/bash

# Fitter arguments
# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
# Reload control samples from tuples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0
Y_SIG=56
Y_PART_RECO=52
Y_COMB=12
Y_JPSILEAK=10
TRIGGER=0
N_TOYS=5
BDT_CUT=0.369
OUTPUT_FOLDER="toy_result/"
CONSTRAINED=0
NDIMS=1
WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4500
MAX_B_MASS=6200

if [[ "$RELOAD_CTRL_SAMPLES" == '1' ]]; 
then 
reply="No"; 
else 
reply="Yes"; 
fi

if [[ "$WANT_HOP_CUT" == '1' ]]; 
then 
replyHOP="Yes"; 
else 
replyHOP="No"; 
fi

echo "================================================================================"
echo "   Starting "${N_TOYS}" toys with the following configuration: "
echo "   "
echo "       # signal events         : "${Y_SIG}
echo "       # part. reco. events    : "${Y_PART_RECO}
echo "       # comb. bkg. events     : "${Y_COMB}
echo "       # J psi Leakage events  : "${Y_JPSILEAK}
echo "       Trigger category (0=ETOS, 1=HTOS, 2=TIS) : "${TRIGGER}" "
echo "       BDT cut : "${BDT_CUT}
echo "       Reloading control samples? "$reply 
echo "       Constraining Part. Reco? "$CONSTRAINED
echo "       # dimensions: "$NDIMS
echo "       HOP cut applied? "$replyHOP
echo "       Visible mass range : "$MIN_B_MASS" < M < "$MAX_B_MASS
echo "================================================================================"

$SULLEYROOT/bin/toystudy $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $TRIGGER $N_TOYS $BDT_CUT $CONSTRAINED $NDIMS $OUTPUT_FOLDER $WANT_HOP_CUT $MIN_B_MASS $MAX_B_MASS
