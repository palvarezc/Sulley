#!/bin/bash

# Fitter arguments
# Trigger category: 0=ETOS, 1=HTOS, 2=TIS
# Reload control samples from tuples: 0=YES, 1=NO
RELOAD_CTRL_SAMPLES=0
Y_SIG=199
Y_PART_RECO=92
Y_COMB=440
Y_JPSILEAK=16
TRIGGER=0
N_TOYS=2
BDT_VAR_NAME="BDTNewu4bR"
BDT_CUT=0.1428  #-0.0187   0.297
OUTPUTDIR="toy_result/"
CONSTRAINED=0
MODE=7
WANT_HOP_CUT=0 #0: no HOP cut, 1: HOP cut
MIN_B_MASS=4880
MAX_B_MASS=6200
N_KEMU=420   #420 for ETOS, 542 for HTOS, 299 for TIS                     #OLD:: 563 for ETOS, 553 for HTOS, 425 for TIS
PPERP_CUT=600

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
echo "       BDT Variable :"$BDT_VAR_NAME
echo "       BDT cut : "${BDT_CUT}
echo "       Reloading control samples? "$reply 
echo "       Constraining Part. Reco? "$CONSTRAINED
echo "       # dimensions: "$MODE
echo "       HOP cut applied? "$replyHOP
echo "       Visible mass range : "$MIN_B_MASS" < M < "$MAX_B_MASS
echo "================================================================================"

$SULLEYROOT/bin/toystudy $RELOAD_CTRL_SAMPLES $Y_SIG $Y_PART_RECO $Y_COMB $Y_JPSILEAK $TRIGGER $N_TOYS $BDT_VAR_NAME  $BDT_CUT $CONSTRAINED $MODE $OUTPUTDIR $WANT_HOP_CUT $MIN_B_MASS $MAX_B_MASS $N_KEMU $PPERP_CUT

