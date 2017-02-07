from ROOT import *
from math import *

gROOT.ProcessLine(".x /home/hep/palvare1/B2Kmumu_Kee/Sulley/src/RooMcorMvisTsallis.cpp++")

minBMass = 4500
maxBMass = 6200

trigStr = "B_plus_ETA"
weightStr = "PIDTrigDataCondWeight_TISOnly"

Nsig = 56
Npartreco = 52
Ncombinatorial = 12

fzero = 0.368
fone = 0.484
ftwo = 1- fzero - fone

superdict = {"signal_zero" : {"filename" : "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/Kee_RedoCalo_signalPresel_HltTrigged.root" ,
                              "cuts"     : " ("+trigStr+"  > 0.9) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5) && B_plus_M > "+str(minBMass)+" && B_plus_M < "+str(maxBMass),
                              "stat" : floor(Nsig*fzero)},

             "signal_one" : {"filename" : "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/Kee_RedoCalo_signalPresel_HltTrigged.root" ,
                             "cuts"     : " ("+trigStr+"  > 0.9) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5) && B_plus_M > "+str(minBMass)+" && B_plus_M < "+str(maxBMass),
                              "stat" : floor(Nsig*fone)},

             "signal_two" : {"filename" : "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/Kee_RedoCalo_signalPresel_HltTrigged.root" ,
                             "cuts"     : " ("+trigStr+"  > 0.9) &&  ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5) && B_plus_M > "+str(minBMass)+" && B_plus_M < "+str(maxBMass),
                              "stat" : floor(Nsig*ftwo)},

             "partreco" : {"filename" : "/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root",
                           "cuts"     : "("+trigStr+"  > 0.9) &&  B_plus_M > "+str(minBMass)+" && B_plus_M < "+str(maxBMass),
                           "stat" : Npartreco},
             }


nbins_mass = int(floor((maxBMass-minBMass)/(40.)))
nbins_corr = int(floor((10000-minBMass)/120.))

B_plus_M = RooRealVar("B_plus_M","B_plus_M",minBMass,maxBMass)
B_plus_M_corr = RooRealVar("B_plus_M_corr","B_plus_M_corr",minBMass,10000)
defaultMBins = RooBinning(nbins_mass, B_plus_M.getMin(), B_plus_M.getMax() ); 
defaultMCorrBins =  RooBinning(nbins_corr, B_plus_M_corr.getMin(), B_plus_M_corr.getMax()); 

B_plus_M.setBinning( defaultMBins);
B_plus_M_corr.setBinning( defaultMCorrBins );

expoConst = RooRealVar("expoConst", "expoConst", -1e-3, -1, 0.1);
T = RooRealVar("T", "T", 97, 0, 200); 
n = RooRealVar("n", "n", 3.5, 1., 5.5);
combPDF = RooMcorMvisTsallis("combPDF", "combPDF", B_plus_M_corr, B_plus_M, T, n, expoConst);

arglist = RooArgList(B_plus_M, B_plus_M_corr)
argset = RooArgSet(B_plus_M, B_plus_M_corr)


fout = TFile("histograms.root","recreate")

# htoy = TH2D("dataToy","dataToy",nbins_mass,minBMass, maxBMass, nbins_corr, minBMass, 10000)
toy_ds = RooDataSet("toy_ds","toy_ds",argset)

for sample in superdict.keys():

    fin = TFile(superdict[sample]["filename"])
    tin = fin.Get("DecayTree") 
    tin.SetBranchStatus("*", 0);
    tin.SetBranchStatus("B_plus_M", 1); tin.SetBranchStatus("B_plus_M_corr", 1);
    tin.SetBranchStatus("B_plus_DTFM_M_zero", 1); 
    tin.SetBranchStatus("e_plus_BremMultiplicity", 1); tin.SetBranchStatus("e_minus_BremMultiplicity", 1); 
    tin.SetBranchStatus(trigStr,1); tin.SetBranchStatus(weightStr,1);    
    h2 = TH2D(sample, sample, nbins_mass,minBMass, maxBMass, nbins_corr, minBMass, 10000)
    tin.Draw("B_plus_M_corr:B_plus_M>>"+sample,"("+superdict[sample]["cuts"]+")*"+weightStr)
    fout.cd("")
    h2.Sumw2()
    for i in range(1,nbins_mass+1):
        for j in range(1,nbins_corr+1):
            value = h2.GetBinContent(i,j)
            if value==0.: h2.SetBinContent(i,j,1e-9)
            print "Adding 1e-9 to a zero bin"

    h2.Write()
    roohist = RooDataHist("roohist","roohist",arglist,h2)
    histpdf = RooHistPdf("histpdf","histpdf",argset,roohist)
    toy_ds_i = histpdf.generate(argset,int(superdict[sample]["stat"]),1,0)
    toy_ds.append(toy_ds_i)


toy_ds_i = combPDF.generate(argset,Ncombinatorial)
toy_ds.append(toy_ds_i)

htoy = toy_ds.createHistogram(B_plus_M,B_plus_M_corr)
htoy.Write("dataToy")

fout.Save()
fout.Close()
a = TBrowser()

