#include "fitter_utils_1Dsimultaneous.h"
#include "RooSimultaneous.h"


FitterUtils1DSimultaneous::FitterUtils1DSimultaneous(int nGenKemu_, int nGenSignal_, int nGenPartReco_, int nGenComb_, int nGenJpsiLeak_, double nGenFracZeroGamma_, double nGenFracOneGamma_, double PPerpCut_, string workspacename_)
   :FitterUtils(nGenSignal_, nGenPartReco_, nGenComb_, nGenJpsiLeak_, nGenFracZeroGamma_, nGenFracOneGamma_, false, workspacename_),
   nGenFracZeroGammaLowPPerp(0),
   nGenFracOneGammaLowPPerp(0),
   nGenFracZeroGammaHighPPerp(0),
   nGenFracOneGammaHighPPerp(0),
   nGenSignalLowPPerp(0),
   nGenSignalZeroGammaLowPPerp(0),
   nGenSignalOneGammaLowPPerp(0),
   nGenSignalTwoGammaLowPPerp(0),
   nGenPartRecoLowPPerp(0),
   nGenCombLowPPerp(0),
   nGenJpsiLeakLowPPerp(0),
   nGenSignalHighPPerp(0),
   nGenSignalZeroGammaHighPPerp(0),
   nGenSignalOneGammaHighPPerp(0),
   nGenSignalTwoGammaHighPPerp(0),
   nGenPartRecoHighPPerp(0),
   nGenCombHighPPerp(0),
   nGenJpsiLeakHighPPerp(0),
   nGenFracSigLowOverTot(0),
   nGenFracSigZeroGammaLowOverTot(0),
   nGenFracSigOneGammaLowOverTot(0),
   nGenFracSigTwoGammaLowOverTot(0),
   nGenFracCombLowOverTot(0),
   nGenFracPartRecoLowOverTot(0),
   nGenFracJpsiLeakLowOverTot(0),
   nGenSignalZeroGamma(0),
   nGenSignalOneGamma(0),
   nGenSignalTwoGamma(0),
   FracCombFromFitLowOverTot(0),
   nCombFromFitLowPPerp(0),
   nCombFromFitHighPerp(0),
   nGenKemu(nGenKemu_),
   PPerpCut(PPerpCut_),
   FracCombRealVarLowOverTot("FracCombFormulaLowOverTot", "FracCombFormulaLowOverTot", -1000, 1000)
{}


void FitterUtils1DSimultaneous::prepare_PDFs(string trigStr, string BDTVar, double BDTcut,
      string signalfile, string partrecofile, string combfile, string JpsiLeakfile,
      double minBMass, double maxBMass,
      string signaltree, string partrecotree, string combtree, string JpsiLeaktree)
{


   //***********Get the datasets
   TFile* fSignal = new TFile(signalfile.c_str());
   TTree* tSignal = (TTree*)fSignal->Get(signaltree.c_str());
   TFile* fPartReco = new TFile(partrecofile.c_str());
   TTree* tPartReco = (TTree*)fPartReco->Get(partrecotree.c_str());
   TFile* fComb = new TFile(combfile.c_str());
   TTree* tComb = (TTree*)fComb->Get(combtree.c_str()); 
   TFile* fJpsiLeak = new TFile(JpsiLeakfile.c_str());
   TTree* tJpsiLeak = (TTree*)fJpsiLeak->Get(JpsiLeaktree.c_str()); 


   //**********Define variables
   RooRealVar trigVar(trigStr.c_str(), trigStr.c_str(), -10, 10);
   RooRealVar BDTRooRealVar(BDTVar.c_str(), BDTVar.c_str(), -1,1);
   RooRealVar B_plus_M("B_plus_M", "M_{visible}", minBMass, maxBMass, "MeV/c^{2}");
   RooRealVar B_plus_M_corr("B_plus_M_corr", "M_{cor}", minBMass, 10000, "MeV/c^{2}");
   RooRealVar B_plus_DTFM_M_zero("B_plus_DTFM_M_zero", "M_{constr}", 0, 20000, "MeV/c^{2}"); 
   RooRealVar e_plus_BremMultiplicity("e_plus_BremMultiplicity","e_plus_BremMultiplicity", -1,2);
   RooRealVar e_minus_BremMultiplicity("e_minus_BremMultiplicity","e_minus_BremMultiplicity", -1,2);
   RooRealVar misPT("misPT", "misPT", 0., 6000.);

   misPT.setRange("lowPPerp", 0., PPerpCut);
   misPT.setRange("highPPerp", PPerpCut, 6000);
   misPT.setRange("default", 0, 6000);

   RooRealVar weightPartReco("weightPartReco", "weightPartReco", 0, 10);
   RooRealVar weightLeakage("weightLeakage", "weightLeakage", 0, 10);
   RooRealVar dataMCWeightee("DataMCWeightee", "DataMCWeightee",0, 30);


   //***********Set only variables needed

   tSignal->SetBranchStatus("*", 0); tSignal->SetBranchStatus("B_plus_M", 1); tSignal->SetBranchStatus("B_plus_M_corr", 1); tSignal->SetBranchStatus("B_plus_DTFM_M_zero", 1); tSignal->SetBranchStatus(BDTVar.c_str(),1);
   tSignal->SetBranchStatus("e_plus_BremMultiplicity", 1); tSignal->SetBranchStatus("e_minus_BremMultiplicity", 1); tSignal->SetBranchStatus(trigStr.c_str()); tSignal->SetBranchStatus("DataMCWeightee",1);
   tSignal->SetBranchStatus("misPT",1);

   tPartReco->SetBranchStatus("*", 0); tPartReco->SetBranchStatus("B_plus_M", 1); tPartReco->SetBranchStatus("B_plus_M_corr", 1); tPartReco->SetBranchStatus("B_plus_DTFM_M_zero", 1);tPartReco->SetBranchStatus(BDTVar.c_str(),1);
   tPartReco->SetBranchStatus("e_plus_BremMultiplicity", 1); tPartReco->SetBranchStatus("e_minus_BremMultiplicity", 1); tPartReco->SetBranchStatus(trigStr.c_str()); tPartReco->SetBranchStatus("weightPartReco",1);
   tPartReco->SetBranchStatus("misPT",1);

   tComb->SetBranchStatus("*", 0); tComb->SetBranchStatus("B_plus_M", 1); tComb->SetBranchStatus("B_plus_M_corr", 1); tComb->SetBranchStatus("B_plus_DTFM_M_zero", 1);tComb->SetBranchStatus(BDTVar.c_str(),1);
   tComb->SetBranchStatus("e_plus_BremMultiplicity", 1); tComb->SetBranchStatus("e_minus_BremMultiplicity", 1); tComb->SetBranchStatus(trigStr.c_str());
   tComb->SetBranchStatus("misPT",1);

   tJpsiLeak->SetBranchStatus("*", 0); tJpsiLeak->SetBranchStatus("B_plus_M", 1); tJpsiLeak->SetBranchStatus("B_plus_M_corr", 1); 
   tJpsiLeak->SetBranchStatus("B_plus_DTFM_M_zero", 1);tJpsiLeak->SetBranchStatus(BDTVar.c_str(),1);
   tJpsiLeak->SetBranchStatus("e_plus_BremMultiplicity", 1); tJpsiLeak->SetBranchStatus("e_minus_BremMultiplicity", 1); tJpsiLeak->SetBranchStatus(trigStr.c_str()); 
   tJpsiLeak->SetBranchStatus("weightLeakage",1);
   tJpsiLeak->SetBranchStatus("misPT", 1);

   //***********Set Binning

   RooBinning defaultMBins(floor((maxBMass-minBMass)/(40.)), B_plus_M.getMin(), B_plus_M.getMax() ); 
   RooBinning defaultMCorrBins(floor((10000-minBMass)/120.), B_plus_M_corr.getMin(), B_plus_M_corr.getMax()); 
   RooBinning broaderMBins(floor((maxBMass-minBMass)/(100.)), B_plus_M.getMin(), B_plus_M.getMax()); 
   RooBinning broaderMCorrBins(floor((10000-minBMass)/240.), B_plus_M_corr.getMin(), B_plus_M_corr.getMax()); 

   B_plus_M.setBinning( defaultMBins);
   B_plus_M_corr.setBinning( defaultMCorrBins );
   B_plus_M.setBinning( broaderMBins, "broaderBins");
   B_plus_M_corr.setBinning( broaderMCorrBins, "broaderBins" );

   B_plus_DTFM_M_zero.setBins(100);

   RooArgSet argset(BDTRooRealVar, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, misPT);
   RooArgSet argsetPartReco(BDTRooRealVar, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightPartReco, misPT);
   RooArgSet argsetLeakage(BDTRooRealVar, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, weightLeakage, misPT);
   RooArgSet argsetSignal(BDTRooRealVar, B_plus_DTFM_M_zero, B_plus_M_corr,  B_plus_M, trigVar, e_plus_BremMultiplicity, e_minus_BremMultiplicity, dataMCWeightee, misPT);

   cout<<"getting the datasets:"<<endl;


   string BDTCutString( ("("+BDTVar+">"+d2s(BDTcut)+")").c_str()  );
   string lowPPerpCut("(misPT<"+d2s(PPerpCut)+")");
   string highPPerpCut("(misPT>"+d2s(PPerpCut)+")");

   RooDataSet* dataSetSignalZeroGammaLowPPerp = new RooDataSet("dataSetSignalZeroGammaLowPPerp", "dataSetSignalZeroGammaLowPPerp", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+lowPPerpCut).c_str()), WeightVar("DataMCWeightee")  );
   RooDataSet* dataSetSignalOneGammaLowPPerp = new RooDataSet("dataSetSignalOneGammaLowPPerp", "dataSetSignalOneGammaLowPPerp", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+lowPPerpCut).c_str()), WeightVar("DataMCWeightee")  );
   RooDataSet* dataSetSignalTwoGammaLowPPerp = new RooDataSet("dataSetSignalTwoGammaLowPPerp", "dataSetSignalTwoGammaLowPPerp", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+lowPPerpCut).c_str()), WeightVar("DataMCWeightee")  );
   RooDataSet* dataSetPartRecoLowPPerp = new RooDataSet("dataSetPartRecoLowPPerp", "dataSetPartRecoLowPPerp",  argsetPartReco, Import(*tPartReco),Cut(("("+trigStr+"  > 0.9) && "+BDTCutString+ " && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+lowPPerpCut).c_str()), WeightVar("weightPartReco"));

   RooDataSet* dataSetJpsiLeakLowPPerp = new RooDataSet("dataSetJpsiLeakLowPPerp", "dataSetJpsiLeakLowPPerp",  argsetLeakage, Import(*tJpsiLeak),Cut(("B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+lowPPerpCut).c_str()), WeightVar("weightLeakage"));
   RooDataSet* dataSetCombLowPPerp = new RooDataSet("dataSetCombLowPPerp", "dataSetCombLowPPerp", tComb, argset, ("("+trigStr+"  > 0.9) && "+BDTCutString+"  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+lowPPerpCut).c_str());



   RooDataSet* dataSetSignalZeroGammaHighPPerp = new RooDataSet("dataSetSignalZeroGammaHighPPerp", "dataSetSignalZeroGammaHighPPerp", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > -0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 0.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+highPPerpCut).c_str()), WeightVar("DataMCWeightee")  );
   RooDataSet* dataSetSignalOneGammaHighPPerp = new RooDataSet("dataSetSignalOneGammaHighPPerp", "dataSetSignalOneGammaHighPPerp", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 0.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 1.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+highPPerpCut).c_str()), WeightVar("DataMCWeightee")  );
   RooDataSet* dataSetSignalTwoGammaHighPPerp = new RooDataSet("dataSetSignalTwoGammaHighPPerp", "dataSetSignalTwoGammaHighPPerp", argsetSignal, Import(*tSignal), Cut(( " ("+trigStr+"  > 0.9) && "+BDTCutString+" && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) > 1.5) && ((e_plus_BremMultiplicity+e_minus_BremMultiplicity) < 2.5) && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+highPPerpCut).c_str()), WeightVar("DataMCWeightee")  );
   RooDataSet* dataSetPartRecoHighPPerp = new RooDataSet("dataSetPartRecoHighPPerp", "dataSetPartRecoHighPPerp",  argsetPartReco, Import(*tPartReco),Cut(("("+trigStr+"  > 0.9) && "+BDTCutString+ " && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+highPPerpCut).c_str()), WeightVar("weightPartReco"));

   RooDataSet* dataSetJpsiLeakHighPPerp = new RooDataSet("dataSetJpsiLeakHighPPerp", "dataSetJpsiLeakHighPPerp",  argsetLeakage, Import(*tJpsiLeak),Cut(("B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+highPPerpCut).c_str()), WeightVar("weightLeakage"));
   RooDataSet* dataSetCombHighPPerp = new RooDataSet("dataSetCombHighPPerp", "dataSetCombHighPPerp", tComb, argset, ("("+trigStr+"  > 0.9) && "+BDTCutString+"  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)+" && "+highPPerpCut).c_str());


   RooDataSet* dataSetKemu = new RooDataSet("dataSetKemu", "dataSetKemu", tComb, argset, ("("+trigStr+"  > 0.9) && "+BDTCutString+"  && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());

   cout<<"Number of zero: LowPPerp:"<< dataSetSignalZeroGammaLowPPerp->sumEntries()<<", HighPPerp: "<<dataSetSignalZeroGammaHighPPerp->sumEntries()<<endl;
   cout<<"Number of one: LowPPerp:"<< dataSetSignalOneGammaLowPPerp->sumEntries()<<", HighPPerp: "<<dataSetSignalOneGammaHighPPerp->sumEntries()<<endl;
   cout<<"Number of two: LowPPerp:"<< dataSetSignalTwoGammaLowPPerp->sumEntries()<<", HighPPerp: "<<dataSetSignalTwoGammaHighPPerp->sumEntries()<<endl;
   cout<<"Number of PartReco: LowPPerp: "<< dataSetPartRecoLowPPerp->sumEntries()<<", HighPPerp: "<<dataSetPartRecoHighPPerp->sumEntries()<<endl;
   cout<<"Number of Jpsi leaking: LowPPerp: "<< dataSetJpsiLeakLowPPerp->sumEntries()<<", HighPPerp: "<<dataSetJpsiLeakHighPPerp->sumEntries()<<endl;
   cout<<"Number of combinatorial events: LowPPerp: "<< dataSetCombLowPPerp->sumEntries()<<", HighPPerp: "<<dataSetCombHighPPerp->sumEntries()<<endl;
   cout<<"Number of Kemu events: "<<dataSetKemu->sumEntries()<<endl;


   cout<<"binning the datasets:"<<endl;

   RooArgSet argset1D(B_plus_M);

   RooDataHist dataHistSignalZeroGammaLowPPerp("dataHistSignalZeroGammaLowPPerp", "dataHistSignalZeroGammaLowPPerp", argset1D, *dataSetSignalZeroGammaLowPPerp); 
   RooDataHist dataHistSignalOneGammaLowPPerp("dataHistSignalOneGammaLowPPerp", "dataHistSignalOneGammaLowPPerp", argset1D, *dataSetSignalOneGammaLowPPerp); 
   RooDataHist dataHistSignalTwoGammaLowPPerp("dataHistSignalTwoGammaLowPPerp", "dataHistSignalTwoGammaLowPPerp", argset1D, *dataSetSignalTwoGammaLowPPerp); 
   RooDataHist dataHistCombLowPPerp("dataHistCombLowPPerp", "dataHistCombLowPPerp", argset1D, *dataSetCombLowPPerp); 
   RooDataHist dataHistPartRecoLowPPerp("dataHistPartRecoLowPPerp", "dataHistPartRecoLowPPerp", argset1D, *dataSetPartRecoLowPPerp); 
   RooDataHist dataHistJpsiLeakLowPPerp("dataHistJpsiLeakLowPPerp", "dataHistJpsiLeakLowPPerp", argset1D);
   dataHistJpsiLeakLowPPerp.add(*dataSetJpsiLeakLowPPerp); 

   RooDataHist dataHistSignalZeroGammaHighPPerp("dataHistSignalZeroGammaHighPPerp", "dataHistSignalZeroGammaHighPPerp", argset1D, "broaderBins"); 
   dataHistSignalZeroGammaHighPPerp.add(*dataSetSignalZeroGammaHighPPerp);
   RooDataHist dataHistSignalOneGammaHighPPerp("dataHistSignalOneGammaHighPPerp", "dataHistSignalOneGammaHighPPerp", argset1D, "broaderBins"); 
   dataHistSignalOneGammaHighPPerp.add(*dataSetSignalOneGammaHighPPerp);
   RooDataHist dataHistSignalTwoGammaHighPPerp("dataHistSignalTwoGammaHighPPerp", "dataHistSignalTwoGammaHighPPerp", argset1D, "broaderBins"); 
   dataHistSignalTwoGammaHighPPerp.add(*dataSetSignalTwoGammaHighPPerp);
   RooDataHist dataHistCombHighPPerp("dataHistCombHighPPerp", "dataHistCombHighPPerp", argset1D, "broaderBins"); 
   dataHistCombHighPPerp.add(*dataSetCombHighPPerp);
   RooDataHist dataHistPartRecoHighPPerp("dataHistPartRecoHighPPerp", "dataHistPartRecoHighPPerp", argset1D, "broaderBins"); 
   dataHistPartRecoHighPPerp.add(*dataSetPartRecoHighPPerp);
   RooDataHist dataHistJpsiLeakHighPPerp("dataHistJpsiLeakHighPPerp", "dataHistJpsiLeakHighPPerp", argset1D);
   dataHistJpsiLeakHighPPerp.add(*dataSetJpsiLeakHighPPerp); 

   //RooDataHist dataHistKemu("dataHistKemu", "dataHistKemu", argset2Kemu, *dataSetKemu);

   //*************** Compute Error on J/psi leak

   double ErrorJpsi(0);
   double nJpsi( dataSetJpsiLeakLowPPerp->sumEntries(("("+trigStr+"  > 0.9) && "+BDTCutString+" && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str()) );
   nJpsi+=dataSetJpsiLeakHighPPerp->sumEntries(("("+trigStr+"  > 0.9) && "+BDTCutString+" && B_plus_M > "+d2s(minBMass)+" && B_plus_M < "+d2s(maxBMass)).c_str());

   if(nJpsi > 0) ErrorJpsi = 1./sqrt(nJpsi);
   RooRealVar fractionalErrorJpsiLeak("fractionalErrorJpsiLeak", "fractionalErrorJpsiLeak", ErrorJpsi);
   cout<<"JPSI LEAK: "<<nJpsi<<endl;
   cout<<"JPSI LEAK fractional Error: "<<ErrorJpsi<<endl;


   //***************Create 2D histogram estimates from data


   cout<<"Preparing the  histPdfs: 1";
   RooHistPdf histPdfSignalZeroGammaLowPPerp("histPdfSignalZeroGammaLowPPerp", "histPdfSignalZeroGammaLowPPerp", argset1D, dataHistSignalZeroGammaLowPPerp,2); cout<<" 2";
   RooHistPdf histPdfSignalOneGammaLowPPerp("histPdfSignalOneGammaLowPPerp", "histPdfSignalOneGammaLowPPerp", argset1D, dataHistSignalOneGammaLowPPerp,2); cout<<" 3";
   RooHistPdf histPdfSignalTwoGammaLowPPerp("histPdfSignalTwoGammaLowPPerp", "histPdfSignalTwoGammaLowPPerp", argset1D, dataHistSignalTwoGammaLowPPerp,2); cout<<" 4";
   RooHistPdf histPdfPartRecoLowPPerp("histPdfPartRecoLowPPerp", "histPdfPartRecoLowPPerp", argset1D, dataHistPartRecoLowPPerp,2); cout<<" 5";
   RooHistPdf histPdfJpsiLeakLowPPerp("histPdfJpsiLeakLowPPerp", "histPdfJpsiLeakLowPPerp", argset1D, dataHistJpsiLeakLowPPerp,2); cout<<" 6";

   RooHistPdf histPdfSignalZeroGammaHighPPerp("histPdfSignalZeroGammaHighPPerp", "histPdfSignalZeroGammaHighPPerp", argset1D, dataHistSignalZeroGammaHighPPerp,2); cout<<" 7";
   RooHistPdf histPdfSignalOneGammaHighPPerp("histPdfSignalOneGammaHighPPerp", "histPdfSignalOneGammaHighPPerp", argset1D, dataHistSignalOneGammaHighPPerp,2); cout<<" 8";
   RooHistPdf histPdfSignalTwoGammaHighPPerp("histPdfSignalTwoGammaHighPPerp", "histPdfSignalTwoGammaHighPPerp", argset1D, dataHistSignalTwoGammaHighPPerp,2); cout<<" 9";
   RooHistPdf histPdfPartRecoHighPPerp("histPdfPartRecoHighPPerp", "histPdfPartRecoHighPPerp", argset1D, dataHistPartRecoHighPPerp,2); cout<<" 10";
   RooHistPdf histPdfJpsiLeakHighPPerp("histPdfJpsiLeakHighPPerp", "histPdfJpsiLeakHighPPerp", argset1D, dataHistJpsiLeakHighPPerp,2); cout<<" 11";

   //RooHistPdf histPdfKemu("histPdfKemu", "histPdfKemu", argset2Kemu, dataHistKemu,2);

   //***************Create combinatorial shape from fit to data

   RooRealVar expoConstCombLowPPerp("expoConstCombLowPPerp", "expoConstCombLowPPerp", -1e-3, -1, 1);
   RooRealVar expoConstCombHighPPerp("expoConstCombHighPPerp", "expoConstCombHighPPerp", -1e-3, -1, 1);

   RooRealVar expoConstKemu("expoConstKemu", "expoConstKemu", -1e-3, -1, 1);
   RooRealVar T("T", "T", 30, 0, 200);
   RooRealVar n("n", "n", 3.5, 2, 10);

   RooExponential pdfCombLowPPerp("pdfCombLowPPerp", "histPdfCombLowPPerp", B_plus_M, expoConstCombLowPPerp);
   pdfCombLowPPerp.fitTo(*dataSetCombLowPPerp); 

   RooExponential pdfCombHighPPerp("pdfCombHighPPerp", "histPdfCombHighPPerp", B_plus_M, expoConstCombHighPPerp);
   pdfCombHighPPerp.fitTo(*dataSetCombHighPPerp); 

   RooPTMVis pdfKemu("pdfKemu", "pdfKemu", misPT, B_plus_M, T, n, expoConstKemu);
   RooFitResult* rfr = pdfKemu.fitTo(*dataSetKemu, RooFit::Save());


   RooRealVar trueExpoConstCombLowPPerp("trueExpoConstCombLowPPerp", "trueExpoConstCombLowPPerp", expoConstCombLowPPerp.getVal(), -1000, 1000); 
   trueExpoConstCombLowPPerp.setError(expoConstCombLowPPerp.getError());
   RooRealVar trueExpoConstCombHighPPerp("trueExpoConstCombHighPPerp", "trueExpoConstCombHighPPerp", expoConstCombHighPPerp.getVal(), -1000, 1000); 
   trueExpoConstCombHighPPerp.setError(expoConstCombHighPPerp.getError());
   RooRealVar trueExpoConstKemu("trueExpoConstKemu", "trueExpoConstKemu", expoConstKemu.getVal(), -1000, 1000);
   trueExpoConstKemu.setError(expoConstKemu.getError());

   RooRealVar trueT("trueT","trueT", T.getVal(), -1000, 1000);
   trueT.setError(T.getError());
   RooRealVar trueN("trueN","trueN", n.getVal(), -1000, 1000);
   trueN.setError(n.getError());

   trueExpoConstCombLowPPerp.setConstant(true);
   trueExpoConstCombHighPPerp.setConstant(true);
   trueExpoConstKemu.setConstant(true);
   trueT.setConstant(true);
   trueN.setConstant(true);

   //***************Create integral to get the fractions for the fit


   RooAbsReal* igLowPPerp = pdfKemu.createIntegral(misPT,RooFit::NormSet(misPT),RooFit::Range("lowPPerp"));
   RooAbsReal* igAllPPerp = pdfKemu.createIntegral(misPT,RooFit::NormSet(misPT),RooFit::Range("default"));
   RooFormulaVar FracCombFormulaLowOverTot("FracCombFormulaLowOverTot", "FracCombFormulaLowOverTot", "@0/@1" , RooArgList(*igLowPPerp, *igAllPPerp));

   FracCombRealVarLowOverTot.setVal(FracCombFormulaLowOverTot.getVal());
   FracCombRealVarLowOverTot.setError(FracCombFormulaLowOverTot.getPropagatedError(*rfr));


   //***************Prepare the different fractions

   nGenSignalZeroGamma = TMath::Nint(nGenFracZeroGamma*nGenSignal);
   nGenSignalOneGamma = TMath::Nint(nGenFracOneGamma*nGenSignal);
   nGenSignalTwoGamma = nGenSignal-nGenSignalOneGamma-nGenSignalZeroGamma ;

   nGenFracSigZeroGammaLowOverTot =  dataSetSignalZeroGammaLowPPerp->sumEntries() / (dataSetSignalZeroGammaLowPPerp->sumEntries() + dataSetSignalZeroGammaHighPPerp->sumEntries() );
   nGenFracSigOneGammaLowOverTot = dataSetSignalOneGammaLowPPerp->sumEntries() / (dataSetSignalOneGammaLowPPerp->sumEntries() + dataSetSignalOneGammaHighPPerp->sumEntries() );
   nGenFracSigTwoGammaLowOverTot = dataSetSignalTwoGammaLowPPerp->sumEntries() / (dataSetSignalTwoGammaLowPPerp->sumEntries() + dataSetSignalTwoGammaHighPPerp->sumEntries() );

   nGenFracSigLowOverTot = (dataSetSignalZeroGammaLowPPerp->sumEntries() + dataSetSignalOneGammaLowPPerp->sumEntries() + dataSetSignalTwoGammaLowPPerp->sumEntries()) / 
      (dataSetSignalZeroGammaLowPPerp->sumEntries() + dataSetSignalZeroGammaHighPPerp->sumEntries()  + dataSetSignalOneGammaLowPPerp->sumEntries() + dataSetSignalOneGammaHighPPerp->sumEntries() + dataSetSignalTwoGammaLowPPerp->sumEntries() + dataSetSignalTwoGammaHighPPerp->sumEntries() );

   nGenSignalZeroGammaLowPPerp = TMath::Nint( nGenFracSigZeroGammaLowOverTot * nGenSignalZeroGamma );
   nGenSignalOneGammaLowPPerp = TMath::Nint( nGenFracSigOneGammaLowOverTot * nGenSignalOneGamma );
   nGenSignalTwoGammaLowPPerp = TMath::Nint( nGenFracSigTwoGammaLowOverTot * nGenSignalTwoGamma ); 
   nGenSignalLowPPerp = nGenSignalZeroGammaLowPPerp + nGenSignalOneGammaLowPPerp + nGenSignalTwoGammaLowPPerp;

   nGenSignalZeroGammaHighPPerp = nGenSignalZeroGamma - nGenSignalZeroGammaLowPPerp;
   nGenSignalOneGammaHighPPerp = nGenSignalOneGamma - nGenSignalOneGammaLowPPerp;
   nGenSignalTwoGammaHighPPerp = nGenSignalTwoGamma - nGenSignalTwoGammaLowPPerp;
   nGenSignalHighPPerp = nGenSignal - nGenSignalLowPPerp;

   nGenFracZeroGammaLowPPerp = nGenSignalZeroGammaLowPPerp / (1.*nGenSignalLowPPerp);
   nGenFracOneGammaLowPPerp = nGenSignalOneGammaLowPPerp / (1.*nGenSignalLowPPerp);
   nGenFracZeroGammaHighPPerp = nGenSignalZeroGammaHighPPerp / (1.*nGenSignalHighPPerp);
   nGenFracOneGammaHighPPerp = nGenSignalOneGammaHighPPerp / (1.*nGenSignalHighPPerp);

   nGenFracCombLowOverTot = dataSetCombLowPPerp->sumEntries() / (dataSetCombLowPPerp->sumEntries() + dataSetCombHighPPerp->sumEntries());
   nGenFracPartRecoLowOverTot = dataSetPartRecoLowPPerp->sumEntries() / (dataSetPartRecoLowPPerp->sumEntries() + dataSetPartRecoHighPPerp->sumEntries());
   nGenFracJpsiLeakLowOverTot = dataSetJpsiLeakLowPPerp->sumEntries() / ( dataSetJpsiLeakLowPPerp->sumEntries() + dataSetJpsiLeakHighPPerp->sumEntries()  );

   nGenCombLowPPerp = TMath::Nint( nGenFracCombLowOverTot * nGenComb );
   nGenPartRecoLowPPerp = TMath::Nint( nGenFracPartRecoLowOverTot * nGenPartReco );
   nGenJpsiLeakLowPPerp = TMath::Nint( nGenFracJpsiLeakLowOverTot * nGenJpsiLeak );

   nGenCombHighPPerp = nGenComb - nGenCombLowPPerp ;
   nGenPartRecoHighPPerp = nGenPartReco - nGenPartRecoLowPPerp;
   nGenJpsiLeakHighPPerp = nGenJpsiLeak - nGenJpsiLeakLowPPerp;



   FracCombFromFitLowOverTot = FracCombRealVarLowOverTot.getVal();
   nCombFromFitLowPPerp = TMath::Nint(FracCombFromFitLowOverTot*nGenComb);
   nCombFromFitHighPerp = nGenComb - nCombFromFitLowPPerp;




   //***************Save everything on a workspace
   RooWorkspace workspace("workspace", "workspace");
   workspace.import(B_plus_DTFM_M_zero);
   workspace.import(B_plus_M);
   workspace.import(B_plus_M_corr);
   workspace.import(misPT);
   workspace.import(expoConstCombLowPPerp);
   workspace.import(expoConstCombHighPPerp);
   workspace.import(expoConstKemu);
   workspace.import(T);
   workspace.import(n);
   workspace.import(trueExpoConstCombLowPPerp);
   workspace.import(trueExpoConstCombHighPPerp);
   workspace.import(trueExpoConstKemu);
   workspace.import(trueT);
   workspace.import(trueN);

   workspace.import(*dataSetKemu);

   workspace.import(*dataSetSignalZeroGammaHighPPerp);
   workspace.import(*dataSetSignalOneGammaHighPPerp);
   workspace.import(*dataSetSignalTwoGammaHighPPerp);
   workspace.import(*dataSetPartRecoHighPPerp);
   workspace.import(*dataSetCombHighPPerp);
   workspace.import(*dataSetJpsiLeakHighPPerp);

   workspace.import(*dataSetSignalZeroGammaLowPPerp);
   workspace.import(*dataSetSignalOneGammaLowPPerp);
   workspace.import(*dataSetSignalTwoGammaLowPPerp);
   workspace.import(*dataSetPartRecoLowPPerp);
   workspace.import(*dataSetCombLowPPerp);
   workspace.import(*dataSetJpsiLeakLowPPerp);

   workspace.import(histPdfSignalZeroGammaHighPPerp);
   workspace.import(histPdfSignalOneGammaHighPPerp);
   workspace.import(histPdfSignalTwoGammaHighPPerp);
   workspace.import(histPdfPartRecoHighPPerp);
   workspace.import(histPdfJpsiLeakHighPPerp);

   workspace.import(histPdfSignalZeroGammaLowPPerp);
   workspace.import(histPdfSignalOneGammaLowPPerp);
   workspace.import(histPdfSignalTwoGammaLowPPerp);
   workspace.import(histPdfPartRecoLowPPerp);
   workspace.import(histPdfJpsiLeakLowPPerp);

   workspace.import(fractionalErrorJpsiLeak);

   workspace.writeToFile(workspacename.c_str());

   //close and erase everything

   fComb->Close();
   fSignal->Close();
   fPartReco->Close();

   delete fComb;
   delete fSignal;
   delete fPartReco;

   delete dataSetSignalZeroGammaHighPPerp; 
   delete dataSetSignalOneGammaHighPPerp;
   delete dataSetSignalTwoGammaHighPPerp;
   delete dataSetPartRecoHighPPerp;
   delete dataSetJpsiLeakHighPPerp;
   delete dataSetCombHighPPerp;

   delete dataSetSignalZeroGammaLowPPerp; 
   delete dataSetSignalOneGammaLowPPerp;
   delete dataSetSignalTwoGammaLowPPerp;
   delete dataSetPartRecoLowPPerp;
   delete dataSetJpsiLeakLowPPerp;
   delete dataSetCombLowPPerp;

   delete dataSetKemu;

   delete rfr;
   delete igLowPPerp;
   delete igAllPPerp;


   cout<<"Workspace ready!"<<endl;
}

void FitterUtils1DSimultaneous::initiateParamsKemu(RooRealVar& expoConstKemu, RooRealVar const& trueExpoConstKemu, RooRealVar& T, RooRealVar const& trueT, RooRealVar& n, RooRealVar const& trueN)
{
   TRandom rand;
   rand.SetSeed();

   T.setVal(rand.Uniform( TMath::Max(0.,trueT.getVal() - 5*trueT.getError()), trueT.getVal() + 5*trueT.getError()));
   T.setRange(TMath::Max(trueT.getVal() - 10*trueT.getError(),0.), trueT.getVal() + 10*trueT.getError());

   n.setVal(rand.Uniform(TMath::Max( trueN.getVal() - 5*trueN.getError(),2.), trueN.getVal() + 5*trueN.getError()));
   n.setRange(TMath::Max(trueN.getVal() - 10*trueN.getError(),2.), trueN.getVal() + 10*trueN.getError()); 

   expoConstKemu.setVal(rand.Uniform( trueExpoConstKemu.getVal() - 5*trueExpoConstKemu.getError(), trueExpoConstKemu.getVal() + 5*trueExpoConstKemu.getError() ) );
   expoConstKemu.setRange( trueExpoConstKemu.getVal() - 10*trueExpoConstKemu.getError(), trueExpoConstKemu.getVal() + 10*trueExpoConstKemu.getError() );  
}


void FitterUtils1DSimultaneous::initiateParamsKee(RooRealVar& nSignal, RooRealVar& nPartReco, RooRealVar& nComb, RooRealVar& nJpsiLeak,
      RooRealVar& fracZeroGammaLowPPerp, RooRealVar const& fracZeroGammaLowPPerpSigma, RooRealVar& fracOneGammaLowPPerp, RooRealVar const& fracOneGammaLowPPerpSigma, 
      RooRealVar& fracZeroGammaHighPPerp, RooRealVar const& fracZeroGammaHighPPerpSigma, RooRealVar& fracOneGammaHighPPerp, RooRealVar const& fracOneGammaHighPPerpSigma, 
      RooRealVar& expoConstCombLowPPerp, RooRealVar const& trueExpoConstCombLowPPerp, RooRealVar& expoConstCombHighPPerp, RooRealVar const& trueExpoConstCombHighPPerp,
      RooRealVar& fracCombLowOverTot, RooRealVar const& fracCombLowOverTotFromKemuFit)
{
   TRandom rand;
   rand.SetSeed();

   double nGenSignal2 = rand.Uniform(nGenSignal-5*sqrt(nGenSignal), nGenSignal+5*sqrt(nGenSignal));
   double nGenPartReco2 = rand.Uniform(TMath::Max(0.,nGenPartReco-5*sqrt(nGenPartReco)), nGenPartReco+5*sqrt(nGenPartReco));
   double nGenComb2 = rand.Uniform(nGenComb-5*sqrt(nGenComb), nGenComb+5*sqrt(nGenComb));
   double nGenJpsiLeak2 = rand.Uniform(TMath::Max(0.,nGenJpsiLeak-5*sqrt(nGenJpsiLeak)), nGenJpsiLeak+5*sqrt(nGenJpsiLeak));

   nSignal.setVal(nGenSignal2);
   nSignal.setRange(TMath::Max(0.,nGenSignal2-10.*sqrt(nGenSignal)) , nGenSignal2+10*sqrt(nGenSignal));
   nPartReco.setVal(nGenPartReco2);
   nPartReco.setRange(TMath::Max(0.,nGenPartReco2-10.*sqrt(nGenPartReco)), nGenPartReco2+10*sqrt(nGenPartReco));
   nComb.setVal(nGenComb2);
   nComb.setRange(TMath::Max(0.,nGenComb2-10.*sqrt(nGenComb)), nGenComb2+10*sqrt(nGenComb));
   nJpsiLeak.setVal(nGenJpsiLeak2);
   nJpsiLeak.setRange(TMath::Max(0., nGenJpsiLeak2-10*sqrt(nGenJpsiLeak)), nGenJpsiLeak2+10*sqrt(nGenJpsiLeak));

   expoConstCombLowPPerp.setVal(rand.Uniform( trueExpoConstCombLowPPerp.getVal() - 5*trueExpoConstCombLowPPerp.getError(), trueExpoConstCombLowPPerp.getVal() + 5*trueExpoConstCombLowPPerp.getError() ) );
   expoConstCombLowPPerp.setRange( trueExpoConstCombLowPPerp.getVal() - 10*trueExpoConstCombLowPPerp.getError(), trueExpoConstCombLowPPerp.getVal() + 10*trueExpoConstCombLowPPerp.getError() );  

   expoConstCombHighPPerp.setVal(rand.Uniform( trueExpoConstCombHighPPerp.getVal() - 5*trueExpoConstCombHighPPerp.getError(), trueExpoConstCombHighPPerp.getVal() + 5*trueExpoConstCombHighPPerp.getError() ) );
   expoConstCombHighPPerp.setRange( trueExpoConstCombHighPPerp.getVal() - 10*trueExpoConstCombHighPPerp.getError(), trueExpoConstCombHighPPerp.getVal() + 10*trueExpoConstCombHighPPerp.getError() );  
  
   fracCombLowOverTot.setVal(rand.Uniform(TMath::Max(0.,fracCombLowOverTotFromKemuFit.getVal() - 5*fracCombLowOverTotFromKemuFit.getError()), TMath::Min(0.99999999,fracCombLowOverTotFromKemuFit.getVal() + 5*fracCombLowOverTotFromKemuFit.getError())));
   fracCombLowOverTot.setRange(fracCombLowOverTotFromKemuFit.getVal() - 10*fracCombLowOverTotFromKemuFit.getError(), fracCombLowOverTotFromKemuFit.getVal() + 10*fracCombLowOverTotFromKemuFit.getError());

   fracZeroGammaLowPPerp.setRange(nGenFracZeroGammaLowPPerp-10*fracZeroGammaLowPPerpSigma.getVal(), nGenFracZeroGammaLowPPerp+10*fracZeroGammaLowPPerpSigma.getVal());
   fracZeroGammaLowPPerp.setVal(rand.Gaus(nGenFracZeroGammaLowPPerp, fracZeroGammaLowPPerpSigma.getVal()));
   fracOneGammaLowPPerp.setRange(nGenFracOneGammaLowPPerp-10*fracOneGammaLowPPerpSigma.getVal(), nGenFracOneGammaLowPPerp+10*fracOneGammaLowPPerpSigma.getVal());
   fracOneGammaLowPPerp.setVal(rand.Gaus(nGenFracOneGammaLowPPerp, fracOneGammaLowPPerpSigma.getVal()));


   fracZeroGammaHighPPerp.setRange(nGenFracZeroGammaHighPPerp-10*fracZeroGammaHighPPerpSigma.getVal(), nGenFracZeroGammaHighPPerp+10*fracZeroGammaHighPPerpSigma.getVal());
   fracZeroGammaHighPPerp.setVal(rand.Gaus(nGenFracZeroGammaHighPPerp, fracZeroGammaHighPPerpSigma.getVal()));
   fracOneGammaHighPPerp.setRange(nGenFracOneGammaHighPPerp-10*fracOneGammaHighPPerpSigma.getVal(), nGenFracOneGammaHighPPerp+10*fracOneGammaHighPPerpSigma.getVal());
   fracOneGammaHighPPerp.setVal(rand.Gaus(nGenFracOneGammaHighPPerp, fracOneGammaHighPPerpSigma.getVal()));
}

void FitterUtils1DSimultaneous::generate()
{
   //***************Get the PDFs from the workspace

   TFile fw(workspacename.c_str(), "UPDATE");   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooRealVar* trueT = workspace->var("trueT");
   RooRealVar *trueN = workspace->var("trueN");
   RooRealVar *trueExpoConstCombLowPPerp = workspace->var("trueExpoConstCombLowPPerp");
   RooRealVar *trueExpoConstCombHighPPerp = workspace->var("trueExpoConstCombHighPPerp");
   RooRealVar* trueExpoConstKemu = workspace->var("trueExpoConstKemu");
   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");


   RooHistPdf *histPdfSignalZeroGammaHighPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGammaHighPPerp");
   RooHistPdf *histPdfSignalOneGammaHighPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalOneGammaHighPPerp");
   RooHistPdf *histPdfSignalTwoGammaHighPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGammaHighPPerp");
   RooHistPdf *histPdfPartRecoHighPPerp = (RooHistPdf *) workspace->pdf("histPdfPartRecoHighPPerp");
   RooHistPdf *histPdfJpsiLeakHighPPerp(0);
   if(nGenJpsiLeak>1) histPdfJpsiLeakHighPPerp = (RooHistPdf *) workspace->pdf("histPdfJpsiLeakHighPPerp");

   RooHistPdf *histPdfSignalZeroGammaLowPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGammaLowPPerp");
   RooHistPdf *histPdfSignalOneGammaLowPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalOneGammaLowPPerp");
   RooHistPdf *histPdfSignalTwoGammaLowPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGammaLowPPerp");
   RooHistPdf *histPdfPartRecoLowPPerp = (RooHistPdf *) workspace->pdf("histPdfPartRecoLowPPerp");
   RooHistPdf *histPdfJpsiLeakLowPPerp(0);
   if(nGenJpsiLeak>1) histPdfJpsiLeakLowPPerp = (RooHistPdf *) workspace->pdf("histPdfJpsiLeakLowPPerp");


   RooExponential* pdfCombLowPPerp = new RooExponential("pdfCombLowPPerp", "pdfCombLowPPerp", *B_plus_M, *trueExpoConstCombLowPPerp);
   RooExponential* pdfCombHighPPerp = new RooExponential("pdfCombHighPPerp", "pdfCombHighPPerp", *B_plus_M, *trueExpoConstCombHighPPerp);

   RooPTMVis* pdfKemu =  new RooPTMVis("pdfKemu", "pdfKemu", *misPT, *B_plus_M, *trueT, *trueN, *trueExpoConstKemu);


   //***************Prepare the GenSpecs

   RooArgSet argset1D(*B_plus_M);

   RooRandom::randomGenerator()->SetSeed();
   RooAbsPdf::GenSpec* GenSpecSignalZeroGammaLowPPerp = histPdfSignalZeroGammaLowPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenSignalZeroGammaLowPPerp));
   RooAbsPdf::GenSpec* GenSpecSignalOneGammaLowPPerp = histPdfSignalOneGammaLowPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenSignalOneGammaLowPPerp));
   RooAbsPdf::GenSpec* GenSpecSignalTwoGammaLowPPerp = histPdfSignalTwoGammaLowPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenSignalTwoGammaLowPPerp));
   RooAbsPdf::GenSpec* GenSpecPartRecoLowPPerp =  histPdfPartRecoLowPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenPartRecoLowPPerp));
   RooAbsPdf::GenSpec* GenSpecCombLowPPerp = pdfCombLowPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenCombLowPPerp));
   RooAbsPdf::GenSpec* GenSpecJpsiLeakLowPPerp(0);
   if(nGenJpsiLeak>1) GenSpecJpsiLeakLowPPerp = histPdfJpsiLeakLowPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenJpsiLeakLowPPerp));

   RooAbsPdf::GenSpec* GenSpecSignalZeroGammaHighPPerp = histPdfSignalZeroGammaHighPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenSignalZeroGammaHighPPerp));
   RooAbsPdf::GenSpec* GenSpecSignalOneGammaHighPPerp = histPdfSignalOneGammaHighPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenSignalOneGammaHighPPerp));
   RooAbsPdf::GenSpec* GenSpecSignalTwoGammaHighPPerp = histPdfSignalTwoGammaHighPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenSignalTwoGammaHighPPerp));
   RooAbsPdf::GenSpec* GenSpecPartRecoHighPPerp =  histPdfPartRecoHighPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenPartRecoHighPPerp));
   RooAbsPdf::GenSpec* GenSpecCombHighPPerp = pdfCombHighPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenCombHighPPerp));
   RooAbsPdf::GenSpec* GenSpecJpsiLeakHighPPerp(0);
   if(nGenJpsiLeak>1) GenSpecJpsiLeakHighPPerp = histPdfJpsiLeakHighPPerp->prepareMultiGen(argset1D, RooFit::Extended(1), NumEvents(nGenJpsiLeakHighPPerp));

   RooArgSet argsetKemu(*misPT, *B_plus_M);

   RooAbsPdf::GenSpec* GenSpecKemu = pdfKemu->prepareMultiGen(argsetKemu, RooFit::Extended(1), NumEvents(nGenKemu));
   

   //***************Generate some datasets

   cout<<"Generating signal Zero Photon HighPPerp"<<endl;
   RooDataSet* dataGenSignalZeroGammaHighPPerp = histPdfSignalZeroGammaHighPPerp->generate(*GenSpecSignalZeroGammaHighPPerp);//(argset, 250, false, true, "", false, true);
   dataGenSignalZeroGammaHighPPerp->SetName("dataGenSignalZeroGammaHighPPerp"); dataGenSignalZeroGammaHighPPerp->SetTitle("dataGenSignalZeroGammaHighPPerp");
   cout<<"Generating signal One Photon HighPPerp"<<endl;
   RooDataSet* dataGenSignalOneGammaHighPPerp = histPdfSignalOneGammaHighPPerp->generate(*GenSpecSignalOneGammaHighPPerp);//(argset, 250, false, true, "", false, true);
   dataGenSignalOneGammaHighPPerp->SetName("dataGenSignalOneGammaHighPPerp"); dataGenSignalOneGammaHighPPerp->SetTitle("dataGenSignalOneGammaHighPPerp");
   cout<<"Generating signal two Photons HighPPerp"<<endl;
   RooDataSet* dataGenSignalTwoGammaHighPPerp = histPdfSignalTwoGammaHighPPerp->generate(*GenSpecSignalTwoGammaHighPPerp);//(argset, 250, false, true, "", false, true);
   dataGenSignalTwoGammaHighPPerp->SetName("dataGenSignalTwoGammaHighPPerp"); dataGenSignalTwoGammaHighPPerp->SetTitle("dataGenSignalTwoGammaHighPPerp");
   cout<<"Generating combinatorial"<<endl;
   RooDataSet* dataGenCombHighPPerp = pdfCombHighPPerp->generate(*GenSpecCombHighPPerp);//(argset, 100, false, true, "", false, true);
   dataGenCombHighPPerp->SetName("dataGenCombHighPPerp"); dataGenCombHighPPerp->SetTitle("dataGenCombHighPPerp");
   cout<<"Generating PartReco HighPPerp"<<endl;
   RooDataSet* dataGenPartRecoHighPPerp = histPdfPartRecoHighPPerp->generate(*GenSpecPartRecoHighPPerp);//argset, 160, false, true, "", false, true);
   dataGenPartRecoHighPPerp->SetName("dataGenPartRecoHighPPerp"); dataGenPartRecoHighPPerp->SetTitle("dataGenPartRecoHighPPerp");
   RooDataSet* dataGenJpsiLeakHighPPerp(0);
   if(nGenJpsiLeak>0)
   {
      cout<<"Generating Leaking JPsi"<<endl;
      dataGenJpsiLeakHighPPerp = histPdfJpsiLeakHighPPerp->generate(*GenSpecJpsiLeakHighPPerp);//argset, 160, false, true, "", false, true);
      dataGenJpsiLeakHighPPerp->SetName("dataGenJpsiLeakHighPPerp"); dataGenJpsiLeakHighPPerp->SetTitle("dataGenJpsiLeakHighPPerp");
   }

   cout<<"Generating signal Zero Photon LowPPerp"<<endl;
   RooDataSet* dataGenSignalZeroGammaLowPPerp = histPdfSignalZeroGammaLowPPerp->generate(*GenSpecSignalZeroGammaLowPPerp);//(argset, 250, false, true, "", false, true);
   dataGenSignalZeroGammaLowPPerp->SetName("dataGenSignalZeroGammaLowPPerp"); dataGenSignalZeroGammaLowPPerp->SetTitle("dataGenSignalZeroGammaLowPPerp");
   cout<<"Generating signal One Photon LowPPerp"<<endl;
   RooDataSet* dataGenSignalOneGammaLowPPerp = histPdfSignalOneGammaLowPPerp->generate(*GenSpecSignalOneGammaLowPPerp);//(argset, 250, false, true, "", false, true);
   dataGenSignalOneGammaLowPPerp->SetName("dataGenSignalOneGammaLowPPerp"); dataGenSignalOneGammaLowPPerp->SetTitle("dataGenSignalOneGammaLowPPerp");
   cout<<"Generating signal two Photons LowPPerp"<<endl;
   RooDataSet* dataGenSignalTwoGammaLowPPerp = histPdfSignalTwoGammaLowPPerp->generate(*GenSpecSignalTwoGammaLowPPerp);//(argset, 250, false, true, "", false, true);
   dataGenSignalTwoGammaLowPPerp->SetName("dataGenSignalTwoGammaLowPPerp"); dataGenSignalTwoGammaLowPPerp->SetTitle("dataGenSignalTwoGammaLowPPerp");
   cout<<"Generating combinatorial"<<endl;
   RooDataSet* dataGenCombLowPPerp = pdfCombLowPPerp->generate(*GenSpecCombLowPPerp);//(argset, 100, false, true, "", false, true);
   dataGenCombLowPPerp->SetName("dataGenCombLowPPerp"); dataGenCombLowPPerp->SetTitle("dataGenCombLowPPerp");
   cout<<"Generating PartReco LowPPerp"<<endl;
   RooDataSet* dataGenPartRecoLowPPerp = histPdfPartRecoLowPPerp->generate(*GenSpecPartRecoLowPPerp);//argset, 160, false, true, "", false, true);
   dataGenPartRecoLowPPerp->SetName("dataGenPartRecoLowPPerp"); dataGenPartRecoLowPPerp->SetTitle("dataGenPartRecoLowPPerp");
   RooDataSet* dataGenJpsiLeakLowPPerp(0);
   if(nGenJpsiLeak>0)
   {
      cout<<"Generating Leaking JPsi"<<endl;
      dataGenJpsiLeakLowPPerp = histPdfJpsiLeakLowPPerp->generate(*GenSpecJpsiLeakLowPPerp);//argset, 160, false, true, "", false, true);
      dataGenJpsiLeakLowPPerp->SetName("dataGenJpsiLeakLowPPerp"); dataGenJpsiLeakLowPPerp->SetTitle("dataGenJpsiLeakLowPPerp");
   }

   RooDataSet* dataGenKemu = pdfKemu->generate(*GenSpecKemu);
   dataGenKemu->SetName("dataGenKemu"); dataGenKemu->SetName("dataGenKemu");

   //*************Saving the generated datasets in a workspace

   RooWorkspace workspaceGen("workspaceGen", "workspaceGen");

   workspaceGen.import(*dataGenSignalZeroGammaHighPPerp);
   workspaceGen.import(*dataGenSignalOneGammaHighPPerp);
   workspaceGen.import(*dataGenSignalTwoGammaHighPPerp);
   workspaceGen.import(*dataGenCombHighPPerp);
   workspaceGen.import(*dataGenPartRecoHighPPerp);
   if(nGenJpsiLeak>1) workspaceGen.import(*dataGenJpsiLeakHighPPerp);

   workspaceGen.import(*dataGenSignalZeroGammaLowPPerp);
   workspaceGen.import(*dataGenSignalOneGammaLowPPerp);
   workspaceGen.import(*dataGenSignalTwoGammaLowPPerp);
   workspaceGen.import(*dataGenCombLowPPerp);
   workspaceGen.import(*dataGenPartRecoLowPPerp);
   if(nGenJpsiLeak>1) workspaceGen.import(*dataGenJpsiLeakLowPPerp);
    
   workspaceGen.import(*dataGenKemu);
   
   workspaceGen.Write("", TObject::kOverwrite);


   //delete workspace;
   fw.Close();


   delete dataGenSignalZeroGammaHighPPerp;
   delete dataGenSignalOneGammaHighPPerp;
   delete dataGenSignalTwoGammaHighPPerp;
   delete dataGenCombHighPPerp;
   delete dataGenPartRecoHighPPerp;
   if(nGenJpsiLeak>0) delete dataGenJpsiLeakHighPPerp;

   delete GenSpecSignalZeroGammaHighPPerp;
   delete GenSpecSignalOneGammaHighPPerp;
   delete GenSpecSignalTwoGammaHighPPerp;
   delete GenSpecCombHighPPerp;
   delete GenSpecPartRecoHighPPerp;
   if(nGenJpsiLeak>0) delete GenSpecJpsiLeakHighPPerp;

   delete histPdfSignalZeroGammaHighPPerp; 
   delete histPdfSignalOneGammaHighPPerp;
   delete histPdfSignalTwoGammaHighPPerp;
   delete histPdfPartRecoHighPPerp;
   if(nGenJpsiLeak>0) delete histPdfJpsiLeakHighPPerp;

   delete dataGenSignalZeroGammaLowPPerp;
   delete dataGenSignalOneGammaLowPPerp;
   delete dataGenSignalTwoGammaLowPPerp;
   delete dataGenCombLowPPerp;
   delete dataGenPartRecoLowPPerp;
   if(nGenJpsiLeak>0) delete dataGenJpsiLeakLowPPerp;

   delete GenSpecSignalZeroGammaLowPPerp;
   delete GenSpecSignalOneGammaLowPPerp;
   delete GenSpecSignalTwoGammaLowPPerp;
   delete GenSpecCombLowPPerp;
   delete GenSpecPartRecoLowPPerp;
   if(nGenJpsiLeak>0) delete GenSpecJpsiLeakLowPPerp;

   delete histPdfSignalZeroGammaLowPPerp; 
   delete histPdfSignalOneGammaLowPPerp;
   delete histPdfSignalTwoGammaLowPPerp;
   delete histPdfPartRecoLowPPerp;
   if(nGenJpsiLeak>0) delete histPdfJpsiLeakLowPPerp;

   delete pdfCombLowPPerp;
   delete pdfCombHighPPerp; 
   delete pdfKemu;

}



void FitterUtils1DSimultaneous::fit(bool wantplot, bool constPartReco,
      double fracPartReco_const,
      ofstream& out, TTree* tKee, TTree* tKemu, bool update, string plotsfile)
{

   //***************Get the PDFs from the workspace

   TFile fw(workspacename.c_str());   
   RooWorkspace* workspace = (RooWorkspace*)fw.Get("workspace");
   RooRealVar *B_plus_M = workspace->var("B_plus_M");
   RooRealVar *misPT = workspace->var("misPT");
   RooRealVar *T = workspace->var("T");
   RooRealVar *n = workspace->var("n");
   RooRealVar *expoConstCombHighPPerp = workspace->var("expoConstCombHighPPerp");
   RooRealVar *expoConstCombLowPPerp = workspace->var("expoConstCombLowPPerp");
   RooRealVar *expoConstKemu = workspace->var("expoConstKemu");
   RooRealVar *fractionalErrorJpsiLeak = workspace->var("fractionalErrorJpsiLeak");

   RooHistPdf *histPdfSignalZeroGammaHighPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGammaHighPPerp");
   RooHistPdf *histPdfSignalOneGammaHighPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalOneGammaHighPPerp");
   RooHistPdf *histPdfSignalTwoGammaHighPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGammaHighPPerp");
   RooHistPdf *histPdfPartRecoHighPPerp = (RooHistPdf *) workspace->pdf("histPdfPartRecoHighPPerp");
   RooHistPdf *histPdfJpsiLeakHighPPerp(0);
   if(nGenJpsiLeak>0) histPdfJpsiLeakHighPPerp = (RooHistPdf *) workspace->pdf("histPdfJpsiLeakHighPPerp");

   RooHistPdf *histPdfSignalZeroGammaLowPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGammaLowPPerp");
   RooHistPdf *histPdfSignalOneGammaLowPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalOneGammaLowPPerp");
   RooHistPdf *histPdfSignalTwoGammaLowPPerp = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGammaLowPPerp");
   RooHistPdf *histPdfPartRecoLowPPerp = (RooHistPdf *) workspace->pdf("histPdfPartRecoLowPPerp");
   RooHistPdf *histPdfJpsiLeakLowPPerp(0);
   if(nGenJpsiLeak>0) histPdfJpsiLeakLowPPerp = (RooHistPdf *) workspace->pdf("histPdfJpsiLeakLowPPerp");


   RooExponential* pdfCombLowPPerp =  new RooExponential("pdfCombLowPPerp", "pdfCombLowPPerp", *B_plus_M, *expoConstCombLowPPerp);
   RooExponential* pdfCombHighPPerp =  new RooExponential("pdfCombHighPPerp", "pdfCombHighPPerp", *B_plus_M, *expoConstCombHighPPerp);
   RooPTMVis* pdfKemu = new RooPTMVis("pdfKemu", "pdfKemu", *misPT, *B_plus_M, *T, *n, *expoConstKemu);

   //*************** get the true variables for randomisation

   RooRealVar *trueT = workspace->var("trueT");
   RooRealVar *trueN = workspace->var("trueN");
   RooRealVar *trueExpoConstCombHighPPerp = workspace->var("trueExpoConstCombHighPPerp");
   RooRealVar *trueExpoConstCombLowPPerp = workspace->var("trueExpoConstCombLowPPerp");
   RooRealVar *trueExpoConstKemu = workspace->var("trueExpoConstKemu");

   //************** Get the generated datasets to fit

   RooWorkspace* workspaceGen = (RooWorkspace*)fw.Get("workspaceGen");

   RooDataSet* dataGenSignalZeroGammaHighPPerp = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGammaHighPPerp");
   RooDataSet* dataGenSignalOneGammaHighPPerp = (RooDataSet*)workspaceGen->data("dataGenSignalOneGammaHighPPerp");
   RooDataSet* dataGenSignalTwoGammaHighPPerp = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGammaHighPPerp");
   RooDataSet* dataGenPartRecoHighPPerp = (RooDataSet*)workspaceGen->data("dataGenPartRecoHighPPerp");
   RooDataSet* dataGenCombHighPPerp = (RooDataSet*)workspaceGen->data("dataGenCombHighPPerp");
   RooDataSet* dataGenJpsiLeakHighPPerp(0);
   if(nGenJpsiLeak>0) dataGenJpsiLeakHighPPerp = (RooDataSet*)workspaceGen->data("dataGenJpsiLeakHighPPerp");

   RooDataSet* dataGenSignalZeroGammaLowPPerp = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGammaLowPPerp");
   RooDataSet* dataGenSignalOneGammaLowPPerp = (RooDataSet*)workspaceGen->data("dataGenSignalOneGammaLowPPerp");
   RooDataSet* dataGenSignalTwoGammaLowPPerp = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGammaLowPPerp");
   RooDataSet* dataGenPartRecoLowPPerp = (RooDataSet*)workspaceGen->data("dataGenPartRecoLowPPerp");
   RooDataSet* dataGenCombLowPPerp = (RooDataSet*)workspaceGen->data("dataGenCombLowPPerp");
   RooDataSet* dataGenJpsiLeakLowPPerp(0);
   if(nGenJpsiLeak>0) dataGenJpsiLeakLowPPerp = (RooDataSet*)workspaceGen->data("dataGenJpsiLeakLowPPerp");

   RooDataSet* dataGenKemu = (RooDataSet*)workspaceGen->data("dataGenKemu");

   if(wantplot)
   {
      //**************Must get the datasets

      RooDataSet* dataSetSignalZeroGammaHighPPerp = (RooDataSet*)workspace->data("dataSetSignalZeroGammaHighPPerp");
      RooDataSet* dataSetSignalOneGammaHighPPerp = (RooDataSet*)workspace->data("dataSetSignalOneGammaHighPPerp");
      RooDataSet* dataSetSignalTwoGammaHighPPerp = (RooDataSet*)workspace->data("dataSetSignalTwoGammaHighPPerp");
      RooDataSet* dataSetPartRecoHighPPerp = (RooDataSet*)workspace->data("dataSetPartRecoHighPPerp");
      RooDataSet* dataSetCombHighPPerp = (RooDataSet*)workspace->data("dataSetCombHighPPerp");
      RooDataSet* dataSetJpsiLeakHighPPerp = (RooDataSet*)workspace->data("dataSetJpsiLeakHighPPerp");

      RooDataSet* dataSetSignalZeroGammaLowPPerp = (RooDataSet*)workspace->data("dataSetSignalZeroGammaLowPPerp");
      RooDataSet* dataSetSignalOneGammaLowPPerp = (RooDataSet*)workspace->data("dataSetSignalOneGammaLowPPerp");
      RooDataSet* dataSetSignalTwoGammaLowPPerp = (RooDataSet*)workspace->data("dataSetSignalTwoGammaLowPPerp");
      RooDataSet* dataSetPartRecoLowPPerp = (RooDataSet*)workspace->data("dataSetPartRecoLowPPerp");
      RooDataSet* dataSetCombLowPPerp = (RooDataSet*)workspace->data("dataSetCombLowPPerp");
      RooDataSet* dataSetJpsiLeakLowPPerp = (RooDataSet*)workspace->data("dataSetJpsiLeakLowPPerp");

      RooDataSet* dataSetKemu = (RooDataSet*)workspace->data("dataSetKemu");

      //**************Plot all the different components

      cout<<"fit2D:   "<<fit2D<<" "<<B_plus_M<<endl;

      PlotShape(*dataSetSignalZeroGammaHighPPerp, *dataGenSignalZeroGammaHighPPerp, *histPdfSignalZeroGammaHighPPerp, plotsfile, "cSignalZeroGammaHighPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalOneGammaHighPPerp, *dataGenSignalOneGammaHighPPerp, *histPdfSignalOneGammaHighPPerp, plotsfile, "cSignalOneGammaHighPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalTwoGammaHighPPerp, *dataGenSignalTwoGammaHighPPerp, *histPdfSignalTwoGammaHighPPerp, plotsfile, "cSignalTwoGammaHighPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetPartRecoHighPPerp, *dataGenPartRecoHighPPerp, *histPdfPartRecoHighPPerp, plotsfile, "cPartRecoHighPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetCombHighPPerp, *dataGenCombHighPPerp, *pdfCombHighPPerp, plotsfile, "cCombHighPPerp", *B_plus_M, *misPT);
      if(nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeakHighPPerp, *dataGenJpsiLeakHighPPerp, *histPdfJpsiLeakHighPPerp, plotsfile, "cJpsiLeakHighPPerp", *B_plus_M, *misPT);

      PlotShape(*dataSetSignalZeroGammaLowPPerp, *dataGenSignalZeroGammaLowPPerp, *histPdfSignalZeroGammaLowPPerp, plotsfile, "cSignalZeroGammaLowPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalOneGammaLowPPerp, *dataGenSignalOneGammaLowPPerp, *histPdfSignalOneGammaLowPPerp, plotsfile, "cSignalOneGammaLowPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetSignalTwoGammaLowPPerp, *dataGenSignalTwoGammaLowPPerp, *histPdfSignalTwoGammaLowPPerp, plotsfile, "cSignalTwoGammaLowPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetPartRecoLowPPerp, *dataGenPartRecoLowPPerp, *histPdfPartRecoLowPPerp, plotsfile, "cPartRecoLowPPerp", *B_plus_M, *misPT);
      PlotShape(*dataSetCombLowPPerp, *dataGenCombLowPPerp, *pdfCombLowPPerp, plotsfile, "cCombLowPPerp", *B_plus_M, *misPT);
      if(nGenJpsiLeak>1) PlotShape(*dataSetJpsiLeakLowPPerp, *dataGenJpsiLeakLowPPerp, *histPdfJpsiLeakLowPPerp, plotsfile, "cJpsiLeakLowPPerp", *B_plus_M, *misPT);

      cout<<dataSetKemu<<" "<<dataGenKemu<<" "<<pdfKemu<<" "<<B_plus_M<<" "<<misPT<<endl;
      fit2D = true;
      PlotShape(*dataSetKemu, *dataGenKemu, *pdfKemu, plotsfile, "cKemu", *B_plus_M, *misPT);
      fit2D = false;
   }


   //***************Merge datasets

   RooDataSet* dataGenTotHighPPerp(dataGenPartRecoHighPPerp);
   dataGenTotHighPPerp->append(*dataGenSignalZeroGammaHighPPerp);
   dataGenTotHighPPerp->append(*dataGenSignalOneGammaHighPPerp);
   dataGenTotHighPPerp->append(*dataGenSignalTwoGammaHighPPerp);
   dataGenTotHighPPerp->append(*dataGenCombHighPPerp);
   if(nGenJpsiLeak>0) dataGenTotHighPPerp->append(*dataGenJpsiLeakHighPPerp);

   RooDataSet* dataGenTotLowPPerp(dataGenPartRecoLowPPerp);
   dataGenTotLowPPerp->append(*dataGenSignalZeroGammaLowPPerp);
   dataGenTotLowPPerp->append(*dataGenSignalOneGammaLowPPerp);
   dataGenTotLowPPerp->append(*dataGenSignalTwoGammaLowPPerp);
   dataGenTotLowPPerp->append(*dataGenCombLowPPerp);
   if(nGenJpsiLeak>0) dataGenTotLowPPerp->append(*dataGenJpsiLeakLowPPerp);

   //**************Create index category and join samples

   RooCategory category("category", "category");
   category.defineType("KeeHighPPerp");
   category.defineType("KeeLowPPerp");

   RooDataSet dataGenSimultaneous("dataGenSimultaneous", "dataGenSimultaneous", RooArgSet(*B_plus_M), Index(category), Import("KeeHighPPerp", *dataGenTotHighPPerp), Import("KeeLowPPerp", *dataGenTotLowPPerp));


   //**************Prepare fitting functions

   RooRealVar nSignal("nSignal", "#signal events", 1.*nGenSignal, nGenSignal-7*sqrt(nGenSignal), nGenSignal+7*sqrt(nGenSignal));
   RooRealVar nPartReco("nPartReco", "#nPartReco", 1.*nGenPartReco, nGenPartReco-7*sqrt(nGenPartReco), nGenPartReco+7*sqrt(nGenPartReco));
   RooRealVar nComb("nComb", "#nComb", 1.*nGenComb, nGenComb-7*sqrt(nGenComb), nGenComb+7*sqrt(nGenComb));
   RooRealVar nJpsiLeak("nJpsiLeak", "#nJpsiLeak", 1.*nGenJpsiLeak, nGenJpsiLeak-7*sqrt(nGenJpsiLeak), nGenJpsiLeak+7*sqrt(nGenJpsiLeak));
   
   
   //**************Construct PDFs for Signal


   RooRealVar fracZeroGammaLowPPerp("fracZeroGammaLowPPerp", "fracZeroGammaLowPPerp", 0.5, 0, 1);
   RooRealVar fracOneGammaLowPPerp("fracOneGammaLowPPerp", "fracOneGammaLowPPerp", 0.5, 0, 1);

   RooAddPdf pdfSignalLowPPerp("pdfSignalLowPPerp", "pdfSignalLowPPerp", RooArgList(*histPdfSignalZeroGammaLowPPerp, *histPdfSignalOneGammaLowPPerp, *histPdfSignalTwoGammaLowPPerp),
      RooArgList(fracZeroGammaLowPPerp, fracOneGammaLowPPerp));
   
  
   RooRealVar fracZeroGammaLowPPerpConstraintMean("fracZeroGammaLowPPerpConstraintMean", "fracZeroGammaLowPPerpConstraintMean", nGenFracZeroGammaLowPPerp);
   RooRealVar fracZeroGammaLowPPerpConstraintSigma("fracZeroGammaLowPPerpConstraintSigma", "fracZeroGammaLowPPerpConstraintSigma", 0.012*nGenFracZeroGammaLowPPerp);
   RooGaussian fracZeroGammaLowPPerpConstraint("fracZeroGammaLowPPerpConstraint", "fracZeroGammaLowPPerpConstraint", fracZeroGammaLowPPerp, fracZeroGammaLowPPerpConstraintMean, fracZeroGammaLowPPerpConstraintSigma); 

   RooRealVar fracOneGammaLowPPerpConstraintMean("fracOneGammaLowPPerpConstraintMean", "fracOneGammaLowPPerpConstraintMean", nGenFracOneGammaLowPPerp);
   RooRealVar fracOneGammaLowPPerpConstraintSigma("fracOneGammaLowPPerpConstraintSigma", "fracOneGammaLowPPerpConstraintSigma", 0.012*nGenFracOneGammaLowPPerp);
   RooGaussian fracOneGammaLowPPerpConstraint("fracOneGammaLowPPerpConstraint", "fracOneGammaLowPPerpConstraint", fracOneGammaLowPPerp, fracOneGammaLowPPerpConstraintMean, fracOneGammaLowPPerpConstraintSigma); 
   
  

   RooRealVar fracZeroGammaHighPPerp("fracZeroGammaHighPPerp", "fracZeroGammaHighPPerp", 0.5, 0, 1);
   RooRealVar fracOneGammaHighPPerp("fracOneGammaHighPPerp", "fracOneGammaHighPPerp", 0.5, 0, 1);

   RooAddPdf pdfSignalHighPPerp("pdfSignalHighPPerp", "pdfSignalHighPPerp", RooArgList(*histPdfSignalZeroGammaHighPPerp, *histPdfSignalOneGammaHighPPerp, *histPdfSignalTwoGammaHighPPerp),
      RooArgList(fracZeroGammaLowPPerp, fracOneGammaLowPPerp));
   
  
   RooRealVar fracZeroGammaHighPPerpConstraintMean("fracZeroGammaHighPPerpConstraintMean", "fracZeroGammaHighPPerpConstraintMean", nGenFracZeroGammaHighPPerp);
   RooRealVar fracZeroGammaHighPPerpConstraintSigma("fracZeroGammaHighPPerpConstraintSigma", "fracZeroGammaHighPPerpConstraintSigma", 0.04*nGenFracZeroGammaHighPPerp);
   RooGaussian fracZeroGammaHighPPerpConstraint("fracZeroGammaHighPPerpConstraint", "fracZeroGammaHighPPerpConstraint", fracZeroGammaHighPPerp, fracZeroGammaHighPPerpConstraintMean, fracZeroGammaHighPPerpConstraintSigma); 

   RooRealVar fracOneGammaHighPPerpConstraintMean("fracOneGammaHighPPerpConstraintMean", "fracOneGammaHighPPerpConstraintMean", nGenFracOneGammaHighPPerp);
   RooRealVar fracOneGammaHighPPerpConstraintSigma("fracOneGammaHighPPerpConstraintSigma", "fracOneGammaHighPPerpConstraintSigma", 0.04*nGenFracOneGammaHighPPerp);
   RooGaussian fracOneGammaHighPPerpConstraint("fracOneGammaHighPPerpConstraint", "fracOneGammaHighPPerpConstraint", fracOneGammaHighPPerp, fracOneGammaHighPPerpConstraintMean, fracOneGammaHighPPerpConstraintSigma); 
  
  
  
   //**************Prepare the yields for everything and constraint

   RooRealVar fracSigLowOverTot("fracSigLowOverTot", "fracSigLowOverTot", nGenFracSigLowOverTot);
   fracSigLowOverTot.setConstant(true);
   RooFormulaVar nSignalLowPPerp("nSignalLowPPerp", "@0*@1", RooArgList(fracSigLowOverTot, nSignal) );
   RooFormulaVar nSignalHighPPerp("nSignalHighPPerp", "(1-@0)*@1", RooArgList(fracSigLowOverTot, nSignal) );
  
   RooRealVar fracPartRecoLowOverTot("fracPartRecoLowOverTot", "fracPartRecoLowOverTot", nGenFracPartRecoLowOverTot);
   fracPartRecoLowOverTot.setConstant(true);
   RooFormulaVar nPartRecoLowPPerp("nPartRecoLowPPerp", "@0*@1", RooArgList(fracPartRecoLowOverTot, nPartReco) );
   RooFormulaVar nPartRecoHighPPerp("nPartRecoHighPPerp", "(1-@0)*@1", RooArgList(fracSigLowOverTot, nPartReco) );
   
   RooRealVar fracJpsiLeakLowOverTot("fracJpsiLeakLowOverTot", "fracJpsiLeakLowOverTot", nGenFracJpsiLeakLowOverTot);
   fracJpsiLeakLowOverTot.setConstant(true);
   RooFormulaVar nJpsiLeakLowPPerp("nJpsiLeakLowPPerp", "@0*@1", RooArgList(fracJpsiLeakLowOverTot, nJpsiLeak));
   RooFormulaVar nJpsiLeakHighPPerp("nJpsiLeakHighPPerp", "(1-@0)*@1", RooArgList(fracJpsiLeakLowOverTot, nJpsiLeak));

   RooRealVar nJpsiLeakConstraintMean("nJpsiLeakConstraintMean", "nJpsiLeakConstraintMean", nGenJpsiLeak);
   RooRealVar nJpsiLeakConstraintSigma("nJpsiLeakConstraintSigma", "nJpsiLeakConstraintSigma", fractionalErrorJpsiLeak->getVal()*nGenJpsiLeak);  
   RooGaussian nJpsiLeakConstraint("nJpsiLeakConstraint", "nJpsiLeakConstraint", nJpsiLeak, nJpsiLeakConstraintMean, nJpsiLeakConstraintSigma);
   
   RooRealVar fracCombLowOverTot("fracCombLowOverTot", "fracCombLowOverTot", 0,1);
   RooFormulaVar nCombLowPPerp("nCombLowPPerp", "@0*@1", RooArgList(fracCombLowOverTot, nComb));
   RooFormulaVar nCombHighPPerp("nCombHighPPerp", "(1-@0)*@1", RooArgList(fracCombLowOverTot, nComb));


   //************** Construt the total PDF for Kee samples

   RooArgList pdfListLowPPerp(pdfSignalLowPPerp, *histPdfPartRecoLowPPerp, *pdfCombLowPPerp);
   RooArgList yieldListLowPPerp(nSignalLowPPerp, nPartRecoLowPPerp, nCombLowPPerp);

   if(nGenJpsiLeak>0)
   {
      pdfListLowPPerp.add(*histPdfJpsiLeakLowPPerp);
      yieldListLowPPerp.add(nJpsiLeakLowPPerp); 
   }
   RooAddPdf totPdfLowPPerp("totPdfLowPPerp", "totPdfLowPPerp", pdfListLowPPerp, yieldListLowPPerp); 

   RooArgList pdfListHighPPerp(pdfSignalHighPPerp, *histPdfPartRecoHighPPerp, *pdfCombHighPPerp);
   RooArgList yieldListHighPPerp(nSignalHighPPerp, nPartRecoHighPPerp, nCombHighPPerp);

   if(nGenJpsiLeak>0)
   {
      pdfListHighPPerp.add(*histPdfJpsiLeakHighPPerp);
      yieldListHighPPerp.add(nJpsiLeakHighPPerp); 
   }
   RooAddPdf totPdfHighPPerp("totPdfHighPPerp", "totPdfHighPPerp", pdfListHighPPerp, yieldListHighPPerp); 

   RooSimultaneous simPdf("simPdf", "simPdf", category);
   simPdf.addPdf(totPdfHighPPerp, "KeeHighPPerp");
   simPdf.addPdf(totPdfLowPPerp, "KeeLowPPerp");


   //**************** Fit Kemu sample

   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-8) ;



   RooAbsReal* nllKemu = pdfKemu->createNLL(*dataGenKemu);
   RooMinuit minuitKemu(*nllKemu);
   minuitKemu.setStrategy(2);

   int migradRes(1);
   int hesseRes(4);

   vector<int> migradResVec;
   vector<int> hesseResVec;

   double edm(10);
   int nrefit(0);

   RooFitResult* fitRes(0);
   vector<RooFitResult*> fitResVec;

   bool hasConverged(false);

   for(int i(0); (i<10) && !hasConverged ; ++i)
   {
      initiateParamsKemu(*expoConstKemu, *trueExpoConstKemu, *T, *trueT, *n, *trueN);

      migradRes = minuitKemu.migrad();
      hesseRes = minuitKemu.hesse();

      fitRes = minuitKemu.save();
      edm = fitRes->edm();

      fitResVec.push_back(fitRes); 
      migradResVec.push_back(migradRes);
      hesseResVec.push_back(hesseRes);

      if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;

      ++nrefit;
      cout<<"Fitting nbr "<<i<<", Kemu, done. Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   if(!hasConverged)
   {
      double minNll(1e20);
      int minIndex(-1);
      for(unsigned int i(0); i<fitResVec.size(); ++i)
      {
         if( fitResVec.at(i)->minNll() < minNll)
         {
            minIndex = i;
            minNll = fitResVec[i]->minNll();
         }
      }

      migradRes = migradResVec.at(minIndex);
      hesseRes = hesseResVec.at(minIndex);
      cout<<"Fit Kemu not converged, choose fit "<<minIndex<<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(tKemu, fitRes,  update, migradRes, hesseRes, hasConverged);


   //prepare integral and combinatorial constraint to the fit

   RooAbsReal* igLowPPerp = pdfKemu->createIntegral(*misPT,RooFit::NormSet(*misPT),RooFit::Range("lowPPerp"));
   RooAbsReal* igAllPPerp = pdfKemu->createIntegral(*misPT,RooFit::NormSet(*misPT),RooFit::Range("default"));
   RooFormulaVar fracCombFormulaLowOverTot("fracCombFormulaLowOverTot", "fracCombFormulaLowOverTot", "@0/@1" , RooArgList(*igLowPPerp, *igAllPPerp));

   RooRealVar fracCombLowOverTotFromFit("fracCombLowOverTotFromFit", "fracCombLowOverTotFromFit", fracCombFormulaLowOverTot.getVal(), 0,1);
   fracCombLowOverTotFromFit.setError(fracCombFormulaLowOverTot.getPropagatedError(*fitRes));
   RooRealVar fracCombLowOverTotConstraintMean("fracCombLowOverTotConstraintMean", "fracCombLowOverTotConstraintMean", fracCombFormulaLowOverTot.getVal());
   RooRealVar fracCombLowOverTotConstraintSigma("fracCombLowOverTotConstraintSigma", "fracCombLowOverTotConstraintSigma", fracCombFormulaLowOverTot.getPropagatedError(*fitRes));
   RooGaussian fracCombLowOverTotConstraint("fracCombLowOverTotConstraint", "fracCombLowOverTotConstraint", fracCombLowOverTot, fracCombLowOverTotConstraintMean, fracCombLowOverTotConstraintSigma);

   

   //delete results of Kemu fit

   for(unsigned int i(0); i<fitResVec.size(); ++i) delete fitResVec.at(i);

   //fit Kee sample

   RooArgSet constraints(fracZeroGammaLowPPerpConstraint, fracOneGammaLowPPerpConstraint,  fracZeroGammaHighPPerpConstraint, fracOneGammaHighPPerpConstraint, nJpsiLeakConstraint,  fracCombLowOverTotConstraint);

   RooAbsReal* nllKee = simPdf.createNLL(dataGenSimultaneous,  Extended(), ExternalConstraints(constraints));;
   RooMinuit minuitKee(*nllKee);
   minuitKee.setStrategy(2);

   migradRes = 1;
   hesseRes = 4;

   migradResVec.clear();
   hesseResVec.clear();

   edm = 10;
   nrefit = 0;

   fitRes = 0;
   fitResVec.clear();

   hasConverged = false;

   for(int i(0); (i<10) && !hasConverged ; ++i)
   {
      initiateParamsKee(nSignal, nPartReco, nComb, nJpsiLeak,
            fracZeroGammaLowPPerp, fracZeroGammaLowPPerpConstraintSigma, fracOneGammaLowPPerp, fracOneGammaLowPPerpConstraintSigma,
            fracZeroGammaHighPPerp, fracZeroGammaHighPPerpConstraintSigma, fracOneGammaHighPPerp, fracOneGammaHighPPerpConstraintSigma,
            *expoConstCombLowPPerp, *trueExpoConstCombLowPPerp, *expoConstCombHighPPerp, *trueExpoConstCombHighPPerp,
            fracCombLowOverTot, fracCombLowOverTotFromFit);

      migradRes = minuitKee.migrad();
      hesseRes = minuitKee.hesse();

      fitRes = minuitKee.save();
      edm = fitRes->edm();

      fitResVec.push_back(fitRes); 
      migradResVec.push_back(migradRes);
      hesseResVec.push_back(hesseRes);

      if( migradRes == 0 && hesseRes == 0 && edm < 1e-4 ) hasConverged = true;

      ++nrefit;
      cout<<"Fitting nbr "<<i<<", Kee, done. Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   if(!hasConverged)
   {
      double minNll(1e20);
      int minIndex(-1);
      for(unsigned int i(0); i<fitResVec.size(); ++i)
      {
         if( fitResVec.at(i)->minNll() < minNll)
         {
            minIndex = i;
            minNll = fitResVec[i]->minNll();
         }
      }

      migradRes = migradResVec.at(minIndex);
      hesseRes = hesseResVec.at(minIndex);
      cout<<"Fit Kee not converged, choose fit "<<minIndex<<". Hesse: "<<hesseRes<<" migrad: "<<migradRes<<" edm: "<<edm<<" minNll: "<<fitRes->minNll()<<endl;
   }


   fillTreeResult(tKee, fitRes,  update, migradRes, hesseRes, hasConverged);


   //*************** output fit status


   int w(12);
   out<<setw(w)<<migradRes<<setw(w)<<hesseRes<<setw(w)<<edm<<setw(w)<<nrefit<<endl;

   if(wantplot)
   {
       plot_fit_result_Kemu(plotsfile, *pdfKemu, *dataGenKemu);
       plot_fit_result_Kee(*B_plus_M, category, plotsfile, simPdf, dataGenSimultaneous);
   }

   fw.Close();
   //delete and return
   delete nllKee;
   delete nllKemu;
   delete workspace;
   delete workspaceGen;
   delete pdfCombLowPPerp;
   delete pdfCombHighPPerp;
   delete pdfKemu;

}


void FitterUtils1DSimultaneous::plot_fit_result_Kee(RooRealVar B_plus_M, RooCategory& category, string plotsfile, RooAbsPdf &totPdf, RooDataSet& dataGenTot)
{

   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit


   RooPlot* frameLowPPerp = B_plus_M.frame();
   dataGenTot.plotOn(frameLowPPerp, Cut("category==category::KeeLowPPerp"));
   totPdf.plotOn(frameLowPPerp, Components("histPdfPartRecoLowPPerp"), Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot),  LineColor(kBlue));
   totPdf.plotOn(frameLowPPerp, Components("histPdfSignalZeroGammaLowPPerp"), Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot), LineColor(kGreen));
   totPdf.plotOn(frameLowPPerp, Components("histPdfSignalOneGammaLowPPerp"), Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot), LineColor(kMagenta));
   totPdf.plotOn(frameLowPPerp, Components("histPdfSignalTwoGammaLowPPerp"), Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot), LineColor(kOrange));
   totPdf.plotOn(frameLowPPerp, Components("histPdfJpsiLeakLowPPerp"), Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot), LineColor(14));
   totPdf.plotOn(frameLowPPerp, Components("pdfCombLowPPerp"), Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot), LineColor(kBlack));
   totPdf.plotOn(frameLowPPerp, Slice(category, "KeeLowPPerp"), ProjWData(category, dataGenTot), LineColor(kRed));

   RooPlot* frameHighPPerp = B_plus_M.frame();
   dataGenTot.plotOn(frameHighPPerp, Cut("category==category::KeeHighPPerp"));
   totPdf.plotOn(frameHighPPerp, Components("histPdfPartRecoHighPPerp"), Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot),  LineColor(kBlue));
   totPdf.plotOn(frameHighPPerp, Components("histPdfSignalZeroGammaHighPPerp"), Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot), LineColor(kGreen));
   totPdf.plotOn(frameHighPPerp, Components("histPdfSignalOneGammaHighPPerp"), Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot), LineColor(kMagenta));
   totPdf.plotOn(frameHighPPerp, Components("histPdfSignalTwoGammaHighPPerp"), Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot), LineColor(kOrange));
   totPdf.plotOn(frameHighPPerp, Components("histPdfJpsiLeakHighPPerp"), Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot), LineColor(14));
   totPdf.plotOn(frameHighPPerp, Components("pdfCombHighPPerp"), Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot), LineColor(kBlack));
   totPdf.plotOn(frameHighPPerp, Slice(category, "KeeHighPPerp"), ProjWData(category, dataGenTot), LineColor(kRed));



   TCanvas cFit("cFitKee", "cFitKee", 600, 800);
   cFit.Divide(1,2);
   cFit.cd(1); frameLowPPerp->Draw();
   cFit.cd(2); frameHighPPerp->Draw();

   cFit.Write();
   f2.Close();
}

void FitterUtils1DSimultaneous::plot_fit_result_Kemu(string plotsfile, RooAbsPdf &totPdf, RooDataSet& dataGenTot)
{

   //**************Prepare TFile to save the plots

   TFile f2(plotsfile.c_str(), "UPDATE");
   //**************Plot the results of the fit

   RooArgSet *var_set = totPdf.getObservables(dataGenTot);
   TIterator *iter = var_set->createIterator();
   RooRealVar *var;

   std::vector<RooPlot*> plots;
   RooPlot* frame;

   while((var = (RooRealVar*) iter->Next()))
   {

      frame = var->frame();
      dataGenTot.plotOn(frame);
      totPdf.plotOn(frame, LineColor(kRed));

      plots.push_back(frame);

   }

   if (!(plots.size())) return;

   TCanvas cFit("cFitKemu", "cFitKemu", 600, 800);
   plots[0]->Draw();

   cFit.Write();
   f2.Close();


}


void FitterUtils1DSimultaneous::display()
{
   cout<<"nGenFracZeroGammaLowPPerp "<<nGenFracZeroGammaLowPPerp<<endl;
   cout<<"nGenFracOneGammaLowPPerp "<<nGenFracOneGammaLowPPerp<<endl;
   cout<<"nGenFracZeroGammaHighPPerp "<<nGenFracZeroGammaHighPPerp<<endl;
   cout<<"nGenFracOneGammaHighPPerp "<<nGenFracOneGammaHighPPerp<<endl;
   cout<<"nGenSignalLowPPerp "<<nGenSignalLowPPerp<<endl;
   cout<<"nGenSignalZeroGammaLowPPerp "<<nGenSignalZeroGammaLowPPerp<<endl;
   cout<<"nGenSignalOneGammaLowPPerp "<<nGenSignalOneGammaLowPPerp<<endl;
   cout<<"nGenSignalTwoGammaLowPPerp "<<nGenSignalTwoGammaLowPPerp<<endl;
   cout<<"nGenPartRecoLowPPerp "<<nGenPartRecoLowPPerp<<endl;
   cout<<"nGenCombLowPPerp "<<nGenPartRecoLowPPerp<<endl;
   cout<<"nGenJpsiLeakLowPPerp "<<nGenJpsiLeakLowPPerp<<endl;
   cout<<"nGenSignalHighPPerp "<<nGenSignalHighPPerp<<endl;
   cout<<"nGenSignalZeroGammaHighPPerp "<<nGenSignalZeroGammaHighPPerp<<endl;
   cout<<"nGenSignalOneGammaHighPPerp "<<nGenSignalOneGammaHighPPerp<<endl;
   cout<<"nGenSignalTwoGammaHighPPerp "<<nGenSignalTwoGammaHighPPerp<<endl;
   cout<<"nGenPartRecoHighPPerp "<<nGenPartRecoHighPPerp<<endl;
   cout<<"nGenCombHighPPerp "<<nGenCombHighPPerp<<endl;
   cout<<"nGenJpsiLeakHighPPerp "<<nGenJpsiLeakHighPPerp<<endl;
   cout<<"nGenFracSigLowOverTot "<<nGenFracSigLowOverTot<<endl;
   cout<<"nGenFracSigZeroGammaLowOverTot "<<nGenFracSigZeroGammaLowOverTot<<endl;
   cout<<"nGenFracSigOneGammaLowOverTot "<<nGenFracSigOneGammaLowOverTot<<endl;
   cout<<"nGenFracSigTwoGammaLowOverTot "<<nGenFracSigTwoGammaLowOverTot<<endl;
   cout<<"nGenFracCombLowOverTot "<<nGenFracCombLowOverTot<<endl;
   cout<<"nGenFracPartRecoLowOverTot "<<nGenFracPartRecoLowOverTot<<endl;
   cout<<"nGenFracJpsiLeakLowOverTot "<<nGenFracJpsiLeakLowOverTot<<endl;
   cout<<"nGenSignalZeroGamma "<<nGenSignalZeroGamma<<endl;
   cout<<"nGenSignalOneGamma "<<nGenSignalOneGamma<<endl;
   cout<<"nGenSignalTwoGamma "<<nGenSignalTwoGamma<<endl;


   cout<<"FracCombRealVarLowOverTot "<<FracCombRealVarLowOverTot.getVal()<<" "<<FracCombRealVarLowOverTot.getError()<<endl;
   cout<<"FracCombFromFitLowOverTot "<<FracCombFromFitLowOverTot<<endl;
   cout<<"nCombFromFitLowPPerp "<<nCombFromFitLowPPerp<<endl;
   cout<<"nCombFromFitHighPerp "<<nCombFromFitHighPerp<<endl;

   cout<<"nGenSignal "<<nGenSignal<<endl;
   cout<<"nGenPartReco "<<nGenPartReco<<endl;
   cout<<"nGenComb "<<nGenComb<<endl;
   cout<<"nGenJpsiLeak "<<nGenJpsiLeak<<endl;
   cout<<"nGenFracZeroGamma "<<nGenFracZeroGamma<<endl;
   cout<<"nGenFracOneGamma "<<nGenFracOneGamma<<endl;


   cout<<"nGenKemu "<<nGenKemu<<endl;
   cout<<"PPerpCut "<<PPerpCut<<endl;

}




