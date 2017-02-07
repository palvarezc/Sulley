#define analysis_histfact_working_cxx
#include <iostream>
#include <stdio.h>

#include "TRandom3.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TIterator.h"
#include "TH3.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooParamHistFunc.h"
#include "RooHist.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"


using namespace std;
using namespace RooFit;
using namespace RooStats ;
using namespace HistFactory ;

TDatime *date = new TDatime();

void analysis_roofit_working()
{

  //Some stylistic things:
  TLatex *t=new TLatex();
  t->SetTextAlign(22);
  t->SetTextSize(0.06);
  t->SetTextFont(132);
  gROOT->ProcessLine("gStyle->SetLabelFont(132,\"xyz\");");
  gROOT->ProcessLine("gStyle->SetTitleFont(132,\"xyz\");");
  gROOT->ProcessLine("gStyle->SetTitleFont(132,\"t\");");
  gROOT->ProcessLine("gStyle->SetTitleSize(0.08,\"t\");");
  gROOT->ProcessLine("gStyle->SetTitleY(0.970);");
  char substr[128];
  //For ToyMC, so I can run multiple copies with different seeds without recompiling
  RooRandom::randomGenerator()->SetSeed(date->Get()%100000);
  cout << date->Get()%100000 << endl; 

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ; //avoid accidental unblinding!

  //Open the ROOT file with the histograms in and renormalise the fit components to 1.

  TFile q("toyMaker/histograms.root");
  TH2 *htemp;

  TString mchistos[4]={"signal_zero","signal_one","signal_two","partreco"};
  double mcN_sigzero,mcN_sigone,mcN_sigtwo , mcN_partreco;
  double *mcnorms[4]={&mcN_sigzero,&mcN_sigone,&mcN_sigtwo ,& mcN_partreco};
  for(int i =0; i < 4; i++){
    q.GetObject(mchistos[i],htemp);
    assert(htemp!=NULL);
    *(mcnorms[i])=1./htemp->Integral();
    cout << "mcN_"+mchistos[i]+" = " << 1./ *(mcnorms[i]) << endl;
  }

  // // Useful later to have the bin max and min for drawing
  // TH2 *JUNK;
  // q.GetObject("partreco",JUNK);
  // double q2_low=JUNK->GetZaxis()->GetXmin();
  // double q2_high=JUNK->GetZaxis()->GetXmax();
  // const int q2_bins = JUNK->GetZaxis()->GetNbins();
  // JUNK->SetDirectory(0);
  q.Close();

  // Many many flags for steering
  /* STEERING OPTIONS */ 
  const bool constrainDstst=true;
  const bool useMinos=true;
  const bool useMuShapeUncerts = false;
  const bool useTauShapeUncerts = false;
  const bool fixshapes = true;
  const bool dofit = true;
  const bool toyMC = false;
  const bool fitfirst = false;
  const bool slowplots = true;
  const bool BBon3d = false;
  const bool blind = false;
  const int  blindingSeed = 110689 ;
 

  double Nsig = 560.;
  double Npartreco = 550.;
  double fpartreco = 1;
  
// double fzero = 1.;
  // double fone = 0;
  // double ftwo = 0;
  
  double fzero = 0.368;
  double fone = 0.484;
  double ftwo = 1- fzero - fone;


  // TFile fin("toyMaker/histograms.root");
  // TH2D *signal_zero_hist = (TH2D*) fin.Get("signal_zero");
  // TH2D *signal_one_hist = (TH2D*) fin.Get("signal_one");
  // TH2D *signal_two_hist = (TH2D*) fin.Get("signal_two");
  // TH2D *partreco_hist = (TH2D*) fin.Get("partreco");
  // TH2D *data_hist = (TH2D*) fin.Get("dataToy");  


  std::string workspacename = "toy_result_notconv/workspaceMode7.root";
  std::string channelName = "B2Kee";
  TFile *fw = new TFile(workspacename.c_str());
  RooWorkspace* workspace = (RooWorkspace*)fw->Get("combined");
  ModelConfig *mc = (ModelConfig*) workspace->obj("ModelConfig");
  RooArgSet *obs = (RooArgSet*) mc->GetObservables();
  RooRealVar *x = (RooRealVar*) obs->find(("obs_x_"+channelName).c_str());
  RooRealVar *y = (RooRealVar*) obs->find(("obs_y_"+channelName).c_str());

  // x->setMin(signal_zero_hist->GetXaxis()->GetXmin());
  // x->setMax(signal_zero_hist->GetXaxis()->GetXmax());
  // y->setMin(signal_zero_hist->GetYaxis()->GetXmin());
  // y->setMax(signal_zero_hist->GetYaxis()->GetXmax());
  

  RooWorkspace* workspaceGen = (RooWorkspace*)fw->Get("workspaceGen");
  RooDataSet* dataGenSignalZeroGamma = (RooDataSet*)workspaceGen->data("dataGenSignalZeroGamma");
  RooDataSet* dataGenSignalOneGamma = (RooDataSet*)workspaceGen->data("dataGenSignalOneGamma");
  RooDataSet* dataGenSignalTwoGamma = (RooDataSet*)workspaceGen->data("dataGenSignalTwoGamma");
  RooDataSet* dataGenPartReco = (RooDataSet*)workspaceGen->data("dataGenPartReco");
   
  RooDataSet* dataGenTot_chan(dataGenPartReco);
  dataGenTot_chan->append(*dataGenSignalZeroGamma);
  dataGenTot_chan->append(*dataGenSignalOneGamma);
  dataGenTot_chan->append(*dataGenSignalTwoGamma);
  RooDataHist *dataGenTot_chan_bin = dataGenTot_chan->binnedClone();
  RooDataHist *data = new RooDataHist("dataGenTot","dataGenTot",RooArgSet(*x, *y),
                                            *dataGenTot_chan_bin);
  


  RooHistPdf *signal_zero_pdf_HF = (RooHistPdf *) workspace->pdf("histPdfSignalZeroGamma");
  RooHistPdf *signal_one_pdf_HF = (RooHistPdf *) workspace->pdf("histPdfSignalOneGamma");
  RooHistPdf *signal_two_pdf_HF = (RooHistPdf *) workspace->pdf("histPdfSignalTwoGamma");
  RooHistPdf *partreco_pdf_HF = (RooHistPdf *) workspace->pdf("histPdfPartReco");
     
  TH2D * signal_zero_hist = (TH2D*) signal_zero_pdf_HF->createHistogram("signal_zero",*x, RooFit::YVar(*y));
  TH2D * signal_one_hist = (TH2D*) signal_one_pdf_HF->createHistogram("signal_one",*x, RooFit::YVar(*y));
  TH2D * signal_two_hist = (TH2D*) signal_two_pdf_HF->createHistogram("signal_two",*x, RooFit::YVar(*y));
  TH2D * partreco_hist = (TH2D*) partreco_pdf_HF->createHistogram("partreco",*x, RooFit::YVar(*y));

  // RooRealVar x("x","x",signal_zero_hist->GetXaxis()->GetXmin(),signal_zero_hist->GetXaxis()->GetXmax());
  // RooRealVar y("y","y",signal_zero_hist->GetYaxis()->GetXmin(),signal_zero_hist->GetYaxis()->GetXmax());
  
  RooArgList varlist(*x,*y);
  RooArgSet varset(*x,*y);

  RooDataHist signal_zero_rhist("signal_zero_rhist","signal_zero_rhist",varlist,signal_zero_hist);
  RooDataHist signal_one_rhist("signal_one_rhist","signal_one_rhist",varlist,signal_one_hist);
  RooDataHist signal_two_rhist("signal_two_rhist","signal_two_rhist",varlist,signal_two_hist);
  RooDataHist partreco_rhist("partreco_rhist","partreco_rhist",varlist,partreco_hist);
  // RooDataHist data("data","data",varlist,data_hist);
  

  RooHistPdf signal_zero_pdf("signal_zero_pdf","signal_zero_pdf",varset,signal_zero_rhist);
  RooHistPdf signal_one_pdf("signal_one_pdf","signal_one_pdf",varset,signal_one_rhist);
  RooHistPdf signal_two_pdf("signal_two_pdf","signal_two_pdf",varset,signal_two_rhist);
  RooHistPdf partreco_pdf("partreco_pdf","partreco_pdf",varset,partreco_rhist);
  
  RooRealVar fzero_par("fzero","fzero",fzero);
  RooRealVar fone_par("fone","fone",fone);
  RooRealVar ftwo_par("ftwo","ftwo",ftwo);
  RooRealVar Nsig_par("Nsig","Nsig",Nsig,1e-6,1e6);
  RooRealVar Npartreco_par("Npartreco","Npartreco",Npartreco,1e-6,1e6);
  
  // RooAddPdf *model_signal = new RooAddPdf("model_signal","model_signal",RooArgList(*signal_zero_pdf,*signal_one_pdf,*signal_two_pdf),
  //                        RooArgList(fzero_par,fone_par));

  // RooAddPdf model_hf("MyModelWithCombinatorial","MyModelWithCombinatorial", RooArgList(*model_signal,*partreco_pdf),
  //                    RooArgList(Nsig_par,Npartreco_par));

  RooAddPdf *model_signal = new RooAddPdf("model_signal","model_signal",RooArgList(signal_zero_pdf,signal_one_pdf,signal_two_pdf),
                         RooArgList(fzero_par,fone_par));

  RooAddPdf model_hf("MyModelWithCombinatorial","MyModelWithCombinatorial", RooArgList(*model_signal,partreco_pdf),
                     RooArgList(Nsig_par,Npartreco_par));


  std::cout<<"Aqui"<<std::endl;
  

  


  //Create the workspace that is going to save the config
  RooWorkspace *w = new RooWorkspace("mywsp");  
  w->import(model_hf);
 
  RooFitResult *toyresult;
  RooAbsReal *nll_hf; 
  RooAbsReal *nll; 

  RooFitResult *result, *result2;

  cerr << "Saving PDF snapshot" << endl;
  RooArgSet allpars(fzero_par,fone_par,Nsig_par,Npartreco_par);
  w->saveSnapshot("TMCPARS",allpars,kTRUE);
  RooRealVar poierror("poierror","poierror",0.00001,0.010);
  TIterator *iter = allpars.createIterator();
  RooAbsArg *tempvar;
  RooArgSet *theVars = (RooArgSet*) allpars.Clone();
  theVars->add(poierror);
  

  nll_hf= model_hf.createNLL(*data,Offset(kTRUE));

  RooMinuit* minuit_hf = new RooMinuit(*nll_hf) ;
  RooArgSet *temp = new RooArgSet();
  nll_hf->getParameters(temp)->Print("V");
  cout << "**********************************************************************" << endl;
  minuit_hf->setErrorLevel(0.5);
  //#ifndef UNBLIND
  minuit_hf->setPrintLevel(1);
  //#endif
  std::cout << "Minimizing the Minuit (Migrad)" << std::endl;
  w->saveSnapshot("TMCPARS",allpars,kTRUE);
  minuit_hf->setStrategy(2);
  minuit_hf->fit("smh");
  RooFitResult *tempResult=minuit_hf->save("TempResult","TempResult");
  
  cout << tempResult->edm() << endl;
  if (useMinos) minuit_hf->minos(RooArgSet(Nsig_par));
  result = minuit_hf->save("Result","Result");
  
  
  
  if(result != NULL)
  {
    printf("Fit ran with status %d\n",result->status());
    
    printf("Stat error on Nsig is %f\n",Nsig_par.getError());
    
    printf("EDM at end was %f\n",result->edm());
      
    result->floatParsInit().Print();
      
    cout << "CURRENT PARAMETERS:" << endl;
    //TIterator *paramiter = mc->GetNuisanceParameters()->createIterator();
    TIterator *paramiter = result->floatParsFinal().createIterator();
    RooRealVar *__temp= (RooRealVar *)paramiter->Next();
    int final_par_counter=0;
    while (__temp!=NULL)
    {
      if(!__temp->isConstant())
      {
        //          if(!(TString(__temp->GetName()).EqualTo(poi->GetName())))
        //          {
        cout << final_par_counter << ": "
             << __temp->GetName() << "\t\t\t = "
             << ((RooRealVar*)result->floatParsFinal().find(__temp->GetName()))->getVal()
             << " +/- "
             << ((RooRealVar*)result->floatParsFinal().find(__temp->GetName()))->getError() << endl;
        //          }
      }
      final_par_counter++;
      __temp=(RooRealVar *)paramiter->Next();
    }
    
    
    result->correlationMatrix().Print();
  }
  
  //Now for the plotting
  RooPlot *m_frame = x->frame(Title("M(B)"));
  RooPlot *mcor_frame = y->frame(Title("M_{corr}(B)"));
  
  const int nframes = 2;
  RooPlot *drawframes[nframes] = {m_frame, mcor_frame};
  
  const int ncomps = 10;
  int colors[ncomps]={kRed,kBlue+1,kViolet,kViolet+1,kViolet+2,kGreen,kGreen+1,kOrange+1,kOrange+2,kOrange+3};
  const int ncomps2 = 8;
  TString names[ncomps2+1] = {"Data","Total Fit"
                              ,"B #rightarrow K e e"
                              ,"Part. Reco."
  };
      
      
  RooHist* mresid;// = mm2_frame->pullHist() ;
  RooHist* mcorresid;// = El_frame->pullHist() ;
  
  RooHist *resids[nframes];

  for (int i = 0; i < nframes; i++){
    data->plotOn(drawframes[i],Name("data"),DataError(RooAbsData::Poisson),MarkerSize(0.6),DrawOption("ZP"));
    model_hf.plotOn(drawframes[i],Name("model_hf"), DrawOption("F"),FillColor(kRed),LineWidth(0));
    resids[i]=drawframes[i]->pullHist();
    model_hf.plotOn(drawframes[i],Name("signal"), DrawOption("F"),FillColor(kBlue),Components("*Zero*,*One*,*Two*"),LineWidth(0));
    model_hf.plotOn(drawframes[i],Name("sigzero_sample"), Components("*Zero*"),LineWidth(1),LineColor(kYellow));
    model_hf.plotOn(drawframes[i],Name("sigone_sample"), LineColor(kViolet),Components("*One*"),LineWidth(1));
    model_hf.plotOn(drawframes[i],Name("sigtwo_sample"), LineColor(kGreen),Components("*Two*"),LineWidth(1));
    model_hf.plotOn(drawframes[i],Name("partreco_sample"), DrawOption("F"),FillColor(kMagenta),Components("*PartReco*"),LineWidth(0));
    
    data->plotOn(drawframes[i],Name("data"),DataError(RooAbsData::Poisson),MarkerSize(0.8),DrawOption("ZP"));
  }
  
  
  mresid=resids[0];
  mcorresid=resids[1];
  
  TCanvas *c1 = new TCanvas("c1","c1",1100,600);
  c1->SetTickx();
  c1->SetTicky();
  TPad * pad1 = new TPad("pad1","pad1",0,0.3,0.33,1.0);
  pad1->SetTickx();
  pad1->SetTicky();
  pad1->SetRightMargin(0.11);
  pad1->SetLeftMargin(0.17);
  pad1->SetTopMargin(0.08);
  pad1->SetBottomMargin(0.13);
  pad1->Draw();
  TPad * pad2 = new TPad("pad2","pad2",0.33,0.3,0.66,1.0);
  pad2->SetTickx();
  pad2->SetTicky();
  pad2->SetRightMargin(0.11);
  pad2->SetLeftMargin(0.17);
  pad2->SetTopMargin(0.08);
  pad2->SetBottomMargin(0.12);
  pad2->Draw();
  TPad * pad3 = new TPad("pad3","pad3",0.66,0.3,0.99,1.0);
  pad3->SetTickx();
  pad3->SetTicky();
  pad3->SetRightMargin(0.11);
  pad3->SetLeftMargin(0.17);
  pad3->SetTopMargin(0.08);
  pad3->SetBottomMargin(0.13);
  pad3->Draw();
  TPad * pad4 = new TPad("pad4","pad4",0,0.,0.33,0.3);
  pad4->SetTickx();
  pad4->SetTicky();
  pad4->SetRightMargin(0.11);
  pad4->SetLeftMargin(0.17);
  pad4->SetTopMargin(0.05);
  pad4->SetBottomMargin(0.13);
  pad4->Draw();
  TPad * pad5 = new TPad("pad5","pad5",0.33,0.,0.66,0.3);
  pad5->SetTickx();
  pad5->SetTicky();
  pad5->SetRightMargin(0.11);
  pad5->SetLeftMargin(0.17);
  pad5->SetTopMargin(0.05);
  pad5->SetBottomMargin(0.13);
  pad5->Draw();
  TPad * pad6 = new TPad("pad6","pad6",0.66,0.,0.99,0.3);
  pad6->SetTickx();
  pad6->SetTicky();
  pad6->SetRightMargin(0.11);
  pad6->SetLeftMargin(0.17);
  pad6->SetTopMargin(0.05);
  pad6->SetBottomMargin(0.13);
  pad6->Draw();

  TLegend * leg = new TLegend(0.57,0.6,0.87,0.91);
  leg->SetLineColor(0);
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);
  leg->AddEntry(m_frame->findObject("data"),"data","LP");
  leg->AddEntry(m_frame->findObject("model"),"ALL","F");
  leg->AddEntry(m_frame->findObject("signal"),"B #rightarrow K e e","F");
  leg->AddEntry(m_frame->findObject("sigzero_sample"),"Zero Brem","l");
  leg->AddEntry(m_frame->findObject("sigone_sample"),"One Brem","l");
  leg->AddEntry(m_frame->findObject("sigtwo_sample"),"Two Brem","l");
  leg->AddEntry(m_frame->findObject("partreco_sample"),"PartReco","F");

  pad1->cd();
  m_frame->SetTitle("");
  m_frame->GetXaxis()->SetLabelSize(0.06);
  m_frame->GetXaxis()->SetTitleSize(0.06);
  m_frame->GetYaxis()->SetLabelSize(0.06);
  m_frame->GetYaxis()->SetTitleSize(0.06);
  m_frame->GetYaxis()->SetTitleOffset(1.3);
  m_frame->GetXaxis()->SetTitleOffset(0.9);
  TGaxis * x1 = (TGaxis*)m_frame->GetXaxis();
  x1->SetMaxDigits(3);
  TString thetitle=m_frame->GetYaxis()->GetTitle();
  thetitle.Replace(0,6,"Candidates");
  m_frame->GetYaxis()->SetTitle(thetitle);
  m_frame->Draw();
  leg->Draw();
  t->DrawLatex(2.9,m_frame->GetMaximum()*0.91,"LHCb Internal");
  pad2->cd();
  mcor_frame->SetTitle("");
  mcor_frame->GetXaxis()->SetLabelSize(0.06);
  mcor_frame->GetXaxis()->SetTitleSize(0.06);
  mcor_frame->GetYaxis()->SetLabelSize(0.06);
  mcor_frame->GetYaxis()->SetTitleSize(0.06);
  mcor_frame->GetYaxis()->SetTitleOffset(1.3);
  mcor_frame->GetXaxis()->SetTitleOffset(0.9);
  thetitle=mcor_frame->GetYaxis()->GetTitle();
  thetitle.Replace(0,6,"Candidates");
  mcor_frame->GetYaxis()->SetTitle(thetitle);
  mcor_frame->Draw();
  t->DrawLatex(1800,mcor_frame->GetMaximum()*0.91,"LHCb Internal");

  RooPlot *m_resid_frame=x->frame(Title("M(B)"));
  RooPlot *mcor_resid_frame=y->frame(Title("M_{corr}(B)"));


  m_resid_frame->addPlotable(mresid,"B");
  mcor_resid_frame->addPlotable(mcorresid,"B");


  pad4->cd();
  m_resid_frame->SetTitle("");
  m_resid_frame->GetXaxis()->SetTitle("");
  m_resid_frame->GetXaxis()->SetLabelOffset(0.03);
  m_resid_frame->GetYaxis()->SetTitle("Pull");
  m_resid_frame->GetYaxis()->SetRangeUser(-3,3);
  m_resid_frame->GetYaxis()->SetLabelSize(0.12);
  m_resid_frame->GetYaxis()->SetTitleSize(0.12);
  m_resid_frame->GetYaxis()->SetTitleOffset(0.5);
  m_resid_frame->GetXaxis()->SetLabelSize(0.12);
  mresid->SetLineWidth(0);
  m_resid_frame->Draw();
  pad5->cd();
  mcor_resid_frame->SetTitle("");
  mcor_resid_frame->GetXaxis()->SetLabelOffset(0.03);
  mcor_resid_frame->GetXaxis()->SetTitle("");
  mcor_resid_frame->GetYaxis()->SetTitle("Pull");
  mcor_resid_frame->GetYaxis()->SetRangeUser(-3,3);
  mcor_resid_frame->GetYaxis()->SetLabelSize(0.12);
  mcor_resid_frame->GetYaxis()->SetTitleSize(0.12);
  mcor_resid_frame->GetYaxis()->SetTitleOffset(0.5);
  mcor_resid_frame->GetXaxis()->SetLabelSize(0.12);
  mcorresid->SetLineWidth(0);
  mcor_resid_frame->Draw();


  TCanvas *c2 = new TCanvas("c2","c2",1100,600);
  c2->Divide(3,1);
  
  TH2D * signal_proj = (TH2D*) model_signal->createHistogram("signal_proj",*x, RooFit::YVar(*y));
  TH2D * bkg_proj = (TH2D*) partreco_pdf.createHistogram("bkg_proj",*x, RooFit::YVar(*y));
  TH2D * data_proj = (TH2D*) data->createHistogram("data_proj",*x, RooFit::YVar(*y));

  c2->cd(1);
  signal_proj->Draw("LEGO2");
  c2->cd(2);
  bkg_proj->Draw("LEGO2");
  c2->cd(3);
  data_proj->Draw("LEGO2");
  


                                                  

  // w->writeToFile("toyMaker/myWorkspace_roofit.root");

  // cout<<"end"<<endl;

    }
