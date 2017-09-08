void compare_partreco(  ) {


  TFile fnew("/vols/lhcb/palvare1/KeeTuples/MC/presel/PartReco_signalPresel_Kmixture.root");
  TTree *tnew = (TTree*) fnew.Get("DecayTree");
  
  TFile fold("/vols/lhcb/th1011/KeeTuples/MC/HltTOS/B_JpsiX_2012RedoCalo_presel_KisKorPi_HltTrigged_rareBkg.root");
  TTree *told = (TTree*) fold.Get("DecayTree");
  
  TH1F hold("hodl","hold",50,4880,6200);
  TH1F hnew("hnew","hnew",50,4880,6200);

  told->Draw("B_plus_M>>hodl","FidelCastro_TIS");
  tnew->Draw("B_plus_M>>hnew","FidelCastro_TIS");
  hnew.Scale(1.*hold.GetSumOfWeights()/hnew.GetSumOfWeights());
  hnew.SetLineColor(kRed);
  hnew.SetMarkerColor(kRed);
  
  TCanvas *cv = new TCanvas("cv","cv");
  hold.DrawCopy();
  hnew.DrawCopy("same");

  cv->Print("jpsix_kxee_partreco_TIS.pdf");
  
}
