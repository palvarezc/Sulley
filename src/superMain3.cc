#include "usefulFunctions.h"
#include "fitAndSplotKeeDataForTraining.h"
// #include "fitAndSplotPsi2SeeDataForTraining.h"
// #include "fitAndSplotKmumuDataDTFForTraining.h"
// #include "fitAndSplotPsi2SmumuDataDTFForTraining.h"


using namespace std;

int main(int argc, char* argv[])
{

   RooRealVar dummy("dummy", "dummy", -1e10, 1e10);
   RooRealVar yield("yield", "yield", -1e10, 1e10);
   lhcbStyle(0.045);


   FitAndSplotKeeDataForTraining eefit;

   int trigCatMin(0);
   int trigCatMax(2);
   int run(1);
 
   
  
   if(argc == 2)
   {
      if(atoi(argv[1]) < 0)
      {
         eefit.initiateHistoYield();

      }
      return 0;
   }

   if(argc == 3)
   {
      trigCatMin = atoi(argv[1]);
      trigCatMax = atoi(argv[2]);
   }

   if(argc == 4)
   {
      trigCatMin = atoi(argv[1]);
      trigCatMax = atoi(argv[2]);
      run = atoi(argv[3]);
   }

   //*********************************
   //fit Kee MC
   //*********************************


   eefit.minBMass_MC = 4880;  
   eefit.maxBMass_MC = 6200;
   eefit.minBMass_data = 4880;  
   eefit.maxBMass_data = 6200;
   eefit.B_plus_DTF_M_cut = 0.;
 
   eefit.B_plus_M_min_data = 0;
   eefit.B_plus_M_min_MC = 0;
   eefit.B_plus_M_max_data = 10000;
   eefit.B_plus_M_max_MC = 10000;
   
   eefit.weightMC = "";
   eefit.turnOffWeightsSignal = false;
   
   vector<double> BDTcuts;
   vector<string> BDTnames;
   
   eefit.workspaceFileName = "/vols/lhcb/palvare1/RK_analysis/FitWorkspaces/WorkspaceKee_Run"+to_string(run)+".root";
     
   if (run==1)
   {
     BDTnames.push_back("BDT1ETOSRun1");
     BDTcuts.push_back(0.25);
     
     BDTnames.push_back("BDT1Run1");
     BDTcuts.push_back(0.29);
   
     BDTnames.push_back("BDT1Run1");
     BDTcuts.push_back(0.29);

     eefit.tupleMCdir = "/vols/lhcb/palvare1/KeeTuples/MC/presel/";
     eefit.tupleMCname = "Kee_RedoCalo_signalPresel.root";

   }
   
   if (run==2)
   {
     BDTnames.push_back("BDT1ETOSRun2");
     BDTcuts.push_back(0.17);
     
     BDTnames.push_back("BDT1Run2");
     BDTcuts.push_back(0.21);
   
     BDTnames.push_back("BDT1Run2");
     BDTcuts.push_back(0.23);

     eefit.tupleMCdir = "/vols/lhcb/palvare1/KeeTuples/MC/presel/";
     eefit.tupleMCname = "Kee_run2_signalPresel.root";
   }
   


   for(int i(trigCatMin); i<trigCatMax; ++i)
   {
     eefit.plotdir = "/home/hep/palvare1/public_html/KeeMC_fit_Run"+to_string(run)+"/";
     eefit.cutStringMC = "(passTrigCat"+i2s(i)+"&& "+BDTnames[i]+">"+d2s(BDTcuts[i])+")";
     // eefit.cutStringMC = "(K_Kst_ProbNNk>0.2 && e_plus_PIDe > 3 && e_minus_PIDe > 3 && passTrigCat"+i2s(i)+")";

     if (i==0) eefit.weightMC = "FidelCastro_ETOS";
     if (i==1) eefit.weightMC = "FidelCastro_HTOS";
     if (i==2) eefit.weightMC = "FidelCastro_TIS";

     // eefit.prepareMCWorkspace(i, -1);
     eefit.prepareMCWorkspace(i, 0);
     eefit.prepareMCWorkspace(i, 1);
     eefit.prepareMCWorkspace(i, 2);

     eefit.fitMCAuto(i, 0, false);
     eefit.fitMCAuto(i, 1, false);
     eefit.fitMCAuto(i, 2, false);

   }


   return 0;
}
