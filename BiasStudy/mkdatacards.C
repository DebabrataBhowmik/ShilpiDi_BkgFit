/* Performs unbinned fitting of invariant mass shapes and creates text datacards
 * for the "combine" tool from the HiggsAnalysis/CombinedLimit package.
 *
 * Minitrees are taken from output/minitrees.root. The following branches are
 * assumed to exist in each minitree: category, hzg_mass, mcwei (for MC), see
 * fill_events().
 *
 * Signal PDFs are evaluated from minitrees corresponding to existing MC
 * productions (i.e. for certain Higgs mass values). At the same time, signal
 * PDFs for intermediate Higgs masses (for which no MC productions exist) are
 * extrapolated from two PDFs corresponding to closest existing MC productions.
 *
 * The PDF for background description as well as the expected background yield
 * are both obtained from real data.
 *
 * TH1 versions of invariant mass shapes, TF1 versions of produced PDFs as well
 * as normalization factors (RooDouble) which adapt MC scales to real data are
 * written into output/for_visualizations.root.
 *
 * Observed mass points (from real data) together with all obtained PDFs are
 * written into per-channel (eeg or mmg) and per-category (cat1, ..., cat4)
 * files output/datacards/for_datacards_*.root. Datacards themselves are
 * produced as plain text files output/datacards/datacard_*.txt.
 *
 * The directory "./output/datacards" is removed every time this macro is
 * executed.
 *
 * Actual configuration is encoded directly into functions mkdatacard(),
 * normfactor_SMHiggs_8TeV() and mkdatacards().
 *
 * Documentation:
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombination
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise?rev=11
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit?rev=111
 * https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions?rev=18
 * https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationInput?rev=11
 * https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationInputUnbinnedDraft?rev=2
 *
 * Usage: root -b -l -q mkdatacards.C
 */

#include "fitting_functions/fitting_functions.cc"

//#include "myfitting.cc"

//#include "fitting_functions/RooGaussStepBernstein.cxx"
#include "RooMsgService.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "RooWorkspace.h"
#include "TROOT.h"
#include "RooDouble.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "RooPlot.h"

#include <fstream>
#include "TChain.h"

#include "RooAbsPdf.h"
#include "mkdatacards.h"
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooAbsData.h>
#include <RooCBShape.h>
#include <RooProduct.h>
#include <RooAddition.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooBernstein.h>
#include <RooExtendPdf.h>
#include <RooPolynomial.h>
#include <RooAbsPdf.h>

#include <RooExponential.h>
#include <RooFFTConvPdf.h>

#include <RooGaussModel.h>
#include <RooTruthModel.h>
#include <RooDecay.h>

#include "RooAddModel.h"

#include <RooNumConvPdf.h>

#include <fitting_functions/RooStepBernstein.h>

#include <RooGenericPdf.h>
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "CalcSigma.C" ///for eff sigma calc

#include "flashggFinalFit/Background/interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

//#include "RooMultiPdf.h"

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)


//double trigEff_ele = 0.98*0.98*0.997; //leg1 * leg2 * DZ
//double trigEff_mu  = 0.95*0.95*0.97; //leg1 * leg2 * DZ


double trigEff_ele = 1; //leg1 * leg2 * DZ
double trigEff_mu  = 1; //leg1 * leg2 * DZ


double xmin = 115; 
double xmax = 170;

//double mylumi =12900; //pb-1
double mylumi =36460; //pb-1


//bool useLumiWeight = true;
//double lumi_weight = 2000*1000./12900;
bool useLumiWeight = false;
double lumi_weight = 1.0;

//double lumi_weight = 1;

double ngen(const char* channel = "eeg", int p=0, double mass=125);
//void getExpectedEvents(int cat, string channel, double xmin, double xmax, double sigma_away=2,bool useRange=false);
void getExpectedEvents(int cat, string scat, string channel, double xmin, double xmax, double sigma_away=2,bool useRange=false); ///if useRange is set to true, then the events are calculated within that range else within the sig_eff
double getSigEff(const char* filename, bool usewei = true);

//______________________________________________________________________________

////get sig eff for all the categories
double getSigEff(const char* filename, bool usewei)
{
  // open file and get requested minitree
  TFile* fi = TFile::Open(filename);
  if (!fi || fi->IsZombie())
    FATAL("TFile::Open() failed");
  
  //TTree* tree = dynamic_cast<TTree*> (fi->Get("t"));
  TTree* tree = dynamic_cast<TTree*> (fi->Get("minitree"));
  if (!tree) FATAL("TFile::Get() failed");
  
  // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;
   
   
   double nev = 0;
   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
      FATAL("TTree::SetBranchAddress() failed");
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
      FATAL("TTree::SetBranchAddress() failed");

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
     FATAL("TTree::SetBranchAddress() failed");
   
   // NOTE: the X axis variable below (invariant mass) is named exactly as "x"
   // since otherwise unbinned fitting will not work. This is because
   // TF1Wrapper::fX is named as "x"
   
   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     double x = hzg_mass;
     
     ///SJ
     if(hzg_mass<xmin || hzg_mass>xmax) continue;
      
     if (usewei) {
	nev += puwei;
      } else {
       nev += 1;
      }
   }

   return nev;
}

///////
RooAbsData* statAn::fill_events(TH1* hMass, const char* filename, double trigEff, int cat, bool usewei)
{
   /* Fills hMass + a RooDataSet with invariant mass values from given minitree.
    * Returns the filled RooDataSet object.
    *
    * treename = name of minitree in minitrees.root to take;
    * cat = if > 0, take only events from this particular category;
    * usewei = (for MC) if true, attribute weights to all events.
    */

  std::cout<<"Inside fill_events; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  cout<<"trigEff is "<<trigEff<<endl;
   // open file and get requested minitree
  TFile* fi = TFile::Open(filename);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   //TTree* tree = dynamic_cast<TTree*> (fi->Get("t"));
   TTree* tree = dynamic_cast<TTree*> (fi->Get("minitree"));
   if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;

   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
      FATAL("TTree::SetBranchAddress() failed");
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
      FATAL("TTree::SetBranchAddress() failed");

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");

   // NOTE: the X axis variable below (invariant mass) is named exactly as "x"
   // since otherwise unbinned fitting will not work. This is because
   // TF1Wrapper::fX is named as "x"
   RooRealVar x("x", "", 125, xmin, xmax);
   RooRealVar w("w", "", 1);  // NOTE: weights are ignored without this variable

   RooRealVar wLumi("wLumi", "", 10e+6);  // NOTE: weights are ignored without this variable
   RooDataSet* rooMass;
   if (usewei)
      rooMass = new RooDataSet("rooMass", "", RooArgSet(x, w), RooFit::WeightVar(w));
   else
      rooMass = new RooDataSet("rooMass", "", RooArgSet(x));

   
   if(useLumiWeight) rooMass = new RooDataSet("rooMass", "", RooArgSet(x, wLumi), RooFit::WeightVar(wLumi));

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if (cat > 0 && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=5) )
       continue;
     
      x = hzg_mass;

      ///SJ
      if(hzg_mass<xmin || hzg_mass>xmax) continue;
      
      if (usewei) {
	hMass->Fill(hzg_mass, puwei*trigEff*lumi_weight);
	//rooMass->add(RooArgSet(x), puwei*trigEff);
	rooMass->add(RooArgSet(x),puwei*lumi_weight);
	//rooMass->add(RooArgSet(x),1);
	//cout<<"For signal trigEff = "<<trigEff<<endl;
	//hMass->Fill(hzg_mass, trigEff);
	//rooMass->add(RooArgSet(x), trigEff);

      } else {
	hMass->Fill(hzg_mass,lumi_weight);
	rooMass->add(RooArgSet(x),lumi_weight);
      }
   }

   delete tree;
   delete fi;

   return rooMass;
}



///////fill_events for LT///////////////

RooAbsData* statAn::fill_events(TH1* hMass, const char* filename1, const char* filename2, double trigEff1, double trigEff2, int cat, bool usewei)
{
   /* Fills hMass + a RooDataSet with invariant mass values from given minitree.
    * Returns the filled RooDataSet object.
    *
    * treename = name of minitree in minitrees.root to take;
    * cat = if > 0, take only events from this particular category;
    * usewei = (for MC) if true, attribute weights to all events.
    */

  std::cout<<"Inside fill_events; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;
  cout<<"file name is "<<filename1<<" "<<filename2<<endl;
  cout<<"trigEff1 trigEff2 "<<trigEff1<<" "<<trigEff2<<endl;

  TChain *tree = new TChain("minitree");
  tree->Add(filename1);
  tree->Add(filename2);
  
  if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;

   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
     {
       //FATAL("TTree::SetBranchAddress() of category failed");
     }
   
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
     {
       //FATAL("TTree::SetBranchAddress() of hzg_mass failed");
     }

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
     {
       // FATAL("TTree::SetBranchAddress() of mcwei failed");
     }
       

   // NOTE: the X axis variable below (invariant mass) is named exactly as "x"
   // since otherwise unbinned fitting will not work. This is because
   // TF1Wrapper::fX is named as "x"
   //RooRealVar x("x", "", 125, 100, 170);
   RooRealVar x("x", "", 125, xmin, xmax);
   RooRealVar w("w", "", 1);  // NOTE: weights are ignored without this variable
   RooDataSet* rooMass;
   if (usewei)
      rooMass = new RooDataSet("rooMass", "", RooArgSet(x, w), RooFit::WeightVar(w));
   else
      rooMass = new RooDataSet("rooMass", "", RooArgSet(x));

   cout<<"Looping over events now"<<endl;
   
   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       {
	 FATAL("TTree::GetEntry() failed");
       }
     
     if ( (cat > 0 && cat!=6789) && category != cat)
       continue;
     
     //if (cat==0 && !(category>=1 && category<=5) )
     if (cat==0 && !(category>=1 && category<=9) ) ///include the eff of lepton tagged events when i qupte the total eff
       continue;


     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;
     
     x = hzg_mass;
     

     double trigEff = 1;
     string fname(tree->GetCurrentFile()->GetName());
     
     int fstr;
     fstr = fname.find("mu",0);
     if(fstr!=string::npos)
       trigEff = trigEff2;
     
     fstr = fname.find("ele",0);
     if(fstr!=string::npos)
       trigEff = trigEff1;
     
     
     ///SJ
     if(hzg_mass<xmin || hzg_mass>xmax) continue;
     
     if (usewei) {
       
       ///try the condition of using +weighted events only

       if(puwei>0){
	 hMass->Fill(hzg_mass, puwei*trigEff);
	 //rooMass->add(RooArgSet(x), puwei*trigEff);
	 rooMass->add(RooArgSet(x),puwei);
	 
       }//if(puwei>0)
     } else {
       hMass->Fill(hzg_mass);
	rooMass->add(RooArgSet(x));
     }
   }
   
   delete tree;
   //delete fi;

   return rooMass;
}

///////end of fill_events for LT///////////////////


/////////////for ID sys
void statAn::getIDsys(const char* filename, double trigEff , double &sys_lep, double &sys_pho, int cat, bool usewei)
{
  std::cout<<"Inside getIDsys; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  sys_lep = 0;
  sys_pho = 0;

  cout<<"trigEff is "<<trigEff<<endl;
   // open file and get requested minitree
  TFile* fi = TFile::Open(filename);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   //TTree* tree = dynamic_cast<TTree*> (fi->Get("t"));
   TTree* tree = dynamic_cast<TTree*> (fi->Get("minitree"));
   if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;
   float idSys_lep;
   float idSys_pho;

   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
      FATAL("TTree::SetBranchAddress() failed");
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
      FATAL("TTree::SetBranchAddress() failed");

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");
   
   if (tree->SetBranchAddress("idSys_lep", &idSys_lep) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");

   if (tree->SetBranchAddress("idSys_pho", &idSys_pho) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if (cat > 0 && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=5) )
       continue;
     
     //cout<<"idSys_lep idSys_pho "<<idSys_lep << " "<<idSys_pho<<endl;
     sys_lep += idSys_lep*idSys_lep*trigEff*trigEff;
     sys_pho += idSys_pho*idSys_pho*trigEff*trigEff;

     //sys_lep += idSys_lep*idSys_lep;
     //sys_pho += idSys_pho*idSys_pho;

   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   //cout<<"final sys_lep sys_pho "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<endl; 
   
}/////end of ID sys




//______________________________________________________________________________
double statAn::normfactor_SMHiggs_8TeV(const char* channel, int p, int mass)
{
   /* Values were calculated by normfactors_and_theoretical_uncertainties.py.
    *
    * channel = either of "eeg" (for electrons) or "mmg" (for muons);
    * p = production process number (0=ggH, 1=qqH, 2=WH, 3=ZH, 4=ttH);
    * mass = value of SM Higgs mass (from the region 120-150GeV with step 1GeV).
    */


  ///with 2.689... fb-1 of lumi
  /*double norm_eeg[2][11] = {
    {7.07762e-05, 6.32692e-05, 5.77296e-05, 5.34782e-05, 5.01108e-05, 4.72691e-05, 5.44929e-05, 6.33605e-05, 7.45368e-05, 8.89834e-05, 1.08677e-04}, 
    {1.18510e-05, 1.05499e-05, 9.61225e-06, 8.90765e-06, 8.36013e-06, 7.90504e-06, 9.19275e-06, 1.07886e-05, 1.28213e-05, 1.54827e-05, 1.91734e-05}};
  */


  ///with 12.9... fb-1 of lumi - 2016 data
  /*double norm_eeg[2][11] = {
    {6.78877e-04, 7.22026e-04, 7.64049e-04, 8.05428e-04, 8.46358e-04, 8.84950e-04, 7.66563e-04, 6.81145e-04, 6.16182e-04, 5.63947e-04, 5.21204e-04}, 
    {5.63821e-05, 6.05048e-05, 6.45843e-05, 6.86704e-05, 7.27813e-05, 7.67489e-05, 8.01014e-05, 8.33385e-05, 8.64425e-05, 8.92720e-05, 9.19544e-05}};
  */





  ///with 36.46... fb-1 of lumi - 2016 data
  /*double norm_eeg[6][11] = { {2.11478e-03, 2.24872e-03, 2.37852e-03, 2.50662e-03, 2.63312e-03, 2.75275e-03, 2.86643e-03, 2.97609e-03, 3.08076e-03, 3.17547e-03, 3.26395e-03},
			     {1.59354e-04, 1.71007e-04, 1.82536e-04, 1.94085e-04, 2.05704e-04, 2.16918e-04, 2.27130e-04, 2.37076e-04, 2.46701e-04, 2.55598e-04, 2.64125e-04},
			     {3.87179e-05, 2.91123e-05, 2.36827e-05, 2.01498e-05, 1.76907e-05, 1.58447e-05, 1.88806e-05, 2.30382e-05, 2.90243e-05, 3.85125e-05, 5.61572e-05},
			     {2.52410e-05, 1.88637e-05, 1.52930e-05, 1.30202e-05, 1.14052e-05, 1.02038e-05, 1.20848e-05, 1.46049e-05, 1.82022e-05, 2.37422e-05, 3.34815e-05},
			     {4.02512e-05, 3.07220e-05, 2.52567e-05, 2.16839e-05, 1.91539e-05, 1.72218e-05, 2.05429e-05, 2.49733e-05, 3.13497e-05, 4.12918e-05, 5.92335e-05},
			     {2.31706e-04, 2.04980e-04, 1.85018e-04, 1.69730e-04, 1.57444e-04, 1.46963e-04, 1.68559e-04, 1.94717e-04, 2.27816e-04, 2.71005e-04, 3.30462e-04}
  };
  */

  ///with 36.46... fb-1 of lumi - 2016 data - on 21st jan
  /*double norm_eeg[6][11] = {
    {2.11478e-03, 2.24872e-03, 2.37852e-03, 2.50662e-03, 2.63312e-03, 2.75275e-03, 2.86643e-03, 2.97609e-03, 3.08076e-03, 3.17547e-03, 3.26395e-03}, 
    {1.59354e-04, 1.71007e-04, 1.82536e-04, 1.94085e-04, 2.05704e-04, 2.16918e-04, 2.27130e-04, 2.37076e-04, 2.46701e-04, 2.55598e-04, 2.64125e-04}, 
    {3.87179e-05, 2.91123e-05, 2.36827e-05, 2.01498e-05, 1.76907e-05, 1.58447e-05, 2.02061e-05, 2.72942e-05, 4.07792e-05, 7.71365e-05, 5.49160e-04}, 
    {2.52410e-05, 1.88637e-05, 1.52930e-05, 1.30202e-05, 1.14052e-05, 1.02038e-05, 1.20848e-05, 1.46049e-05, 1.82022e-05, 2.37422e-05, 3.34815e-05}, 
    {4.02512e-05, 3.07220e-05, 2.52567e-05, 2.16839e-05, 1.91539e-05, 1.72218e-05, 2.05429e-05, 2.49733e-05, 3.13497e-05, 4.12918e-05, 5.92335e-05}, 
    {2.31706e-04, 2.97737e-04, 3.99534e-04, 5.80330e-04, 9.91324e-04, 2.86481e-03, 2.98677e-03, 3.09848e-03, 3.20597e-03, 3.30575e-03, 3.39988e-03}
  };
  */

  /*double norm_eeg[6][11] = {
    {2.11478e-03, 2.24872e-03, 2.37852e-03, 2.50662e-03, 2.63312e-03, 2.75275e-03, 2.86643e-03, 2.97609e-03, 3.08076e-03, 3.17547e-03, 3.26395e-03},
    {1.59354e-04, 1.71007e-04, 1.82536e-04, 1.94085e-04, 2.05704e-04, 2.16918e-04, 2.27130e-04, 2.37076e-04, 2.46701e-04, 2.55598e-04, 2.64125e-04},
    {3.87179e-05, 2.91123e-05, 2.36827e-05, 2.01498e-05, 1.76907e-05, 1.58447e-05, 1.88806e-05, 2.30382e-05, 2.90243e-05, 3.85125e-05, 5.61572e-05},
    {2.52410e-05, 1.88637e-05, 1.52930e-05, 1.30202e-05, 1.14052e-05, 1.02038e-05, 1.20848e-05, 1.46049e-05, 1.82022e-05, 2.37422e-05, 3.34815e-05},
    {4.02512e-05, 3.07220e-05, 2.52567e-05, 2.16839e-05, 1.91539e-05, 1.72218e-05, 2.05429e-05, 2.49733e-05, 3.13497e-05, 4.12918e-05, 5.92335e-05},
    {2.31706e-04, 2.04980e-04, 1.85018e-04, 1.69730e-04, 1.57444e-04, 1.46963e-04, 1.70590e-04, 2.00069e-04, 2.38749e-04, 2.91748e-04, 3.69797e-04}};
  */

    double norm_eeg[6][11] = {
    {7.04918e-04, 6.65316e-04, 6.32514e-04, 6.05239e-04, 5.82129e-04, 5.61145e-04, 6.32915e-04, 7.16845e-04, 8.16361e-04, 9.35264e-04, 1.08215e-03}, 
    {1.59354e-04, 1.22036e-04, 1.01120e-04, 8.77570e-05, 7.84911e-05, 7.15310e-05, 8.64925e-05, 1.06814e-04, 1.36073e-04, 1.81726e-04, 2.64125e-04}, 
    {3.87179e-05, 2.91123e-05, 2.36827e-05, 2.01498e-05, 1.76907e-05, 1.58447e-05, 1.88806e-05, 2.30382e-05, 2.90243e-05, 3.85125e-05, 5.61572e-05}, 
    {2.52410e-05, 1.88637e-05, 1.52930e-05, 1.30202e-05, 1.14052e-05, 1.02038e-05, 1.20848e-05, 1.46049e-05, 1.82022e-05, 2.37422e-05, 3.34815e-05}, 
    {4.02512e-05, 3.07220e-05, 2.52567e-05, 2.16839e-05, 1.91539e-05, 1.72218e-05, 2.05429e-05, 2.49733e-05, 3.13497e-05, 4.12918e-05, 5.92335e-05}, 
    {2.31706e-04, 2.04980e-04, 1.85018e-04, 1.69730e-04, 1.57444e-04, 1.46963e-04, 1.70590e-04, 2.00069e-04, 2.38749e-04, 2.91748e-04, 3.69797e-04}};


   int index = mass - 120;
   // if (index < 0 || index > 30)
   if (index < 0 || index > 11)
      FATAL("wrong mass given");

   if (TString(channel).Contains("ee"))
      return norm_eeg[p][index];

   return norm_eeg[p][index];
}



//______________________________________________________________________________
void statAn::mkdatacard(const char* channel, int cat)
{
   /* Does everything for particular channel and particular event category.
    *
    * channel = either of "eeg" (for electrons) or "mmg" (for muons);
    * cat = 1, 2, 3 or 4.
    *
    * NOTE: we follow naming conventions of
    * https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
    *
    * NOTE: in this function, memory is leaked (via "new") for coding brevity.
    */

  //double lumi = 2.7*1000;///pb-1
  double lumi =2689.8; //pb-1

   // cat1, cat2, cat3 or cat4
   TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);

   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("hzg_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   //TFile* fo = new TFile("output/for_visualizations.root", "RECREATE");
   /*if (!fo || fo->IsZombie())
      FATAL("TFile::Open() failed");
   */
   
   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_8TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");

   //
   // evaluate PDF of background from real data
   //

   TH1D hDataObs("hDataObs", "", 60, 110, 170);
   hDataObs.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");

   // false = do not add per-event weights
   //TString filename = TString::Format("%s", "input/output_job_2photon_2012abcd_Jan22rereco.root");
   TString filename;
   if(schannel=="eeg") filename = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
   if(schannel=="mmg") filename = TString::Format("%s", "minitree_mu_data_doubleEG_SJ_out.root");
   RooAbsData* dataObs = fill_events(&hDataObs, filename, 1, cat, false);

   //TF1RooGaussStepBernstein bgrfit(110, 170, 4); // 5 = polynomial degree
   //TF1RooPolynomial bgrfit(110, 170, 5);

   ///SJ
   //TF1RooPolynomial bgrfit(110, 170, 5);
   //TF1RooPolynomial bgrfit(110, 170, 3);
   TF1RooPolynomial bgrfit(110, 170, 4);
   //TF1RooGaussStepBernstein bgrfit(110, 170, 4); // 5 = polynomial degree
   //TF1RooGaussStepBernstein bgrfit(110, 170, 5); // 5 = polynomial degree

   //TF1RooPolynomial bgrfit(110, 170, 5);
   RooFitResult *r = bgrfit.FitTo(&hDataObs, dataObs);
   cout<<"Dumping RooFit result ..."<<endl;
   r->Dump();
   // save TH1 and TF1 objects into root file for later visualization
   wd->cd();
   hDataObs.Write("", TObject::kOverwrite);
   bgrfit.GetTF1()->Write("fit_hDataObs", TObject::kOverwrite);
   
   // rename the X axis variable and add observed mass points into the workspace
   dataObs->SetName( ("data_obs_"+catString) );
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass")))
      FATAL("RooWorkspace::import() failed");

   // - rename parameters of the background PDF;
   // - declare parameters of the background PDF to be constants (otherwise the
   //   parameters will be considered as freely floating even without using the
   //   "flatParam" directives below)
   for (int i = 0; i < bgrfit.GetNPar(); i++) {
      const char* name = bgrfit.GetPar(i)->GetName();
      bgrfit.GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s",
                                                channel, cats, name));
      bgrfit.GetPar(i)->setConstant(true);
   }

   // - make normalization utilized by the combine tool to be equal to expected
   //   number of events written into datacards;
   // - add (extended version of) background PDF into the workspace under the
   //   name "pdf_bgr";
   // - add suffix "_bgr" to all PDF's parameters (and to subPDFs, if any);
   // - connect the PDF's X axis variable to the X axis variable of data_obs
   bgrfit.GetPar(0)->setVal(1);
   bgrfit.GetExtPdf()->SetName("pdf");

   /*   if (wspace.import(*bgrfit.GetExtPdf(), RooFit::RenameAllNodes("bgr"),
                                          RooFit::RenameAllVariablesExcept("bgr", "x"),
                                          RooFit::RenameVariable("x", "CMS_hzg_mass")))
   */
   //if (wspace.import(*bgrfit.GetExtPdf(), RooFit::RenameAllNodes("bgr"),

   //SJ
   bgrfit.GetPdf()->SetName("pdf");
   if (wspace.import(*bgrfit.GetPdf(), RooFit::RenameAllNodes(("bgr_"+catString) ),
		     RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		     RooFit::RenameVariable("x", "CMS_hzg_mass")))

      FATAL("RooWorkspace::import() failed");

   // total number of events observed/expected in the region [100, 190]
   int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
   double expectation = bgrfit.GetTF1()->GetParameter(0);

   cout<<"EXPECTATION FROM FITTING IS "<<expectation<<endl;

   ///SJ
   RooRealVar bkgNorm("norm","",expectation);
   bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
   if ( wspace.import(bkgNorm) )
     FATAL("RooWorkspace::import() failed");
   

   r->SetName(("pdf_bgr_"+catString+"_fitresult"));
   if ( wspace.import(*r) )
     FATAL("RooWorkspace::import() fit result failed");

   //RooDouble(expectation).Write("bgr_pdf_norm", TObject::kOverwrite);

   // important consistency check (inconsistency may happen e.g. if by a mistake
   // the mass region in mkminitrees.cc is different from [100, 190]
   if (observation != dataObs->sumEntries())
      FATAL("observation != dataObs->sumEntries()");

   //
   // evaluate signal PDFs for existing MC productions
   //
   
   // - references to signal PDFs for existing MC productions
   // - expected signal yields for all mass points

   cout<<"DOING SIGNAL FIT NOW"<<endl;
   TF1Wrapper* sigfit[2][7]; // [proc][mass]
   double    expected[2][31];

   // Higgs production processes & Higgs masses for existing MC productions
   //const char* proc[2] = {"ggH", "VBFH"};
   //int mass_exist[7] = {120, 125, 130, 135, 140, 145, 150};

   const char* proc[2] = {"ggH", "VBF"};
   int mass_exist[3] = {120, 125, 130};
   

   double sigeff_All[10][100];

   for (int p = 0; p < 2; p++)
   //for (int p = 1; p < 2; p++)
     //for (int m = 0; m < 7; m++) {
     for (int m=0; m<3; m++) {

       if(m>=3) continue;
       
       cout<<"At teh very beginning, m is taking value now "<<m<<endl;
       int mass = mass_exist[m];
       
       cout<<"At the beginning, m, mass_exist is and m is "<<m<<" " <<mass_exist[m]<< " " <<mass<<endl;
       TString hname = TString::Format("hMass_%s_%d", proc[p], mass);
       
       TH1D hMass(hname, "", 120, 110, 170);
       hMass.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
       hMass.SetYTitle("Counts");
       hMass.Sumw2();

         //TString filename = TString::Format("input/output_job_hzg_eeg_dalitz_%s_%d.root",
	 
       TString filename;
       if(channel=="eeg") {
	 filename = TString::Format("minitree_ele_sig%s_%d_out.root",
					    proc[p], mass);
       }


       if(channel=="mmg") {
	 filename = TString::Format("minitree_mu_sig%s_%d_out.root",
					    proc[p], mass);
       }


	 cout<<"proc and mass : Requesting to open "<<proc[p]<<" "<<mass<<" "<<filename<<endl;

       double trigEff = 1;
       if( schannel == "eeg" ) trigEff = trigEff_ele;
       if( schannel == "mmg" ) trigEff = trigEff_mu;

       RooAbsData* rooMass = fill_events(&hMass, filename, trigEff, cat);
	 
     
         // true = use nuisance parameters for the energy scale and resolution
         TF1Wrapper* sfit = new TF1RooCBGaussian(110, 170, true);
	 

         // set initial values and limits
         sfit->GetTF1()->SetParameter(0, hMass.Integral(1, hMass.GetNbinsX()));
         sfit->GetTF1()->SetParameter(1, mass);
         sfit->GetTF1()->SetParLimits(1, 0.85 * mass, 1.15 * mass);

         sfit->FitTo(&hMass, rooMass);
	 
         // save TH1 and TF1 objects into root file for later visualization
         wd->cd();
         hMass.Write("", TObject::kOverwrite);
         sfit->GetTF1()->Write("fit_" + hname, TObject::kOverwrite);

         // save the normalization factor which adapts the MC scale to real
         // data; evaluate expected signal yield
         double norm = normfactor_SMHiggs_8TeV(channel, p, mass); ////returns the xsec
	 //double norm = 1;

	 
         RooDouble(norm).Write("norm_" + hname, TObject::kOverwrite);
         //expected[p][m*5] = norm * sfit->GetTF1()->GetParameter(0); 

	 ///SJ
	 /////for upper limit on xsec////
	 /*double sigEff = hMass.Integral()/ngen(channel,p,mass);
	 double sigExp = sigEff * lumi;
	 expected[p][m*5] = sigExp;
	 sigeff_All[p][m*5] = sigEff;
	 */

	 ///for upper limit on xsec/SMxsec
	 double sigEff = hMass.Integral()/ngen(channel,p,mass);
	 double sigExp = sigEff * lumi;
	 expected[p][m*5] = sigExp;
	 sigeff_All[p][m*5] = sigEff;

	 cout<<"Mass is, Signal eff "<<mass<<" " <<sigEff<<endl;
	 ///SJ

         // unique suffix
         //TString sfx = TString::Format("sig_%s_%d", proc[p], mass[m]);
	 TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);

         // - rename parameters of the PDF;
         // - declare parameters of the PDF to be constants (otherwise the
         //   parameters will be considered as freely floating, and the combine
         //   tool will produce weird results)
         for (int i = 0; i < sfit->GetNPar(); i++) {
            const char* name = sfit->GetPar(i)->GetName();
            sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
                                                     channel, cats, name, sfx.Data()));
            sfit->GetPar(i)->setConstant(true);
         }

	 rooMass->SetName(TString::Format("signaldata_%s_%d_cat%d",proc[p], mass, cat));
	 wspace.import(*rooMass,RooFit::RenameVariable("x", "CMS_hzg_mass"));

         // set names for the nuisance parameters of energy scale and resolution
         sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
         sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));

         // - add signal PDF into the workspace under the name "pdf_sig_...";
         // - add suffix "_sig_..." to all subPDFs;
         // - connect the PDF's X axis variable to the X axis variable of data_obs
         // - connect the nuisance parameters of different signal PDFs together
         sfit->GetPdf()->SetName("pdf");
         if (wspace.import(*sfit->GetPdf(), RooFit::RenameAllNodes(sfx),
                                            RooFit::RenameVariable("x", "CMS_hzg_mass")))
            FATAL("RooWorkspace::import() failed");

         //sigfit[p][m] = sfit;
	 sigfit[p][5*m] = sfit; ///SJ changed
      } // process and mass loops

   //
   // extrapolate signal PDFs for intermediate mass points
   //
  


   
   // use 1 GeV steps
   for (int p = 0; p < 2; p++) // 2 processes
   //for (int p = 1; p < 2; p++) // 2 processes
     //for (int m = 0; m < 6; m++) 
     for (int m = 0; m < 2; m++) 
       for (int k = 1; k <= 4; k++) {
            int mass = 120 + m * 5 + k;

            // neighbours
            //double m1 = mass_exist[m];
            //double m2 = mass_exist[m + 1];

	    ///SJ changed
            double m1 = mass_exist[m];
            double m2 = mass_exist[(m + 1)];
	    
	    cout<<"m1, m2 and their eff are "<<m1<<" "<<m2 <<" Eff "<<sigeff_All[p][5*m]<<" " <<sigeff_All[p][5*(m+1)]<<endl;
            // proportions of first/second neighbour
            double a = (m2 - mass)/(m2 - m1);
            double b = 1 - a;

            // true = use nuisance parameters for the energy scale and resolution
            TF1Wrapper* sfit = new TF1RooCBGaussian(110, 170, true);

            // set values of parameters from linear extrapolation
	    ////y2-y/y2-y1 = x2-x/x2-x1; lets call RHS as a
	    ///gives y = y2(1-a) + ay1
            for (int i = 0; i < sfit->GetNPar() - 2; i++) {
	      //double val1 = sigfit[p][m]    ->GetTF1()->GetParameter(i);
	      //double val2 = sigfit[p][m + 1]->GetTF1()->GetParameter(i);

	      ///SJ changed
               double val1 = sigfit[p][5*m]    ->GetTF1()->GetParameter(i);
               double val2 = sigfit[p][ 5*(m + 1)]->GetTF1()->GetParameter(i);

               sfit->GetPar(i)->setVal(a * val1 + b * val2);
            }

            // evaluate expected signal yield
            double norm = normfactor_SMHiggs_8TeV(channel, p, mass);
            //expected[p][m * 5 + k] = norm * sfit->GetPar(0)->getVal();
       
	    ///extrapolate the eff - SJ
	    expected[p][m * 5 + k] = expected[p][5*m] * a + b * expected[p][5*(m+1)];
	    
	    //cout<<"Setting the element ("<<p<<","<<(m * 5 + k)<<")"
	    sigeff_All[p][m * 5 + k] = expected[p][m * 5 + k]/lumi;
       
	    

            // unique suffix
	    //TString sfx = TString::Format("sig_%s_%d", proc[p], mass);
            TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);

            // - rename parameters of the PDF;
            // - declare parameters of the PDF to be constants (otherwise the
            //   parameters will be considered as freely floating, and the combine
            //   tool will produce weird results)
            for (int i = 0; i < sfit->GetNPar(); i++) {
               const char* name = sfit->GetPar(i)->GetName();
               sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
                                                      channel, cats, name, sfx.Data()));
               sfit->GetPar(i)->setConstant(true);
            }

            // set names for the nuisance parameters of energy scale and resolution
            sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
            sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));

            // - add signal PDF into the workspace under the name "pdf_sig_...";
            // - add suffix "_sig_..." to all subPDFs;
            // - connect the PDF's X axis variable to the X axis variable of data_obs
            // - connect the nuisance parameters of different signal PDFs together
            sfit->GetPdf()->SetName("pdf");
            if (wspace.import(*sfit->GetPdf(), RooFit::RenameAllNodes(sfx),
                                               RooFit::RenameVariable("x", "CMS_hzg_mass")))
               FATAL("RooWorkspace::import() failed");
	       } // extrapolation
   

   wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_8TeV.root",
                                      channel, cats));

   // close output/for_visualizations.root
   delete fo;


   cout<<"WRITING CARD NOW "<<endl;
   //
   // produce datacards
   //
   
   // cuttent UTC time
   TString timestamp(TTimeStamp().AsString());
   timestamp.Remove(timestamp.Last(':'));

   // name of final state
   TString binString = TString::Format("%s_%s_8TeV", channel, cats);
   const char* bin = binString.Data();

   //for (int m = 0; m < 31; m++) {
   for (int m = 0; m < 11; m++) {
      int mass = 120 + m;
      
      cout<<"mass is "<<mass<<endl;
      TString datacard;
      datacard += "# Datacard for the HZg analysis for limit setting\n";
      datacard += "# National Central University, Taiwan\n";
      datacard += "# " + timestamp + " (UTC)\n";
      datacard += TString::Format("# Usage: combine -U -M Asymptotic -m %d datacard.txt\n", mass);
      datacard += "#\n";
      datacard += "imax 1  # number of final states\n";
      datacard += "jmax *  # number of yields given below minus one\n";
      datacard += "kmax *  # number of sources of systematical uncertainties (nuisance parameters)\n";
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("shapes  *         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      //      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_bgr\n", channel, cats);
      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_bgr_%s\n", channel, cats,cats);
      datacard += TString::Format("shapes  data_obs  * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:data_obs_%s\n", channel, cats,cats);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %s\n", bin);
      datacard += TString::Format("bin            %s\n", cats);
      datacard += TString::Format("observation    %d\n", observation);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %-15s %-15s %-15s\n", bin, bin, bin);
      datacard += TString::Format("bin            %-15s %-15s %-15s\n", cats, cats, cats);
      datacard +=                 "process        VBFH            ggH             bgr\n";
      datacard +=                 "process        -1              0               1\n";
      datacard += TString::Format("rate           %-15f %-15f %-15f\n",
				  expected[1][m], expected[0][m], expectation);
      //expected[1][m], expected[0][m], 1);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
   
      // theoretical uncertainties
      /*
      ifstream in("theoretical_uncertainties_SM_Higgs.txt");
      int count = 0; // simple error protection
      char line[4096];

      while (in.good()) {
         in.getline(line, 4096);
         if (!in.eof() && in.fail())
            FATAL("ifstream::getline() failed");

         TString s = line;

         if (s.BeginsWith(TString::Format("%.1f ", (double)mass))) {
            s.Remove(0, s.First(':') + 1); // remove e.g. "155.0   :"
            datacard += s + "\n";
            count++;
         }
      }
      if (count != 7) FATAL("theoretical uncertainties: line count != 7");
      */
   
      // CMS uncertainties.
      //
      // NOTE: naming conventions are from
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
      // TODO: updated numbers?

      datacard += "lumi_8TeV          lnN        1.026          1.026            -\n";
      /*datacard += "eff_8TeV           lnN        1.100          1.100            -\n";
   
      
      // electron vs muon efficiency uncertainties, including uncertainties of
      // triggering efficiencies (it is assumed that the two are not correlated)
      if (TString(channel).Contains("ee"))
         datacard += "CMS_eff_e_8TeV     lnN        1.022          1.022          1.022          1.022          1.022            -\n";
      else
         datacard += "CMS_eff_m_8TeV     lnN        1.038          1.038          1.038          1.038          1.038            -\n";

      // photon efficiency uncertainties: ECAL barrel vs ECAL endcaps
      if (cat == 1 || cat == 2 || cat == 3)
         datacard += "CMS_eff_g_8TeV     lnN        1.006          1.006          1.006          1.006          1.006            -\n";
      else
         datacard += "CMS_eff_g_8TeV     lnN        1.010          1.010          1.010          1.010          1.010            -\n";

      // uncertainty on photon's R9 cut
      if (cat == 1 || cat == 2)
         datacard += "CMS_eff_R9_EB_8TeV lnN        1.050          1.050          1.050          1.050          1.050            -\n";

      float UEPS[] = {1.026, 1.035, 1.018, 1.021};
      datacard += TString::Format("CMS_UEPS_8TeV      lnN          -              -              -            %.3f          1.002            -\n", UEPS[cat-1]);

      float JEC[] = {1.028, 1.032, 1.022, 1.022};
      datacard += TString::Format("CMS_JEC_8TeV       lnN          -              -              -            %.3f          1.001            -\n", JEC[cat-1]);

      float JER[] = {1.010, 1.011, 1.011, 0.014};
      datacard += TString::Format("CMS_JER_8TeV       lnN          -              -              -            %.3f          1.001            -\n", JER[cat-1]);

      if (TString(channel).Contains("ee"))
         datacard += "pdf_PU_mu_8TeV     lnN        1.008          1.008          1.008          1.008          1.008            -\n";
      else
         datacard += "pdf_PU_mu_8TeV     lnN        1.004          1.004          1.004          1.004          1.004            -\n";

      datacard += "------------------------------------------------------------------------------------------------------------\n";
      */


      /*
      // declare uncertainties for the energy scale and resolution
      // (as described by the nuisance parameters above)
      // TODO: errors are set by hand
      datacard += TString::Format("CMS_scale_%s                                param          1     0.05\n", channel);
      datacard += TString::Format("CMS_res_%s                                  param          1     0.01\n", channel);
      
      
      // declare parameters of the background PDF to be freely floating
      for (int i = 0; i < bgrfit.GetNPar(); i++) {
	TString nameStr = TString::Format("%s_bgr", bgrfit.GetPar(i)->GetName());
         datacard += TString::Format("%-44s flatParam\n", nameStr.Data());
      }
      */


      cout<<"Writing data card finally"<<endl;
      // make datacard file
      FILE* out = fopen(TString::Format("output/datacards/datacard_hzg_%s_%s_8TeV_%d.txt",
                                        //channel, cats, mass[m]).Data(), "w");
					channel, cats, mass).Data(), "w");
      if (!out) FATAL("fopen() failed");
      if (fputs(datacard.Data(), out) < 0) FATAL("fputs() failed");
      if (fclose(out) != 0) FATAL("fclose() failed");

      cout<<"Wrote finally"<<endl;
   } // mass loop
   

   cout<<"=========================For cat========================="<<cat<<endl;
   for (int p = 0; p < 2; p++) {// 2 processes
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 

	 //cout<<sigeff_All[p][m * 5 + k]<<", ";
	 cout<<"mass is "<<mass<<" and eff is "<<sigeff_All[p][m * 5 + k]<<endl;
	 
       }//for (int k = 1; k <= 4; k++)
     }
     cout<<"}"<<endl;
   }

   ///Draw the data and the fit to it with error bands
   int W = 800;
   int H = 600;
   
   int H_ref = 600;
   int W_ref = 800;
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;
   
   TCanvas* c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   TFile *fout = TFile::Open(TString::Format("output/datacards/for_datacards_hzg_%s_%s_8TeV.root",channel, cats));
   
   if(fout==NULL) cout<<"Root file is null "<<endl;
   
   RooWorkspace *ws = (RooWorkspace *)fout->Get("hzg_workspace"); 
   //RooRealVar *mass = ws->var("CMS_hzg_mass");
   RooRealVar x("x", "", 125, 110, 170);
   //RooRealVar x("x", "", 125, 130, 170);
   RooPlot *plot = x.frame();
   //RooAbsData *data = ws->data(Form("data_obs_cat%s",cats)); 

   //RooAbsPdf *extbkg = ws->pdf(Form("pdf_bgr_cat%s",cats));

   if(ws==NULL)
     cout<<"WARNING!!!ws "<<endl;

   if(&x==NULL)
     cout<<"WARNING!!!mass "<<endl;

   if(plot==NULL)
     cout<<"WARNING!!!plot "<<endl;

   if(dataObs==NULL)
     cout<<"WARNING!!!data "<<endl;

   if(bgrfit.GetPdf()==NULL)
     cout<<"WARNING!!!extbkg "<<endl;
       
   dataObs->plotOn(plot);
   bgrfit.GetPdf()->plotOn(plot);
   bgrfit.GetExtPdf()->plotOn(plot,RooFit::LineColor(kRed));


   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   r->Print();

   ///error bands
   ///1sigma
   /*bgrfit.GetPdf()->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*r,1),RooFit::FillColor(kOrange-6));
   bgrfit.GetPdf()->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*r,2),RooFit::FillColor(kOrange-10));
   */

   bgrfit.GetExtPdf()->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*r,1),RooFit::FillColor(kYellow-4));
   bgrfit.GetExtPdf()->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*r,2),RooFit::FillColor(kGreen-4));
      
   char *outfilename = new char[100];
   char dirName[100] = "plots";
   

   string scat(cats);
   plot->Draw();
   sprintf(outfilename,"%s/%s.gif",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   sprintf(outfilename,"%s/%s.C",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
}


///fitting function
//void mkdatacards::mybkgfit(string model, int poldeg){
void statAn::mybkgfit(double xmin, double xmax, int icat){

  rfit = NULL;

  //RooRealVar  fX("x", "", xmin, xmax); 


   PdfModelBuilder pdfsModel;
   pdfsModel.setObsVar(fX);

  RooExtendPdf *bgrfit_ext = new RooExtendPdf();

  cout<<"Model is "<<model<<endl;
  if(model=="RooGaussStepBernstein"){

    cout<<"inside RooGaussStepBernstein"<<endl;
    //mean = new RooRealVar(Form("mean_%d",poldeg), "mean", 10, -100, 100);
    mean = new RooRealVar(Form("mean_%d_model%s",poldeg,model.c_str()), "mean", 0);

    //if(poldeg==4 && icat==4) mean = new RooRealVar(TString::Format("mean_%d_model%s",poldeg,model.c_str()), "mean", 90.0,0, 500.);
 
    
    if(icat==1) sigma = new RooRealVar(Form("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 3, 0.01, 100);
    else sigma = new RooRealVar(Form("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 6, 0.01, 100);

    if(poldeg==4 && icat==4) sigma = new RooRealVar(TString::Format("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 2., 0, 500);
    
    stepval = new RooRealVar(Form("stepval_%d_model%s",poldeg,model.c_str()), "Stepval",70, 10, 500);

    if(icat==4 && poldeg==4) stepval = new RooRealVar(Form("stepval_%d_model%s",poldeg,model.c_str()), "Stepval",70, 10, 500);
    
    coeflist = new RooArgList;
    
    
    for(int ipol=0; ipol<=poldeg; ipol++){
      pol[ipol] = new RooRealVar(TString::Format("coef%d_%d_model%s",ipol,poldeg,model.c_str()), TString::Format("coef%d",ipol),5 , 0,50);
      //pol[ipol] = new RooRealVar(TString::Format("coef%d_%d_model%s",ipol,poldeg,model.c_str()), TString::Format("coef%d",ipol),5 , -50,50);


      if(ipol==0){
	
	if(icat==1) pol[ipol]->setVal(5);
	else  pol[ipol]->setVal(1);
	pol[ipol]->setConstant(kTRUE);
      }

      sqpol[ipol] = new RooFormulaVar(TString::Format("sq_coef%d_%d_model%s",ipol,poldeg,model.c_str()),"@0*@1",RooArgList(*pol[ipol],*pol[ipol]));

      
      coeflist->add(*pol[ipol]);
      //coeflist->add(*sqpol[ipol]);
    }    
    

    

    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+10);
    
    cout<<"now defining bgrfit1"<<endl;
    RooGaussStepBernstein *bgrfit1 = new RooGaussStepBernstein(Form("GaussStepBernstein_%d",poldeg), "", *fX, *mean, *sigma, *stepval,*coeflist);
    
    bgrfit_ext = new  RooExtendPdf(Form("fit_%d_model%s",poldeg,model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    

    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(-1),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    
    if(rfit==NULL) cout<<"inside the fitting function, r is NULL"<<endl;
    
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }
  
  
  else if(model=="RooPolynomial"){
    

    coeflist = new RooArgList;
    cout<<"inside RooPolynomial"<<endl;
    cout<<"Power is "<<poldeg<<endl;
    RooAbsPdf *bgrfit1 = pdfsModel.getBernstein(Form("Bernstein_%d_model%s",poldeg,model.c_str()),poldeg);
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+6);    
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+10);    
    
    ///used till 11april, 2017 - v5 of the AN
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+10);    

    //used from 11th april, 2017
    //fPar[0] = new RooRealVar(TString::Format("norm_%d_model%s",poldeg,model.c_str()), "", 200, 0, 1e+10);    

    fPar[0] = new RooRealVar(TString::Format("norm_%d_model%s",poldeg,model.c_str()), "", 500, 0, 1e+10);    

    //RooExtendPdf bgrfit("fit", "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(Form("fit_%d_model%s",poldeg,model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    if(rfit==NULL) cout<<"inside the fitting function, r is NULL"<<endl;

    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }//if(model=="RooPolynomial")


  ///Gaus * exp
  else if(model=="GausconvExp"){
    

    cout<<"inside GausconvExp"<<endl;
    
    mean = new RooRealVar(TString::Format("mean_model%s",model.c_str()), "mean", 120., 90., 150.);
    sigma = new RooRealVar(TString::Format("sigma_model%s",model.c_str()), "sigma", 1., 0.01, 10);
    
    lambda = new RooRealVar(TString::Format("lambda_model%s",model.c_str()), "slope", 5, 0., 50.);
    
    RooGaussModel *gauss = new RooGaussModel(TString::Format("gauss_model%s",model.c_str()),"gaussian PDF", *fX, *mean, *sigma);

    RooDecay *bgrfit1 = new RooDecay(TString::Format("expo_model%s",model.c_str()), "exponential PDF", *fX, *lambda, *gauss, RooDecay::SingleSided);

    // Construct landau (x) gauss
    //cout<<"now defining bgrfit1"<<endl;
    //RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("exg_model%s",model.c_str()),"landau (X) gauss",*fX,*expo,*gauss) ;
      
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit("fit", "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_model%s",model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(-1),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }//if(model=="GausconvExp")

  ///sech * Exp
  else if(model=="SechconvExp"){
    
    cout<<"inside SechconvExp"<<endl;

    
    
    //mean = new RooRealVar("mean", "mean", 115.);
    mean = new RooRealVar(TString::Format("mean_model%s",model.c_str()), "mean", 100.0,0, 500.);
    sigma = new RooRealVar(TString::Format("sigma_model%s",model.c_str()), "sigma", 10., 0, 500);
    
    lambda = new RooRealVar(TString::Format("lambda_model%s",model.c_str()), "slope", 6, 0.0, 100.);
    //alphabkg = new RooRealVar("alphabkg", "alphaSec", 100, 0., 200.);
    alphabkg = new RooRealVar(TString::Format("alphabkg_model%s",model.c_str()), "alphaSec", 85.0, 0., 1000.);
    

    
    //RooGenericPdf *turnOn = new RooGenericPdf(TString::Format("turnOn_model%s",model.c_str()),"Turn On", "exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)", RooArgList(*fX, *mean, *sigma));
    //RooGenericPdf *turnOn = new RooGenericPdf(TString::Format("turnOn_model%s",model.c_str()),"Turn On", "2*exp(-(@0-@1)/@2)/(1.0+exp(-2*(@0-@1)/@2))", RooArgList(*fX, *mean, *sigma));
    RooGenericPdf *turnOn = new RooGenericPdf(TString::Format("turnOn_model%s",model.c_str()),"Turn On", "exp(-(@0-@1)/@2)/(1.0+exp(-2*(@0-@1)/@2))", RooArgList(*fX, *mean, *sigma));
    //RooGenericPdf *bgrfit1 = new RooGenericPdf(TString::Format("bgrfit1_model%s",model.c_str()),"Turn On", "2*exp(-(@0-@1)/@2)/(1.0+exp(-2*(@0-@1)/@2))", RooArgList(*fX, *mean, *sigma));
    //RooGenericPdf *bgrfit1 = new RooGenericPdf(TString::Format("bgrfit1_model%s",model.c_str()),"Turn On", "exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)", RooArgList(*fX, *mean, *sigma));
    
    RooGenericPdf *tail = new RooGenericPdf(TString::Format("tail_model%s",model.c_str()), "tail PDF", "1e-20 + (@0 > @1)*(exp(-@0/@2))", RooArgList(*fX, *alphabkg, *lambda));
    //RooGenericPdf *tail = new RooGenericPdf(TString::Format("tail_model%s",model.c_str()), "tail PDF", "(exp(-@0/@2))", RooArgList(*fX, *lambda));
    //RooGenericPdf *bgrfit1 = new RooGenericPdf(TString::Format("bgrfit1_model%s",model.c_str()), "tail PDF", "1e-20 + (@0 > @1)*(exp(-@0/@2))", RooArgList(*fX, *alphabkg, *lambda));

    
    //fX->setBins(1000,"fft");
    // Construct sech (x) Exp
    cout<<"now defining bgrfit1"<<endl;
    //RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("SechxExp_model%s",model.c_str()),"Sech (X) Exp",*fX,*tail,*turnOn) ;
    RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("SechxExp_model%s",model.c_str()),"Sech (X) Exp",*fX,*turnOn,*tail) ;
    //RooNumConvPdf *bgrfit1 = new RooNumConvPdf(TString::Format("SechxExp_model%s",model.c_Str()),"Sech (X) Exp",*fX,*tail,*turnOn);
    bgrfit1->setShift(0,21.);
  
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit(TString::Format("fit_model%s",model.c_str()), "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_model%s",model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }//if(model=="SechconvExp")

  
  /////
  if(model=="GausPow") {


    cout<<"inside GausPow"<<endl;
    

    mean = new RooRealVar(TString::Format("mean_model%s",model.c_str()), "mean", 100., -250., 200.);
    //mean = new RooRealVar("mean", "mean", 0.);
    sigma = new RooRealVar(TString::Format("sigma_model%s",model.c_str()), "sigma", 10, 0.01, 100);
    alphabkg = new RooRealVar(TString::Format("alphabkg_model%s",model.c_str()), "alphabkg", 105., 0., 500.);
    betabkg = new RooRealVar(TString::Format("betabkg_model%s",model.c_str()), "betabkg", 1, 0., 100.);
    

    /*sigma = new RooRealVar(TString::Format("sigma_model%s",model.c_str()), "sigma", 55, 0.01, 300);
    alphabkg = new RooRealVar(TString::Format("alphabkg_model%s",model.c_str()), "alphabkg", 68.5, 0., 500.);
    betabkg = new RooRealVar(TString::Format("betabkg_model%s",model.c_str()), "betabkg", 6, 0., 50.);
    */
			     
    RooGaussModel *gauss = new RooGaussModel(TString::Format("gauss_model%s",model.c_str()),"gaussian PDF", *fX, *mean, *sigma);
    RooGenericPdf *tail = new RooGenericPdf(TString::Format("tail_model%s",model.c_str()), "Power PDF", "1e-20 + (@0 > @1)*((@0)^(-@2))", RooArgList(*fX,*alphabkg,*betabkg));
    //RooGenericPdf *tail = new RooGenericPdf(TString::Format("tail_model%s",model.c_str()), "Power PDF", "((@0)^(-@1))", RooArgList(*fX,*betabkg));

    // Construct Gaus (x) Powerlaw
    //cout<<"now defining bgrfit1"<<endl;



    //fX->setBins(10000);
    //RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("gxp_model%s",model.c_str()),"PowerLaw (X) gauss",*fX,*tail,*gauss) ;
    RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("gxp_model%s",model.c_str()),"PowerLaw (X) gauss",*fX,*gauss,*tail) ;
    //bgrfit1->setInterpolationOrder(3);
    bgrfit1->setShift(0,40);

    /*RooRealVar *gm1frac = new RooRealVar(TString::Format("gm1frac_model%s",model.c_str()),"fraction of gm1",0.5,0,1) ;
    cout<<"defining rooaddmodel"<<endl;
    RooAddModel *bgrfit1 = new RooAddModel(TString::Format("gxp_model%s",model.c_str()),"PowerLaw (X) gauss",RooArgList(*tail,*tail), *gm1frac) ;
    cout<<"defined rooaddmodel"<<endl;
    */
    
    
    //RooNumConvPdf *bgrfit1 = new RooNumConvPdf(TString::Format("gxp_model%s",model.c_str()),"PowerLaw (X) gauss",*fX,*tail,*gauss) ;
      
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit("fit", "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_model%s",model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


    


  }//if(model=="GausPow")



  ///Sech * Bern
  else if(model=="SechBern"){
    
    cout<<"inside SechBern degree is "<<poldeg<<endl;
    
    //mean = new RooRealVar("mean", "mean",0.);
    if(poldeg==4) mean = new RooRealVar(TString::Format("mean_%d_model%s",poldeg,model.c_str()), "mean", 105.0,0, 500.);
    else mean = new RooRealVar(TString::Format("mean_%d_model%s",poldeg,model.c_str()), "mean", 100.0,0, 500.);

    //if(poldeg==4 && icat==4) mean = new RooRealVar(TString::Format("mean_%d_model%s",poldeg,model.c_str()), "mean", 90.0,0, 500.);

    if(poldeg==3) sigma = new RooRealVar(TString::Format("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 7., 0, 500);
    if(poldeg==2) sigma = new RooRealVar(TString::Format("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 10., 0, 500);
    if(poldeg==4) sigma = new RooRealVar(TString::Format("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 2., 0, 500);

    if(poldeg==4 && icat==4) sigma = new RooRealVar(TString::Format("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 1., 0, 500);
    if(poldeg==3 && icat==2) sigma = new RooRealVar(TString::Format("sigma_%d_model%s",poldeg,model.c_str()), "sigma", 1., 0, 500);
    
    stepbkg = new RooRealVar( TString::Format("stepbkg_%d_model%s",poldeg,model.c_str()),"stepbkg",0.1,0.,10.);

    if(icat==4 && poldeg==4) stepbkg = new RooRealVar(Form("stepbkg_%d_model%s",poldeg,model.c_str()), "Stepval",0.1, 0, 500);
    //if(icat==2 && (poldeg==3||poldeg==4) ) stepbkg = new RooRealVar(Form("stepbkg_%d_model%s",poldeg,model.c_str()), "Stepval",0.1, 0, 500);
    if(icat==2 && poldeg==4 ) stepbkg = new RooRealVar(Form("stepbkg_%d_model%s",poldeg,model.c_str()), "Stepval",0.1, 0, 500);
    if(icat==2 && poldeg==3) stepbkg = new RooRealVar(Form("stepbkg_%d_model%s",poldeg,model.c_str()), "Stepval",0.1, 0, 500);

    coeflist = new RooArgList;
    
    for(int ipol=0; ipol<=poldeg; ipol++){
      pol[ipol] = new RooRealVar(TString::Format("coef%d_%d_model%s",ipol,poldeg,model.c_str()), TString::Format("coef%d",ipol),5 , 0,100);

      if(ipol==0){
	if(icat==1) pol[ipol]->setVal(5);
	else  pol[ipol]->setVal(3);
	pol[ipol]->setConstant(kTRUE);
      }
      
      coeflist->add(*pol[ipol]);
    }//for(int ipol=0; ipol<=poldeg; ipol++)

    
    
    RooGenericPdf *turnOn = new RooGenericPdf(TString::Format("turnOn_%d_model%s",poldeg,model.c_str()),"Turn On", "exp(-(@0-@1)/@2)/(1.0+exp(-2*(@0-@1)/@2))", RooArgList(*fX, *mean, *sigma));


    RooStepBernstein *tail = new RooStepBernstein(TString::Format("tail_%d_model%s",poldeg,model.c_str()),"bernstein pol", *fX,*stepbkg,*coeflist);

    // Construct sech (x) Exp
    cout<<"now defining bgrfit1"<<endl;
    //RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("SechxBern_%d_model%s",poldeg,model.c_str()),"Sech (X) Bern",*fX,*tail,*turnOn) ;
    RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("SechxBern_%d_model%s",poldeg,model.c_str()),"Sech (X) Bern",*fX,*turnOn,*tail) ;
    if(icat==1) bgrfit1->setShift(0,20.);
    else bgrfit1->setShift(0,15.);

    if(icat==4 && poldeg==4) bgrfit1->setShift(0,5.);
  
    //fPar[0] = new RooRealVar(TString::Format("norm_%d_model%s",poldeg,model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_%d_model%s",poldeg,model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit(TString::Format("fit_%d_model%s",poldeg,model.c_str()), "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_%d_model%s",poldeg,model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }//if(model=="SechBern3")


  /////sech * Pow
  else if(model=="SechPow"){
    
    cout<<"inside SechPow"<<endl;
    
    mean = new RooRealVar(TString::Format("mean_model%s",model.c_str()), "mean", 100.0,0, 500.);
    sigma = new RooRealVar(TString::Format("sigma_model%s",model.c_str()), "sigma", 10., 0, 500);

    alphabkg = new RooRealVar(TString::Format("alphabkg_model%s",model.c_str()), "alphabkg", 105., 0., 500.);
    betabkg = new RooRealVar(TString::Format("betabkg_model%s",model.c_str()), "betabkg", 1, 0., 100.);

    
    RooGenericPdf *turnOn = new RooGenericPdf(TString::Format("turnOn_%d_model%s",poldeg,model.c_str()),"Turn On", "exp(-(@0-@1)/@2)/(1.0+exp(-2*(@0-@1)/@2))", RooArgList(*fX, *mean, *sigma));
    
    RooGenericPdf *tail = new RooGenericPdf(TString::Format("tail_model%s",model.c_str()), "Power PDF", "1e-20 + (@0 > @1)*((@0)^(-@2))", RooArgList(*fX,*alphabkg,*betabkg));
    
    // Construct sech (x) Pow
    cout<<"now defining bgrfit1"<<endl;
    //RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("SechxBern_%d_model%s",poldeg,model.c_str()),"Sech (X) Bern",*fX,*tail,*turnOn) ;
    RooFFTConvPdf *bgrfit1 = new RooFFTConvPdf(TString::Format("SechxBern_%d_model%s",poldeg,model.c_str()),"Sech (X) Bern",*fX,*turnOn,*tail) ;
    bgrfit1->setShift(0,15.);
  
    //fPar[0] = new RooRealVar(TString::Format("norm_%d_model%s",poldeg,model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_%d_model%s",poldeg,model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit(TString::Format("fit_%d_model%s",poldeg,model.c_str()), "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_%d_model%s",poldeg,model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }//if(model=="SechPow")

  
  ////Exp
  else if(model=="Exp"){
    

    cout<<"inside Exp"<<endl;
    cout<<"Power is "<<poldeg<<endl;

    RooAbsPdf *bgrfit1 = pdfsModel.getExponentialSingle(Form("expo_model%s",model.c_str()),poldeg);

    //lambda = new RooRealVar(TString::Format("lambda_model%s",model.c_str()), "slope", -1, -100, 0.);
    //RooExponential *bgrfit1 = new RooExponential(TString::Format("expo_model%s",model.c_str()), "exponential PDF", *fX, *lambda);

    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit("fit", "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_model%s",model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }//if(model=="Exp")


  ///power law
  if(model=="Pow") {


    cout<<"inside Pow"<<endl;
    
    cout<<"Power is "<<poldeg<<endl;

    RooAbsPdf *bgrfit1 = pdfsModel.getPowerLawSingle(Form("bgrfit1_model%s",model.c_str()),poldeg);
    
      
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit("fit", "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_model%s",model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


    


  }//if(model=="Pow")


  ///Laurent
  if(model=="Laurent") {


    cout<<"inside Laurent"<<endl;
    cout<<"Power is "<<poldeg<<endl;

    RooAbsPdf *bgrfit1 = pdfsModel.getLaurentSeries(Form("bkgfit1_model%s",model.c_str()),poldeg); 
    
      
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+6);
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200., 0., 1e+10);
    
    //RooExtendPdf bgrfit("fit", "", bgrfit1, *fPar[0], "our_window");
    bgrfit_ext = new  RooExtendPdf(TString::Format("fit_model%s",model.c_str()), "", *bgrfit1, *fPar[0], "our_window");
    
 
    
    
    cout<<"Now fitting with bgrfit"<<endl;
    
    rfit = bgrfit_ext->fitTo(*dataObs,
			     RooFit::PrintLevel(3),
			     //RooFit::Range(xmin, xmax),
			     RooFit::NumCPU(4),
			     //RooFit::SumW2Error(kTRUE),
			     RooFit::Save(true)
			     
			     );
    
    cout<<"Now fitted "<<endl;
    cout<<"pdf set to bgrfit"<<endl;

    bgrfit = bgrfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    bgrfit->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


    


  }//if(model=="Pow")


  else{
    cout<<"Problem!!!!Model given is "<<model<<endl;
    cout<<"Error!!! provide either RooGaussStepBernstein or RooPolynomial as model"<<endl;
    //return 0;
    
  }
  
  
  
}

std::string statAn::getbkgfitParName(int ipar){


  string name[50];  
  if(model=="RooGaussStepBernstein"){   
    cout<<"inside RooGaussStepBernstein"<<endl;
    

    name[0]= Form("mean_%d",poldeg);
    name[1] = Form("sigma_%d",poldeg);
    name[2] = Form("stepval_%d",poldeg);
    
    for(int ipol=1; ipol<=poldeg; ipol++){  
      name[ipol+2] = TString::Format("coef%d",ipol);
    }
  }//if(model=="RooGaussStepBernstein")
  

  
  if(model=="RooPolynomial"){   
    cout<<"inside RooPolynomial"<<endl;
    
    for(int ipol=0; ipol<=poldeg; ipol++){  
      name[ipol] = TString::Format("coef%d",ipol);
    }
  }//if(model=="RooPolynomial")
  

  if(model=="GausconvExp"){

    cout<<"inside GausconvExp"<<endl;
    
    name[0]= Form("mean_%d",poldeg);
    name[1] = Form("sigma_%d",poldeg);
    name[2] = Form("lambda_%d",poldeg);
    
  }//

  return name[ipar];
}

int statAn::getNbkgPar(){

  if(model=="RooGaussStepBernstein"){
    
    return poldeg+3;
  }


  if(model=="RooPolynomial"){
    
    return poldeg;
  }


  if(model=="GausconvExp"){
    
    return 3;
  }

}

/////signal fit function

void statAn::mysigfunc(string model){


  cout<<"Model is "<<model<<endl;
  if(model=="RooCBxGaus"){

    cout<<"inside RooCBxGaus"<<endl;
    
    pdf1 = new RooCBShape ("subpdf1", "", *fX, *sigmean, *sigsigma, *alpha, *power);
    //pdf2 = new RooGaussian("subpdf2", "", *fX, *mean2, *sigma2);
    pdf2 = new RooGaussian("subpdf2", "", *fX, *sigmean, *sigma2);
    sigfit1 = new RooAddPdf("CBGaussian", "", *pdf1, *pdf2, *frac);


    cout<<"----values inside mysigfunc-----"<<endl;
    //cout<<"meana : sigmaa : alpha : power : mean2 : sigma2 : frac "<<sigmean->getVal()<<" "<<sigsigma->getVal()<<" "<<alpha->getVal()<<" "<<power->getVal()<<" "<<mean2->getVal()<<" "<<sigma2->getVal()<<" "<<frac->getVal()<<endl;
    cout<<"meana : sigmaa : alpha : power : mean2 : sigma2 : frac "<<sigmean->getVal()<<" "<<sigsigma->getVal()<<" "<<alpha->getVal()<<" "<<power->getVal()<<" "<<sigma2->getVal()<<" "<<frac->getVal()<<endl;
    
  }// if(model=="RooCBxGaus"



}



void statAn::mysigfit(string model, double xmin, double xmax, bool nuisance, int np, int nm, int icat, string schannel){

  RooExtendPdf *sigfit_ext = new RooExtendPdf();
  
  //get the idea about mean and sigma of the dataset
  double mean_ini = rooMassSig[np][nm]->mean(*fX);
  double sigma_ini = rooMassSig[np][nm]->sigma(*fX);
  
  cout<<"mean_ini and sigma_ini "<<mean_ini<<" "<<sigma_ini<<endl;

  mean_ini = 125;
  sigma_ini = 2;

  cout<<"Model is "<<model<<endl;
  if(model=="RooCBxGaus"){

    cout<<"inside RooCBxGaus"<<endl;
    




    mean1   = new RooRealVar(Form("sig_mean1_chan%d_m%d_cat%d",np,nm,icat),"", 125, xmin, xmax);
    //mean1   = new RooRealVar(Form("sig_mean1_chan%d_m%d_cat%d",np,nm,icat),"", mean_ini, xmin, xmax);
    //sigma1  = new RooRealVar(Form("sig_sigma1%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 8);
    sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 10);
    //if( (icat==1 && np==0 && nm==0) )     sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 6, 0.5, 10);
    //if(icat==1 && np==0) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.3, 20);
    //if(icat==1 && np==0) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 20);
    
    ///25thjan
    //if( (icat==1 && (np==0||np==1) ) || (icat==2&&np==0) || (icat==3 && np==1) || (icat==4&&np==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 20);

    //if( (icat==1 && (np==0||np==1) ) || (icat==3 && np==1) || (icat==4&&np==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 20);
    //if( (icat==1 && (np==0||np==1) ) || (icat==3 && np==1) || (icat==4&&np==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 25);
    if( (icat==1 && (np==0||np==1) ) || (icat==3 && np==1) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 25);

    if( (icat==1 && np==1  && nm==10 && schannel=="mmg") || (icat==4&&np==0&&nm==10&&schannel=="eeg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.1, 10);


    if( (icat==1 && np==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 20);

    //if( (icat==1 && (np==1) ) || (icat==4&&np==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 20);
    //if( (icat==2&&np==0)  && schannel=="eeg" && nm == 5 ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.3, 0.1, 5);
    /*if( (icat==2&&np==0)  && nm == 5 ){
      cout<<"DONE THE CHANGE FOR cat : mass : np "<<icat<<" "<<nm<<" "<<np<<endl;
      sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.01, 8);

    }
    */


    ///22jan  
    //if( (icat==1 && np==0 && nm==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 20);
    
    //if( (icat==3 && np==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 1, 20);
    

    //sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", sigma_ini, 0.1, 10);
    //sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 20);
    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.5, 10);
    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.01, 20);
    alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.01, 50);

    //if( (icat==2&&np==0)  && nm == 5 )  alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",0.1, 0.01, 50);
    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.01, 50);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.5, 10);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.01, 20);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.001, 30);

    ///26th jan
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.001, 70);
    power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.001, 200);

    if( (icat==3&&np==0)  && schannel=="eeg" && nm == 10)  power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.001, 30);
    //if( (icat==2&&np==0)  && schannel=="eeg" && nm == 5)  power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.01, 50);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.01, 30);
    //delta21 = new RooRealVar(Form("delta21_chan%d_m%d_cat%d",np,nm,icat), "",0, -50, 50); // mean2 = mean1 + delta21
    //delta21 = new RooRealVar(Form("delta21_chan%d_m%d_cat%d",np,nm,icat), "",0, -10, 10); // mean2 = mean1 + delta21
    //s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 1, 30);       // sigma2/sigma1
    //s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);       // sigma2/sigma1
    //s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 0.1, 30);       // sigma2/sigma1

    s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 0.1, 30);       // sigma2/sigma1

    
    ///25th jan
    //if( (icat==3 && np==0 && nm==10) || icat==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 2, 30);       // sigma2/sigma1
    //if( (icat==3 && np==0 && nm==10) || icat==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 2, 30);
    if( (icat==3 && np==0 && nm==10) ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 2, 30);
    
    //if( (icat==1 && np==0 && nm==0) ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",0.2, 0.1, 30);       // sigma2/sigma1
    //s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 0.1, 20);       // sigma2/sigma1
    frac    = new RooRealVar(Form("frac_chan%d_m%d_cat%d",np,nm,icat), "",0.9, 0, 1);

    //mean2  = new RooAddition(Form("sig_mean2_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*mean1, *delta21));
    sigma2 = new RooProduct(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*sigma1, *s21));
    
    sigmean  = mean1;
    sigsigma = sigma1;
    
    
    if (nuisance) {
    
      scale = new RooRealVar(Form("scale_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.5, 1.5);
      resol = new RooRealVar(Form("resolution_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.5, 1.5);
      
      scale->setConstant(true);
      resol->setConstant(true);
      
      sigmean   = new RooProduct(Form("sig_mean1_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*mean1, *scale));
      //mean2  = new RooProduct(Form("sig_mean2_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*mean2, *scale));
      sigsigma  = new RooProduct(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*sigma1, *resol));
      sigma2 = new RooProduct(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*sigma2, *resol));
    }
    
    
    
    
    mysigfunc(model);
    

    //fPar_sig[0] = new RooRealVar( Form("norm_chan%d_m%d_cat%d",np,nm,icat), "", 200, 0, 1e+6);    
    fPar_sig[0] = new RooRealVar( Form("norm_chan%d_m%d_cat%d",np,nm,icat), "", 200, 0, 1e+10);    
    RooExtendPdf *sigfit_ext = new  RooExtendPdf("sigfit", "", *sigfit1, *fPar_sig[0], "our_window");
    
    
    cout<<"inside mysigfit np and nm "<<np<<" "<<nm<<endl;
    
    rsig[np][nm] = sigfit_ext->fitTo(*rooMassSig[np][nm],
				     RooFit::PrintLevel(2),
				     RooFit::Range(xmin, xmax),
				     RooFit::NumCPU(4),
				     RooFit::SumW2Error(kTRUE),
				     RooFit::Save(true)
				     //	     ,RooFit::Extended(true)
				     );


    
    
    if(rsig[np][nm]==NULL) cout<<"inside the signal fitting function, sfit is NULL"<<endl;
    
    
    cout<<"Now signal fitted "<<endl;
    cout<<"pdf set to signalfit"<<endl;

    sigfit[np][nm] = sigfit1;
    cout<<"====================Inside hte function printing bgrfit1===================="<<endl;
    sigfit[np][nm]->Print();
    cout<<"-----------------------Printed the function---------------------------"<<endl;


  }
  
  
  else{
    
    cout<<"Error!!! provide either RooCBxGaus model"<<endl;
    //return 0;
    
  }
  
  
  
}










void statAn::deleteSigPar(bool nuisance){

  delete mean1;
  delete sigma1;
  delete alpha;
  delete power;
  //delete delta21;
  delete s21;
  delete frac;

  //delete mean2;
  delete sigma2;
  
  delete sigmean;
  delete sigsigma;
    
    
  if (nuisance) {
    
    delete scale;
    delete resol;
      
    }
    
  delete fPar_sig[0];

}



/////////////////////////fitting function


///////// interpolate
void statAn::siginterpolate(string model, int np, int nm, int nk, double a, double b, bool nuisance, int icat){



  if(model=="RooCBxGaus"){
    
    cout<<"inside siginterpolate with model "<<model<<endl;
    RooRealVar* mean1_fitresult1 = (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("sig_mean1_chan%d_m%d_cat%d",np,nm,icat));
    //cout<<"Found sig1_mean1"<<endl;
    RooRealVar* mean1_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("sig_mean1_chan%d_m%d_cat%d",np,nm+5,icat));

    //cout<<"Found sig2_mean1"<<endl;

    double mean1_fr =  a * mean1_fitresult1->getVal() + b * mean1_fitresult2->getVal();
    
    cout<<"mean1, mean2 and calculated mean1_fr"<<mean1_fitresult1->getVal() <<" "<<mean1_fitresult2->getVal()<< " "<<mean1_fr<<endl;
    //cout<<"mean1, mean2 and calculated mean1_fr"<<mean1_fitresult1->getVal() <<" "<<mean1_fitresult2->getVal()<< " "<<mean1_fr<<endl;

    
    //RooRealVar* delta21_fitresult1 = (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("delta21_chan%d_m%d_cat%d",np,nm,icat));
    //cout<<"Found delta21_1"<<endl;

    //RooRealVar* delta21_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("delta21_chan%d_m%d_cat%d",np,nm+5,icat));
    //cout<<"Found delta21_2"<<endl;
    
    //double delta21_fr =  a * delta21_fitresult1->getVal() + b * delta21_fitresult2->getVal();
    //cout<<"delta1, delta2 and calculated delta1_fr"<<delta21_fitresult1->getVal() <<" "<<delta21_fitresult2->getVal()<< " "<<delta21_fr<<endl;
    
    //double delta21_fr = 0;
    

        
    RooRealVar* sigma1_fitresult1 = (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat));
    //cout<<"Found sig1_sigma1"<<endl;
    RooRealVar* sigma1_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm+5,icat));
    //cout<<"Found sig1_sigma2"<<endl;
    double sigma1_fr =  a * sigma1_fitresult1->getVal() + b * sigma1_fitresult2->getVal();
    //cout<<"calculated sigma1_fr"<<endl;

    cout<<"sigma1, sigma2 and calculated sigma1_fr"<<sigma1_fitresult1->getVal() <<" "<<sigma1_fitresult2->getVal()<< " "<<sigma1_fr<<endl;


    RooRealVar* s21_fitresult1 = (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("s21_chan%d_m%d_cat%d",np,nm,icat));
    RooRealVar* s21_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("s21_chan%d_m%d_cat%d",np,nm+5,icat));
    double s21_fr =  a * s21_fitresult1->getVal() + b * s21_fitresult2->getVal();
    //cout<<"calculated s21_fr"<<endl;

    cout<<"s21_1, s21_2 and calculated s21_fr"<<s21_fitresult1->getVal() <<" "<<s21_fitresult2->getVal()<< " "<<s21_fr<<endl;
    
    RooRealVar* alpha_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("alpha_chan%d_m%d_cat%d",np,nm,icat));
    RooRealVar* alpha_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("alpha_chan%d_m%d_cat%d",np,nm+5,icat));
    double alpha_fr =  a * alpha_fitresult1->getVal() + b * alpha_fitresult2->getVal();

    //cout<<"Found alpha_fr"<<endl;
    cout<<"alpha_1, alpha_2 and calculated alpha_fr"<<alpha_fitresult1->getVal() <<" "<<alpha_fitresult2->getVal()<< " "<<alpha_fr<<endl;

    RooRealVar* power_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("power_chan%d_m%d_cat%d",np,nm,icat));
    RooRealVar* power_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("power_chan%d_m%d_cat%d",np,nm+5,icat));
    double power_fr =  a * power_fitresult1->getVal() + b * power_fitresult2->getVal();
    
    //cout<<"Found power_fr"<<endl;
    cout<<"power_1, power_2 and calculated power_fr"<<power_fitresult1->getVal() <<" "<<power_fitresult2->getVal()<< " "<<power_fr<<endl;

    RooRealVar* frac_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("frac_chan%d_m%d_cat%d",np,nm,icat));
    RooRealVar* frac_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("frac_chan%d_m%d_cat%d",np,nm+5,icat));
    double frac_fr =  a * frac_fitresult1->getVal() + b * frac_fitresult2->getVal();

    //cout<<"Found everything"<<endl;

    cout<<"frac_1, frac_2 and calculated frac_fr"<<frac_fitresult1->getVal() <<" "<<frac_fitresult2->getVal()<< " "<<frac_fr<<endl;

    /*
    ///try mean2 and sigma2////
    ///sig_mean2_chan1_m10_cat1_sig_VBF_130_cat1
    RooRealVar* mean2_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("sig_mean2_chan%d_m%d_cat%d",np,nm,icat));
    RooRealVar* mean2_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("sig_mean2_chan%d_m%d_cat%d",np,nm+5,icat));
    double mean2_fr =  a * mean2_fitresult1->getVal() + b * mean2_fitresult2->getVal();

    RooRealVar* sigma2_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm,icat));
    RooRealVar* sigma2_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm+5,icat));
    double sigma2_fr =  a * sigma2_fitresult1->getVal() + b * sigma2_fitresult2->getVal();
    */

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
    //sigmean  = new RooRealVar(Form("sig_mean1_chan%d_m%d_cat%d",np,nm+nk,icat),"", mean1_fr);
    //sigsigma  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm+nk,icat),"", sigma1_fr);


    mean1  = new RooRealVar(Form("sig_mean1_chan%d_m%d_cat%d",np,nm+nk,icat),"", mean1_fr);

    sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm+nk,icat),"", sigma1_fr);


    sigmean = mean1;
    sigsigma = sigma1;

    //delta21 = new RooRealVar(Form("delta21_chan%d_m%d_cat%d",np,nm+nk,icat), "",delta21_fr); 
    s21 = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm+nk,icat), "",s21_fr); 
    alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm+nk,icat), "", alpha_fr);
    power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm+nk,icat), "",power_fr);

    frac    = new RooRealVar(Form("frac_chan%d_m%d_cat%d",np,nm+nk,icat), "",frac_fr);

    //mean2  = new RooAddition(Form("sig_mean2_chan%d_m%d_cat%d",np,nm+nk,icat), "", RooArgList(*mean1, *delta21));
    sigma2 = new RooProduct(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm+nk,icat), "", RooArgList(*sigma1, *s21));

    //mean2  = new RooAddition(Form("sig_mean2_chan%d_m%d_cat%d",np,nm+nk,icat), "", RooArgList(*sigmean, *delta21));
    //sigma2 = new RooProduct(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm+nk,icat), "", RooArgList(*sigsigma, *s21));

    
    
    if(nuisance){

      
      RooRealVar* scale_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("scale_chan%d_m%d_cat%d",np,nm,icat));
      RooRealVar* scale_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("scale_chan%d_m%d_cat%d",np,nm+5,icat));
      
      //cout<<"pointer to scales "<< scale_fitresult1<<" "<< scale_fitresult2<<endl;
      double scale_fr =  a * scale_fitresult1->getVal() + b * scale_fitresult2->getVal();
      
      cout<<"scale_1, scale_2 and calculated scale_fr"<<scale_fitresult1->getVal() <<" "<<scale_fitresult2->getVal()<< " "<<scale_fr<<endl;
    
      RooRealVar* resolution_fitresult1 =  (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("resolution_chan%d_m%d_cat%d",np,nm,icat));
      RooRealVar* resolution_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("resolution_chan%d_m%d_cat%d",np,nm+5,icat));
      
      //cout<<"pointer to resolution "<< resolution_fitresult1<<" "<< resolution_fitresult2<<endl;
      double resolution_fr =  a * resolution_fitresult1->getVal() + b * resolution_fitresult2->getVal();
      
      cout<<"resolution_1, resolution_2 and calculated resolution_fr"<<resolution_fitresult1->getVal() <<" "<<resolution_fitresult2->getVal()<< " "<<resolution_fr<<endl;

      scale = new RooRealVar(Form("scale_chan%d_m%d_cat%d",np,nm+nk,icat), "",scale_fr);
      resol = new RooRealVar(Form("resolution_chan%d_m%d_cat%d",np,nm+nk,icat), "",resolution_fr);

      sigmean   = new RooProduct(Form("sig_mean1_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*mean1, *scale));
      sigsigma  = new RooProduct(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat), "", RooArgList(*sigma1, *resol));
      //mean2  = new RooProduct(Form("sig_mean2_chan%d_m%d_cat%d",np,nm+nk,icat), "", RooArgList(*mean2, *scale));
      sigma2 = new RooProduct(Form("sig_sigma2_chan%d_m%d_cat%d",np,nm+nk,icat), "", RooArgList(*sigma2, *resol));
    }
    
    
    mysigfunc(model);


  }//if(model=="RooCBxGaus")
}




///set the signal parameters constant at this point
void statAn::setSigParConst(string model, bool nuisance){
  
  if(model=="RooCBxGaus"){
    
    cout<<"inside RooCBxGaus"<<endl;
    
    mean1->setConstant(true);
    sigma1->setConstant(true);
    alpha->setConstant(true);
    power->setConstant(true);
    //delta21->setConstant(true);
    s21->setConstant(true);
    frac->setConstant(true);

    /*
    mean2->setConstant(true);
    sigma2->setConstant(true);
    
    sigmean->setConstant(true);
    sigsigma->setConstant(true);
    */
    
    if (nuisance) {
    
      scale->setConstant(true);
      resol->setConstant(true);
    }
    
  }//if(model=="RooCBxGaus")


}

//////////////////////////////////////

////simple mkdatacards for chk
void statAn::simplemkdatacard(const char* channel, int cat)
{

  string sigmodel = "RooCBxGaus";
  bool _nuisance = false;
  //bool _nuisance = true;
  
  //double lumi = 2.7*1000;///pb-1
  double lumi = 2689.8; //pb-1
  // cat1, cat2, cat3 or cat4
  TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);

   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("hzg_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   
   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_8TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");

   //
   // evaluate PDF of background from real data
   //

   //double xmin = 110;
   double xmin = 100;
   double xmax = 170;

   TH1D hDataObs("hDataObs", "", 60, 100, 170);
   hDataObs.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");
   

   TString filename;
   if(schannel=="eeg") filename = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
   if(schannel=="mmg") filename = TString::Format("%s", "minitree_mu_data_doubleMu_SJ_out.root");
   dataObs = fill_events(&hDataObs, filename, 1, cat, false);
   
   ///SJ
  
   // make background function

   rfit = new RooFitResult();
   
   bgrfit = NULL;
   //poldeg = 4;
   poldeg = 3;
   
   //RooGaussStepBernstein *bgrfit = new RooGaussStepBernstein();
   fX = new RooRealVar("x", "", xmin, xmax); 
   //bgrfit = mybkgfit("RooGaussStepBernstein", poldeg);
   mybkgfit(xmin, xmax,cat);
   //////try importing here//////////////////


   ///set the parameters constant  - is it the right way? - CHECK
   mean->setConstant(true);
   sigma->setConstant(true);
   stepval->setConstant(true);
   for(int ipol=1; ipol<=poldeg; ipol++){

     pol[ipol]->setConstant(true); 
   }


   cout<<"===DONE WITH THE FITTING"<<endl;
   bgrfit->Print();

   //////////////////////////////////////////////
   
   cout<<"in mkdatacards, address of r is "<<rfit<<endl;
   if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;

   cout<<"Dumping RooFit result ..."<<endl;
   rfit->Dump();
   // save TH1 and TF1 objects into root file for later visualization
   wd->cd();
   // rename the X axis variable and add observed mass points into the workspace
   dataObs->SetName( ("data_obs_"+catString) );
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass")))
     FATAL("RooWorkspace::import() failed");
   

   //SJ
   cout<<"now setting bgrfit pdf to name of pdf "<<endl;
   //bgrfit->GetPar(0)->setConstant(true);
   bgrfit->SetName("pdf");

   cout<<"address of bgrfit "<<bgrfit<<endl;
   
   cout<<"set the name to name of pdf "<<endl;
   if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		     RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		     RooFit::RenameVariable("x", "CMS_hzg_mass")))
   
      
   cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;

   
   cout<<"now saved the bkg pdf"<<endl;
   
   
   rfit->SetName(("pdf_bgr_"+catString+"_fitresult"));
  if ( wspace.import(*rfit) )
     FATAL("RooWorkspace::import() fit result failed");


  /*
  RooArgList bkgpar_list;
   bkgpar_list.add(*mean);
   bkgpar_list.add(*sigma);
   bkgpar_list.add(*stepval);
   bkgpar_list.add(*coeflist);
   RooArgSet bkgpar_set(bkgpar_list);
   

   cout<<"Values of parameters after fixing ..."<<endl;
   cout<<"Mean : sigma : stepval : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<stepval->getVal()<<endl;
  */
   
   // total number of events observed/expected in the region [100, 190]
   int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
   //double expectation = bgrfit->getValV(&bkgpar_set);
   RooArgSet obs(*fX);
   //double expectation = bgrfit->getValV(&obs);
   

   RooAbsReal* intbkgfit1 = bgrfit->createIntegral(*fX) ; 
   //double expectation = intbkgfit1->getVal();
   //double expectation = bgrfit1->getVal();
   double expectation = fPar[0]->getVal();
   cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
   
   cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;



   RooRealVar bkgNorm("norm","",expectation);
   bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
   if ( wspace.import(bkgNorm) )
     FATAL("RooWorkspace::import() failed");
   



  /////start with the signal now /////////

   cout<<"DOING SIGNAL FIT NOW"<<endl;
   double    expected[nprocess][nmass];
   const char* proc[nprocess] = {"ggH", "VBF"};
   int mass_exist[nsig] = {120, 125, 130};
   double sigeff_All[nprocess][nmass];
   double sigNraw[nprocess][nmass];
   double sigeff_ggf_vbf[nmass];

   ///systematcs ID
   double sigSysonExp_lep[nprocess][nmass];
   double sigSysonExp_pho[nprocess][nmass];

   double sigSysonEff_lep[nprocess][nmass];
   double sigSysonEff_pho[nprocess][nmass];


   for (int p = 0; p < nprocess; p++){
     for (int m=0; m<nsig; m++) {
       
       cout<<"========================At teh very beginning, m is taking value now================"<<m<<endl;
       
       if((int)m >= (int)3){
	 cout<<"----inside >=3 loop, Will continue"<<endl;
	 continue;
	 cout<<"No, wait!!!! not continued"<<endl;
       }
       

       rooMassSig[p][m*5] = NULL;
       sigfit[p][m*5] = NULL;
       rsig[p][m*5] = NULL;
       
       
       int mass = mass_exist[m];
       
       cout<<"At the beginning, m, mass_exist is and m is "<<m<<" " <<mass_exist[m]<< " " <<mass<<endl;
       TString hname = TString::Format("hMass_%s_%d", proc[p], mass);
       
       TH1D hMass(hname, "", 120, xmin,xmax);
       hMass.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
       hMass.SetYTitle("Counts");
       hMass.Sumw2();
       
       
       TString filename;
       if(channel=="eeg") {
	 filename = TString::Format("minitree_ele_sig%s_%d_out.root",
				    proc[p], mass);
       }
       

       if(channel=="mmg") {
	 filename = TString::Format("minitree_mu_sig%s_%d_out.root",
				    proc[p], mass);
       }
       
       
       cout<<"proc and mass : Requesting to open "<<proc[p]<<" "<<mass<<" "<<filename<<endl;

       double trigEff = 1;
       if( schannel == "eeg" ) trigEff = trigEff_ele;
       if( schannel == "mmg" ) trigEff = trigEff_mu;
  
       cout<<"Calling for fill_events for isgn, trigEff is "<<trigEff<<endl;
       RooAbsData* rooMass = fill_events(&hMass, filename, trigEff, cat);
       
       rooMassSig[p][5*m] = rooMass;
       
       cout<<"=================Channel, cat "<<channel<<" "<<cat<<endl;
       cout<<"===========fitting for mass =================="<<mass<<endl;
       mysigfit(sigmodel,xmin, xmax, _nuisance, p, 5*m, cat, schannel);
       
       setSigParConst(sigmodel, _nuisance);
     
       ///set the parameters const at this point for these mass points
       
    
       // save the normalization factor which adapts the MC scale to real
       // data; evaluate expected signal yield
       double norm = normfactor_SMHiggs_8TeV(channel, p, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double norm = 1;
       
     
       ///SJ
       double sigEff = hMass.Integral()/ngen(channel,p,mass);
       
       double sigExp = sigEff * lumi; /// use this for limit on xsec*BR
       //double sigExp = hMass.Integral()*norm;///use this for limit on xsec/SMxsec
       
       ///systematics on ID
       double sys_lep = 0;
       double sys_pho = 0;
       getIDsys(filename,  trigEff, sys_lep, sys_pho, cat, true);
       //sigSysonExp_lep[p][m*5] = sqrt(sys_lep)*norm/sigExp;
       //sigSysonExp_pho[p][m*5] = sqrt(sys_pho)*norm/sigExp;
       

       sigSysonExp_lep[p][m*5] = sqrt(sys_lep)/hMass.Integral();
       sigSysonExp_pho[p][m*5] = sqrt(sys_pho)/hMass.Integral();
       
       ///sys on Eff is not always right. Only when sigExp=sigEff*lumi, then it is right. 
       ///Only the above sigSysonExp_lep numbers are right. These are on Eff   as well. So sigSysonExp_lep can be used on Exp and Eff
       sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen(channel,p,mass)*sigExp);
       sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen(channel,p,mass)*sigExp);
       
       
       cout<<"sqrt(lepsys) and sqrt(phosys) ; sigExp : "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<hMass.Integral()<<endl;
       
       cout<<" p : mass  : lepID sys EXP : EFF "<<p<<" "<<mass<<" "<<sigSysonExp_lep[p][m*5]<<" "<<sigSysonEff_lep[p][m*5] <<endl;
       cout<<" p : mass  : phoID sys EXP : EFF "<<p<<" "<<mass<<" "<<sigSysonExp_pho[p][m*5]<<" "<<sigSysonEff_pho[p][m*5] <<endl;
       

       expected[p][m*5] = sigExp;
       sigeff_All[p][m*5] = sigEff;
       
       sigNraw[p][m*5] = hMass.Integral();
       
       cout<<"Mass is, Signal eff "<<mass<<" " <<sigEff<<endl;

       cout<<"===p is, channel, cat, expectation=== "<<p<<" "<<channel<<" "<<cat <<" "<<expected[p][m*5]<<endl;
       ///SJ
       
       // unique suffix
       //TString sfx = TString::Format("sig_%s_%d", proc[p], mass[m]);
       TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);

       // - rename parameters of the PDF;
       // - declare parameters of the PDF to be constants (otherwise the
       //   parameters will be considered as freely floating, and the combine
       //   tool will produce weird results)
       

       /*for (int i = 0; i < sfit->GetNPar(); i++) {
	 const char* name = sfit->GetPar(i)->GetName();
	 sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
	 channel, cats, name, sfx.Data()));
	 sfit->GetPar(i)->setConstant(true);
         }
       */
     

       ///CHECK - set signal parameters to const
       
       rooMassSig[p][m*5]->SetName(TString::Format("signaldata_%s_%d_cat%d",proc[p], mass, cat));
       wspace.import(*rooMassSig[p][m*5],RooFit::RenameVariable("x", "CMS_hzg_mass"));
       
       
       cout<<"Inside the signal mass loop for filling the data set"<<endl;
       cout<<"=========Printing the dataset information========="<<endl;
       rooMassSig[p][m*5]->Print();
       cout<<"=========Printed the dataset information========="<<endl;

       // set names for the nuisance parameters of energy scale and resolution
       /*sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
       sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));
       */

       
         // - add signal PDF into the workspace under the name "pdf_sig_...";
         // - add suffix "_sig_..." to all subPDFs;
         // - connect the PDF's X axis variable to the X axis variable of data_obs
         // - connect the nuisance parameters of different signal PDFs together

       cout<<"setting the signal pdf name "<<endl;
       
       sigfit[p][5*m]->SetName("pdf");
       //sfit->SetName(Form("pdf_sig_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*sigfit[p][5*m], RooFit::RenameAllNodes(sfx),
			 RooFit::RenameVariable("x", "CMS_hzg_mass")))
	 FATAL("RooWorkspace::import() failed");
       
       
       cout<<"saved the signal pdf"<<endl;

       rsig[p][m*5]->SetName(Form("fitresult_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*rsig[p][m*5]))
	 FATAL("RooWorkspace::import() failed");
       
       cout<<"saved the signal fit result"<<endl;
       



       //delete sfit;
       
     } // process and mass loops
   }//for (int p = 0; p < nprocess; p++)

   cout<<"nprocess : nsig : "<<nprocess<<" "<<nsig<<endl;

   
   for (int m = 0; m < 3; m++) {

     int mass = mass_exist[m];   
     double sigEff = (sigNraw[0][m * 5] + sigNraw[1][m * 5])/( ngen(channel,0,mass) + ngen(channel,1,mass) );   
     sigeff_ggf_vbf[m*5] = sigEff;
   }


  
   ///interpolation
   // use 1 GeV steps
   for (int p = 0; p < nprocess; p++){ // 2 processes
   //for (int p = nprocess-1; p >=0; p--){ // 2 processes
     for (int mm = 0; mm< nsig-1; mm++){ 
       for (int k = 1; k <= 4; k++) {
	 
	 int mass = 120 + mm * 5 + k;
	 
	 cout<<""<<endl;
	 cout<<"p : mm : k : mass "<<p<<" "<<mm<<" "<<k<<" "<<mass<<endl;
	 
	 ///SJ changed
	 double m1 = mass_exist[mm];
	 double m2 = mass_exist[(mm + 1)];
	 
	 cout<<"m1, m2 and their eff are "<<m1<<" "<<m2 <<" Eff "<<sigeff_All[p][5*mm]<<" " <<sigeff_All[p][5*(mm+1)]<<endl;
	 // proportions of first/second neighbour
	 double a = (m2 - mass)/(m2 - m1);
	 double b = 1 - a;
	 
	 //RooRealVar* par1_fitresult = (RooRealVar*) fitresult->floatParsFinal()->find("par1"); 
	 
	 
	 // set values of parameters from linear extrapolation
	 ////y2-y/y2-y1 = x2-x/x2-x1; lets call RHS as a
	 ///gives y = y2(1-a) + ay1
	 //for (int i = 0; i < sfit->GetNPar() - 2; i++) {
	 //double val1 = sigfit[p][m]    ->GetTF1()->GetParameter(i);
	 //double val2 = sigfit[p][m + 1]->GetTF1()->GetParameter(i);
	 	    
       
	 // evaluate expected signal yield
	 double norm = normfactor_SMHiggs_8TeV(channel, p, mass);
	 //expected[p][m * 5 + k] = norm * sfit->GetPar(0)->getVal();
	 
	 ///extrapolate the eff - SJ
	 expected[p][mm * 5 + k] = expected[p][5*mm] * a + b * expected[p][5*(mm+1)];
	 
	 //cout<<"Setting the element ("<<p<<","<<(m * 5 + k)<<")"
	 //sigeff_All[p][mm * 5 + k] = expected[p][mm * 5 + k]/lumi;
	 sigeff_All[p][mm * 5 + k] = sigeff_All[p][mm*5] * a + b * sigeff_All[p][5*(mm+1)];

	 sigeff_ggf_vbf[mm * 5 + k] = sigeff_ggf_vbf[mm*5] * a + b * sigeff_ggf_vbf[5*(mm+1)];

	 //systematics
	 sigSysonExp_lep[p][mm * 5 + k] = sigSysonExp_lep[p][mm*5] * a + b * sigSysonExp_lep[p][5*(mm+1)];
	 sigSysonExp_pho[p][mm * 5 + k] = sigSysonExp_pho[p][mm*5] * a + b * sigSysonExp_pho[p][5*(mm+1)];

	 sigSysonEff_lep[p][mm * 5 + k] = sigSysonEff_lep[p][mm*5] * a + b * sigSysonEff_lep[p][5*(mm+1)];
	 sigSysonEff_pho[p][mm * 5 + k] = sigSysonEff_pho[p][mm*5] * a + b * sigSysonEff_pho[p][5*(mm+1)];

	 


	 cout<<"========starting to interpolate now for mass "<<mass<<"================="<<endl;
	 cout<<"====masses taken are "<<m1<<" and "<<m2<<endl;
	 siginterpolate(sigmodel, p, mm*5, k, a,b, _nuisance, cat);
	 sigfit[p][5*mm+k] = sigfit1; //filled in mysigfunc
	 
	 setSigParConst(sigmodel, _nuisance);
       
	 cout<<"========END OF interpolate ================="<<endl;

	 // - add signal PDF into the workspace under the name "pdf_sig_...";
	 // - add suffix "_sig_..." to all subPDFs;
	 // - connect the PDF's X axis variable to the X axis variable of data_obs
	 // - connect the nuisance parameters of different signal PDFs together
	 
	 TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);
	 sigfit[p][5*mm+k]->SetName("pdf");
	 if (wspace.import(*sigfit[p][5*mm+k], RooFit::RenameAllNodes(sfx),
			   RooFit::RenameVariable("x", "CMS_hzg_mass")))
	   FATAL("RooWorkspace::import() failed");
	 

	 /*
	 rsig[p][mm*5+k]->SetName(Form("fitresult_%s_%d_cat%d",proc[p], mass, cat));
	 if (wspace.import(*rsig[p][mm*5+k]))
	   FATAL("RooWorkspace::import() failed");
	 */

	 // unique suffix
	 
	 
	 // - rename parameters of the PDF;
	 // - declare parameters of the PDF to be constants (otherwise the
	 //   parameters will be considered as freely floating, and the combine
	 //   tool will produce weird results)
	 //for (int i = 0; i < sfit->GetNPar(); i++) {
	 //const char* name = sfit->GetPar(i)->GetName();
	 //sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
	 //channel, cats, name, sfx.Data()));
	 //sfit->GetPar(i)->setConstant(true);
	 //}
	 
	 // set names for the nuisance parameters of energy scale and resolution
	 //sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
	 //sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));


	 //deleteSigPar(_nuisance);
	 
       } // extrapolation
     }//for (int mm = 0; mm< nsig; mm++)
   }//for (int p = 0; p < nprocess; p++)


            


   

  //////end of the signal now//////////////





   
   

   cout<<"Writing hte work space now"<<endl;
   wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_8TeV.root",
                                      channel, cats));
   
   




   // close output/for_visualizations.root
   delete fo;
   

      cout<<"WRITING CARD NOW "<<endl;
   //
   // produce datacards
   //
   
   // cuttent UTC time
   TString timestamp(TTimeStamp().AsString());
   timestamp.Remove(timestamp.Last(':'));

   // name of final state
   TString binString = TString::Format("%s_%s_8TeV", channel, cats);
   const char* bin = binString.Data();

   //for (int m = 0; m < 31; m++) {
   for (int m = 0; m < 11; m++) {
      int mass = 120 + m;
      


      sigSysonExp_lep[0][m] = 1 + sigSysonExp_lep[0][m];
      sigSysonExp_pho[0][m] = 1 + sigSysonExp_pho[0][m];

      sigSysonExp_lep[1][m] = 1 + sigSysonExp_lep[1][m];
      sigSysonExp_pho[1][m] = 1 + sigSysonExp_pho[1][m];
      
      cout<<"mass is "<<mass<<endl;
      TString datacard;
      datacard += "# Datacard for the HZg analysis for limit setting\n";
      datacard += "# National Central University, Taiwan\n";
      datacard += "# " + timestamp + " (UTC)\n";
      datacard += TString::Format("# Usage: combine -U -M Asymptotic -m %d datacard.txt\n", mass);
      datacard += "#\n";
      datacard += "imax 1  # number of final states\n";
      datacard += "jmax *  # number of yields given below minus one\n";
      datacard += "kmax *  # number of sources of systematical uncertainties (nuisance parameters)\n";
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //      datacard += TString::Format("shapes  *         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_$PROCESS_%d\n", channel, cats, mass);
      //datacard += TString::Format("shapes  *         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      
      //      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_bgr\n", channel, cats);
      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_bgr_%s\n", channel, cats,cats);
      datacard += TString::Format("shapes  data_obs  * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:data_obs_%s\n", channel, cats,cats);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %s\n", bin);
      datacard += TString::Format("bin            %s\n", cats);
      //datacard += TString::Format("observation    %d\n", observation);
      datacard += TString::Format("observation    %d\n", -1);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %-15s %-15s %-15s\n", bin, bin, bin);
      datacard += TString::Format("bin            %-15s %-15s %-15s\n", cats, cats, cats);
      datacard +=                 "process        VBFH            ggH             bgr\n";
      datacard +=                 "process        -1              0               1\n";
      datacard += TString::Format("rate           %-15f %-15f %-15f\n",
				  //			  expected[1][m], expected[0][m], expectation);
				  expected[1][m], expected[0][m], 1.0);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
   
      // theoretical uncertainties
      
      ifstream in("theoretical_uncertainties_SM_Higgs.list");
      int count = 0; // simple error protection
      char line[10000];

      while (in.good()) {
	in.getline(line, 10000);
	if (!in.eof() && in.fail())
	   FATAL("ifstream::getline() failed");

         TString s = line;
	 cout<<""<<endl;
	 cout<<"line is "<<line<<endl;
		 
	 //cout<<TString::Format("%.1f ", (double)mass)<<endl;
         //if(s.BeginsWith(TString::Format("%.1f ", (double)mass))) {
	 if(s.BeginsWith(TString::Format("%.0f ", (double)mass))) {  ///needs to be changed when we have intervals of 0.5 GeV
	   cout<<"added this line "<<endl;
            s.Remove(0, s.First(':') + 1); // remove e.g. "155.0   :"
            datacard += s + "\n";
            count++;
         }
      }
      cout<<"count is "<<count<<endl;
      //if (count != 7) FATAL("theoretical uncertainties: line count != 7"); /// 7 when other unceratinties for ZH WH etc are ehre
      if (count != 5) FATAL("theoretical uncertainties: line count != 7");

   
      // CMS uncertainties.
      //
      // NOTE: naming conventions are from
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
      // TODO: updated numbers?

      datacard += "lumi_8TeV          lnN        1.026          1.026            -\n";
      //datacard += "eff_8TeV           lnN        1.100          1.100            -\n";
   
      
      // electron vs muon efficiency uncertainties, including uncertainties of
      // triggering efficiencies (it is assumed that the two are not correlated)
      if (TString(channel).Contains("ee"))
         datacard += "CMS_eff_e_8TeV     lnN        1.0001          1.0001          1.0001         -\n";
      else
         datacard += "CMS_eff_m_8TeV     lnN        1.0001          1.0001          1.0001         -\n";


      datacard += TString::Format("lepID    lnN       %-15f %-15f -\n",  sigSysonExp_lep[1][m], sigSysonExp_lep[0][m] );
      datacard += TString::Format("phoID    lnN       %-15f %-15f -\n",  sigSysonExp_pho[1][m], sigSysonExp_pho[0][m] );

      /*
      // photon efficiency uncertainties: ECAL barrel vs ECAL endcaps
      if (cat == 1 || cat == 2 || cat == 3)
         datacard += "CMS_eff_g_8TeV     lnN        1.006          1.006          1.006             -\n";
      else
         datacard += "CMS_eff_g_8TeV     lnN        1.010          1.010          1.010             -\n";

      // uncertainty on photon's R9 cut
      if (cat == 1 || cat == 2)
         datacard += "CMS_eff_R9_EB_8TeV lnN        1.050          1.050          1.050             -\n";

      float UEPS[] = {1.026, 1.035, 1.018, 1.021};
      datacard += TString::Format("CMS_UEPS_8TeV      lnN          -              -              -            %.3f          1.002            -\n", UEPS[cat-1]);

      float JEC[] = {1.028, 1.032, 1.022, 1.022};
      datacard += TString::Format("CMS_JEC_8TeV       lnN          -              -              -            %.3f          1.001            -\n", JEC[cat-1]);

      float JER[] = {1.010, 1.011, 1.011, 0.014};
      datacard += TString::Format("CMS_JER_8TeV       lnN          -              -              -            %.3f          1.001            -\n", JER[cat-1]);

      if (TString(channel).Contains("ee"))
         datacard += "pdf_PU_mu_8TeV     lnN        1.008          1.008          1.008          1.008          1.008            -\n";
      else
         datacard += "pdf_PU_mu_8TeV     lnN        1.004          1.004          1.004          1.004          1.004            -\n";

      datacard += "------------------------------------------------------------------------------------------------------------\n";
      */


      /*
      // declare uncertainties for the energy scale and resolution
      // (as described by the nuisance parameters above)
      // TODO: errors are set by hand
      datacard += TString::Format("CMS_scale_%s                                param          1     0.05\n", channel);
      datacard += TString::Format("CMS_res_%s                                  param          1     0.01\n", channel);
      */
      
      /*// declare parameters of the background PDF to be freely floating
      for (int i = 0; i < getNbkgPar(); i++) {
	TString nameStr = TString::Format("%s_bgr", getbkgfitParName(i).c_str());
         datacard += TString::Format("%-44s flatParam\n", nameStr.Data());
      }
      */


      cout<<"Writing data card finally"<<endl;
      // make datacard file
      FILE* out = fopen(TString::Format("output/datacards/datacard_hzg_%s_%s_8TeV_%d.txt",
                                        //channel, cats, mass[m]).Data(), "w");
					channel, cats, mass).Data(), "w");
      if (!out) FATAL("fopen() failed");
      if (fputs(datacard.Data(), out) < 0) FATAL("fputs() failed");
      if (fclose(out) != 0) FATAL("fclose() failed");

      cout<<"Wrote finally"<<endl;
   } // mass loop
   

   cout<<"=========================For cat========================="<<cat<<endl;
   for (int p = 0; p < 2; p++) {// 2 processes
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 

	 //cout<<sigeff_All[p][m * 5 + k]<<", ";
	 cout<<"mass is "<<mass<<" and eff is "<<sigeff_All[p][m * 5 + k]<<endl;
	 
       }//for (int k = 1; k <= 4; k++)
     }
     cout<<"}"<<endl;
   }

   ////////////////////////////////////////////////////////////////
   ////total ggF + VBF eff here///////
   cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 
	 cout<<sigeff_ggf_vbf[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;
     

   ///Draw the data and the fit to it with error bands
   int W = 800;
   int H = 600;
   
   int H_ref = 600;
   int W_ref = 800;
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;
   
   TCanvas* c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   cout<<"Fetching result now and plotting......"<<endl;
   //RooRealVar *mass = ws->var("CMS_hzg_mass");
   RooPlot *plot = fX->frame();
          
   cout<<"PLotting the data function now"<<endl;
   dataObs->plotOn(plot);

   cout<<"PLotting the bgrfit function now"<<endl;
   cout<<"address of bgrfit "<<bgrfit<<endl;


   dataObs->plotOn(plot);
   bgrfit->plotOn(plot);
   //bgrfit_ext->plotOn(plot,RooFit::LineColor(2));


   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   rfit->Print();

   cout<<"==================PRINTED r BEFORE PLOTTING===================== "<<endl;


   ///error bands
   ///1sigma
   /*bgrfit->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*rfit,2),RooFit::FillColor(kGreen-4));
   bgrfit->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*rfit,1),RooFit::FillColor(kYellow-4));
   */

   dataObs->plotOn(plot);
   bgrfit->plotOn(plot);
   //bgrfit_ext->plotOn(plot,RooFit::LineColor(2));
   
      
 
   char *outfilename = new char[100];
   char dirName[100] = "plots";
   

   string scat(cats);
   plot->Draw();
   sprintf(outfilename,"%s/%s.gif",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   sprintf(outfilename,"%s/%s.C",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   //GGF
   RooPlot *plotsig_ggf[3];
   int isig = 0;
   for(int imass=120; imass<=130; imass=imass+5){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     plotsig_ggf[isig] = fX->frame();
     cout<<"signal is "<<Form("signaldata_ggf_%d_cat%d",imass,cat) <<endl;
     rooMassSig[0][5*isig]->plotOn(plotsig_ggf[isig]); 
     cout<<"===========Printing the dataset information============="<<endl;
     rooMassSig[0][5*isig]->Print();
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[0][5*isig]->plotOn(plotsig_ggf[isig]);
     cout<<"got the signal"<<endl;
     ///save now
     
     TCanvas *c = new TCanvas("c","c",50,50,W,H);
     
     c->SetLeftMargin( L/W );
     c->SetRightMargin( R/W );
     c->SetTopMargin( T/H );
     c->SetBottomMargin( B/H );
     
     gStyle->SetOptStat(0);
     
     c->SetFillColor(0);
     c->SetBorderMode(0);
     c->SetBorderSize(2);
     c->SetTickx(1);
     c->SetTicky(1);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetLogy();

     //plot->SetMinimum(0);
     plotsig_ggf[isig]->Draw();
     c->Print(Form("plots/sig_ggf_mass%d_cat%d_%s.gif",imass,cat,channel));
     
     c->Print(Form("plots/sig_ggf_mass%d_cat%d.C",imass,cat));
     
     isig++;
   }


   ///VBF
   RooPlot *plotsig_vbf[3];
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+5){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     plotsig_vbf[isig] = fX->frame();
     cout<<"signal is "<<Form("signaldata_vbf_%d_cat%d",imass,cat) <<endl;
     rooMassSig[1][5*isig]->plotOn(plotsig_vbf[isig]); 
     cout<<"===========Printing the dataset information============="<<endl;
     rooMassSig[1][5*isig]->Print();
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[1][5*isig]->plotOn(plotsig_vbf[isig]);
     cout<<"got the signal"<<endl;
     ///save now
     
     TCanvas *c = new TCanvas("c","c",50,50,W,H);
     
     c->SetLeftMargin( L/W );
     c->SetRightMargin( R/W );
     c->SetTopMargin( T/H );
     c->SetBottomMargin( B/H );
     
     gStyle->SetOptStat(0);
     
     c->SetFillColor(0);
     c->SetBorderMode(0);
     c->SetBorderSize(2);
     c->SetTickx(1);
     c->SetTicky(1);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetLogy();

     //plot->SetMinimum(0);
     plotsig_vbf[isig]->Draw();
     c->Print(Form("plots/sig_vbf_mass%d_cat%d_%s.gif",imass,cat,channel));
     
     c->Print(Form("plots/sig_vbf_mass%d_cat%d.C",imass,cat));
     
     isig++;
   }

   
   ///all the signals on one


   int icol[] = {1,2,3,4,5,6,7,8,9,40,41};
   //ggF
   RooPlot *plotsig_ggf_all = fX->frame();;
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+1){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     
     cout<<"signal is "<<Form("signaldata_ggf_%d_cat%d",imass,cat) <<endl;
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[0][isig]->plotOn(plotsig_ggf_all,RooFit::LineColor(icol[isig]));
     cout<<"got the signal"<<endl;
     ///save now
     
     
     isig++;
   }

   c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetLogy();
   
   //plot->SetMinimum(0);
   plotsig_ggf_all->Draw();
   c->Print(Form("plots/sig_ggf_massAll_%s_cat%d.gif",channel,cat));
   c->Print(Form("plots/sig_ggf_massAll_%s_cat%d.C",channel,cat));



   ///VBF
   RooPlot *plotsig_vbf_all = fX->frame();;
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+1){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     
     cout<<"signal is "<<Form("signaldata_vbf_%d_cat%d",imass,cat) <<endl;
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[1][isig]->plotOn(plotsig_vbf_all,RooFit::LineColor(icol[isig]));
     cout<<"got the signal"<<endl;
     ///save now
     
     
     isig++;
   }

   c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetLogy();
   
   //plot->SetMinimum(0);
   plotsig_vbf_all->Draw();
   c->Print(Form("plots/sig_vbf_massAll_%s_cat%d.gif",channel,cat));
   c->Print(Form("plots/sig_vbf_massAll_%s_cat%d.C",channel,cat));
   
   
}

////end of simple mkdatacards for chk



//////////////////////////////////////BIAS STUDY///////////////////////////////////////////////
void statAn::biasStudymkdatacards(const char* channel, int cat)
{

  // Load the combine Library 
  //gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  
  string sigmodel = "RooCBxGaus";
  bool _nuisance = false;
  //bool _nuisance = true;
  
  //double lumi = 2.7*1000;///pb-1
  double lumi = mylumi; //pb-1
  // cat1, cat2, cat3 or cat4
  TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);

   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("hzg_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   
   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_8TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");


   TH1D hDataObs("hDataObs", "", 60,xmin, xmax);
   hDataObs.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");
   


   fX = new RooRealVar("x", "", xmin, xmax); 

   TString filename;

if(cat!=6789){
     if(schannel=="eeg") filename = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
     if(schannel=="mmg") filename = TString::Format("%s", "minitree_mu_data_doubleMu_SJ_out.root");
     dataObs = fill_events(&hDataObs, filename, 1, cat, false);
   }

   if(cat==6789){
     TString filename1 = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
     TString filename2 = TString::Format("%s", "minitree_mu_data_doubleMu_SJ_out.root");
     dataObs = fill_events(&hDataObs, filename1, filename2, 1, 1, cat, false);
   }
   
   dataObs->SetName( ("data_obs_"+catString) );
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass")))
     FATAL("RooWorkspace::import() failed");
   
   ///SJ
  
   // make background function



   fX->setBins(1000);
   //fX->setBins(1000,"fft");
   

   RooCategory cat_bias("pdf_index","Index of Pdf which is active");
   RooArgList mypdfs;

   //RooAbsPdf *multipdf;
   //poldeg = 4;
   //int totNdeg = 5;
   //int totNdeg = 2;

   model = "RooPolynomial";
   
   if(schannel=="eeg")
     {
       if(cat>=1 && cat<=3) poldeg = 3;
       
       if(cat>=4 && cat<=5) poldeg = 2;
       
       if(cat==10) poldeg = 1;
     }


   if(schannel=="mmg")
     {
       if(cat==1 || cat==4) poldeg = 3;
       if(cat==2) poldeg = 4;
       if(cat==3) poldeg = 2;
       if(cat==5 || cat==10) poldeg = 1;
     }
   

   if(schannel=="eeg_mmg"){

     poldeg = 1;
   }

   cout<<""<<endl;
   cout<<"====cat, channel, "<<cat<<" "<<schannel<<endl;
   cout<<"Mode and power "<<model<<" "<<poldeg<<endl;
   rfit = new RooFitResult();
   
   bgrfit = NULL;
   
   mybkgfit(xmin, xmax,cat);

     ////try 8th aug
     /*if(ideg==5)
       //bgrfit->forceNumIntegral(true);
       bgrfit->forceNumInt(true);
     */
     ////try 8th aug
	 
     /*//SJ
     cout<<"now setting bgrfit pdf to name of pdf "<<endl;
     //bgrfit->SetName(Form("pdf_%d",poldeg));
     bgrfit->SetName(Form("pdf"));
     mypdfs.add(*bgrfit);
	  
   
     
     //mypdfs.add(*bgrfit);

     cout<<"===DONE WITH THE FITTING"<<endl;
     bgrfit->Print();
     
     //////////////////////////////////////////////
     
     cout<<"in mkdatacards, address of r is "<<rfit<<endl;
     if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
     
     cout<<"Dumping RooFit result ..."<<endl;
     rfit->Dump();
     
     
     
     cout<<"address of bgrfit "<<bgrfit<<endl;
     
     cout<<"set the name to name of pdf "<<endl;
     if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		       RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		       RooFit::RenameVariable("x", "CMS_hzg_mass")))
       
       
       cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
     
     
     cout<<"now saved the bkg pdf"<<endl;
   
     
     //rfit->SetName(Form("pdf_bgr_%s_fitresult_%d",cats,poldeg));
     rfit->SetName(Form("pdf_bgr_%s_fitresult",cats));
     if ( wspace.import(*rfit) )
       FATAL("RooWorkspace::import() fit result failed");
     
     cout<<"Values of parameters after fixing ..."<<endl;
     //cout<<"Mean : sigma : stepval : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<stepval->getVal()<<endl;

     
     // total number of events observed/expected in the region [100, 190]
     int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
     //double expectation = bgrfit->getValV(&bkgpar_set);
     RooArgSet obs(*fX);
     //double expectation = bgrfit->getValV(&obs);
     
     
     RooAbsReal* intbkgfit1 = bgrfit->createIntegral(*fX) ; 
     //double expectation = intbkgfit1->getVal();
     //double expectation = bgrfit1->getVal();
     double expectation = fPar[0]->getVal();
     cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
     
     cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;

     cout<<"from histogram .... "<<hDataObs.Integral()<<endl;

     RooRealVar bkgNorm("norm","",expectation);
     //bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
     //bkgNorm.SetName(Form("pdf_%d_bgr_%s_norm",poldeg,cats));
     bkgNorm.SetName(Form("pdf_bgr_%s_norm",cats));
     if ( wspace.import(bkgNorm) )
       FATAL("RooWorkspace::import() failed");
     
     */
   
   
   int observation = 0;
   RooArgSet obs(*fX); 
   RooAbsReal* intbkgfit1;
   double expectation = 0;
   
   RooArgList bkgpar_list;
   RooRealVar bkgNorm("norm","",expectation);

   /// do it now for Exp

   
   cout<<"Doing Exp now for cat "<<cat<<endl;
   model = "Exp";
   poldeg = 1;

   if(schannel=="eeg")
     {
       poldeg = 1;
     }


   if(schannel=="mmg")
     {
       if(cat==1 || cat==3) poldeg = 3;
       else poldeg = 1;
     }

   if(schannel=="eeg_mmg")
     {
       poldeg = 1;
     }
   
   cout<<"Mode and power "<<model<<" "<<poldeg<<endl;

   rfit = new RooFitResult();
   bgrfit = NULL;
   mybkgfit(xmin, xmax,cat);
   
   cout<<"now setting bgrfit pdf to name of pdf "<<endl;
   bgrfit->SetName("pdf_Exp");
   mypdfs.add(*bgrfit);
   
   cout<<"===DONE WITH THE FITTING"<<endl;
   bgrfit->Print();
   
   
   cout<<"in mkdatacards, address of r is "<<rfit<<endl;
   if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
   
   cout<<"Dumping RooFit result ..."<<endl;
   rfit->Dump();
   
   
     
     cout<<"address of bgrfit "<<bgrfit<<endl;
     
     cout<<"set the name to name of pdf "<<endl;
     if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		       RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		       RooFit::RenameVariable("x", "CMS_hzg_mass")))
       cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
     
     
     cout<<"now saved the bkg pdf"<<endl;
     
     
     //rfit->SetName(("pdf_bgr_"+catString+"_fitresult"));
     rfit->SetName(Form("pdf_bgr_%s_fitresult_Exp",cats));
     if ( wspace.import(*rfit) )
       FATAL("RooWorkspace::import() fit result failed");
     
     
     
     cout<<"Values of parameters after fixing ..."<<endl;
        
     
     // total number of events observed/expected in the region [100, 190]
     observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
     
          
     intbkgfit1 = bgrfit->createIntegral(*fX) ; 
     //double expectation = intbkgfit1->getVal();
     //double expectation = bgrfit1->getVal();
     expectation = fPar[0]->getVal();
     cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
     
     cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;



     //RooRealVar bkgNorm("norm","",expectation);
     bkgNorm.SetName(Form("pdf_Exp_bgr_%s_norm",cats));
     if ( wspace.import(bkgNorm) )
       FATAL("RooWorkspace::import() failed");
      

   /// end of Exp

   

     ////Pow 
     cout<<"Doing Pow now "<<endl;
     model = "Pow";
     //poldeg = 1;
     if(schannel=="eeg")
       {
	 poldeg = 1;
       }
     
     
     if(schannel=="mmg")
       {
	 poldeg = 1;
       }
     
   if(schannel=="eeg_mmg")
     {
       poldeg = 1;
     }

     cout<<"Mode and power "<<model<<" "<<poldeg<<endl;
     rfit = new RooFitResult();
     bgrfit = NULL;
     mybkgfit(xmin, xmax,cat);
     
     
     cout<<"now setting bgrfit pdf to name of pdf "<<endl;
     bgrfit->SetName("pdf_Pow");
     mypdfs.add(*bgrfit);
     
     cout<<"===DONE WITH THE FITTING"<<endl;
     bgrfit->Print();
     
     
     cout<<"in mkdatacards, address of r is "<<rfit<<endl;
     if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
     
     cout<<"Dumping RooFit result ..."<<endl;
     rfit->Dump();
     
     cout<<"address of bgrfit "<<bgrfit<<endl;
     
     cout<<"set the name to name of pdf "<<endl;
     if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		       RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		       RooFit::RenameVariable("x", "CMS_hzg_mass")))
       
       
       cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
     
     
     cout<<"now saved the bkg pdf"<<endl;
     
     
     rfit->SetName(Form("pdf_bgr_%s_fitresult_Pow",cats));
     if ( wspace.import(*rfit) )
       FATAL("RooWorkspace::import() fit result failed");
     
     
     
     
     // total number of events observed/expected in the region [100, 190]
     observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
     //obs = *fX;
     
     intbkgfit1 = bgrfit->createIntegral(*fX) ; 
     //double expectation = intbkgfit1->getVal();
     //double expectation = bgrfit1->getVal();
     expectation = fPar[0]->getVal();
     cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
     
     cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;

     
     RooRealVar bkgNorm_gausPow("norm","",expectation);
     bkgNorm_gausPow.SetName(Form("pdf_Pow_bgr_%s_norm",cats));
     if ( wspace.import(bkgNorm_gausPow) )
       FATAL("RooWorkspace::import() failed");
     ///end of Pow


     ///Laurent
     cout<<"Doing Laurent now "<<endl;
     model = "Laurent";
     //poldeg = 1;

     if(schannel=="eeg")
       {
	 poldeg = 1;
	 if(cat==4) poldeg = 2;
	 
       }
     
     
     if(schannel=="mmg")
       {
	 poldeg = 1;
       }
     
   if(schannel=="eeg_mmg")
     {
       poldeg = 1;
     }
   
     cout<<"Mode and power "<<model<<" "<<poldeg<<endl;
     //if(cat==3) poldeg = 2;
     rfit = new RooFitResult();
     bgrfit = NULL;
     mybkgfit(xmin, xmax,cat);
     
     
     cout<<"now setting bgrfit pdf to name of pdf "<<endl;
     bgrfit->SetName("pdf_Laurent");
     mypdfs.add(*bgrfit);
     
     cout<<"===DONE WITH THE FITTING"<<endl;
     bgrfit->Print();
     
     
     cout<<"in mkdatacards, address of r is "<<rfit<<endl;
     if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
     
     cout<<"Dumping RooFit result ..."<<endl;
     rfit->Dump();
     
     cout<<"address of bgrfit "<<bgrfit<<endl;
     
     cout<<"set the name to name of pdf "<<endl;
     if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		       RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		       RooFit::RenameVariable("x", "CMS_hzg_mass")))
       
       
       cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
     
     
     cout<<"now saved the bkg pdf"<<endl;
     
     
     rfit->SetName(Form("pdf_bgr_%s_fitresult_Laurent",cats));
     if ( wspace.import(*rfit) )
       FATAL("RooWorkspace::import() fit result failed");
     
     
     
     
     // total number of events observed/expected in the region [100, 190]
     observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
     //obs = *fX;
     
     intbkgfit1 = bgrfit->createIntegral(*fX) ; 
     //double expectation = intbkgfit1->getVal();
     //double expectation = bgrfit1->getVal();
     expectation = fPar[0]->getVal();
     cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
     
     cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;

     
     RooRealVar bkgNorm_Laurent("norm","",expectation);
     bkgNorm_Laurent.SetName(Form("pdf_Laurent_bgr_%s_norm",cats));
     if ( wspace.import(bkgNorm_Laurent) )
       FATAL("RooWorkspace::import() failed");
     ////end of Laurent


     ///////start of adding other Bernstein Pols//////////////
     model = "RooPolynomial";
     
     int totdeg = 3;
     if(schannel=="eeg")
       {
	 if(cat==1) totdeg == 3;
	 
	 if(cat>=2 && cat<=4) totdeg = 2;
	 if(cat==5) totdeg = 1;
	 if(cat==6789) totdeg = 1;
       }
     
     
     if(schannel=="mmg")
       {
	 if(cat==1 || cat==2) totdeg = 4;
	 if(cat==3||cat==4) totdeg = 2;
	 if(cat==5) totdeg = 1;
	 if(cat==6789) totdeg = 1;
       }
     
     
     int inideg = 2;
     int findeg = 5;
     if(cat==5 || cat==6789 || cat==10){ 

       inideg = 1;
       findeg = 4;
     }

     cout<<"cat inideg findeg "<<cat<<" "<<inideg<<" "<<findeg<<endl;
     //for( int ideg=totdeg+1; ideg<=5; ideg++){
     for( int ideg=inideg; ideg<=findeg; ideg++){
       rfit = new RooFitResult();
       
       poldeg = ideg;
       
       bgrfit = NULL;
       
       mybkgfit(xmin, xmax,cat);
       //SJ

       cout<<"doine now Pol of order "<<ideg<<endl;
       cout<<"now setting bgrfit pdf to name of pdf "<<endl;
       //bgrfit->SetName(Form("pdf_%d",poldeg));
       bgrfit->SetName(Form("pdf_Bern%d",poldeg));
       mypdfs.add(*bgrfit);
       
       cout<<"===DONE WITH THE other Pols fitting FITTING"<<endl;
       bgrfit->Print();
       
       //////////////////////////////////////////////
       
       cout<<"in mkdatacards, address of r is "<<rfit<<endl;
       if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
       
       cout<<"Dumping RooFit result ..."<<endl;
       rfit->Dump();
       
       
       
       cout<<"address of bgrfit "<<bgrfit<<endl;
       
       cout<<"set the name to name of pdf "<<endl;
       if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
			 RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
			 RooFit::RenameVariable("x", "CMS_hzg_mass")))
	 
	 cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
       
       
       cout<<"now saved the bkg pdf"<<endl;
       
       
       //rfit->SetName(Form("pdf_bgr_%s_fitresult_%d",cats,poldeg));
       rfit->SetName(Form("pdf_bgr_%s_fitresult_Bern%d",cats,poldeg));
       if ( wspace.import(*rfit) )
	 FATAL("RooWorkspace::import() fit result failed");
       
       
       cout<<"Values of parameters after fixing ..."<<endl;
       
       
       // total number of events observed/expected in the region [100, 190]
       observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
       //double expectation = bgrfit->getValV(&bkgpar_set);
       intbkgfit1 = bgrfit->createIntegral(*fX) ; 
       //double expectation = intbkgfit1->getVal();
       //double expectation = bgrfit1->getVal();
       expectation = fPar[0]->getVal();
       cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
       
       cout<<"Bern order : from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<poldeg<<" "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;
       
       cout<<"from histogram .... "<<hDataObs.Integral()<<endl;
       
       bkgNorm.SetName(Form("pdf_Bern%d_bgr_%s_norm",poldeg,cats));
       if ( wspace.import(bkgNorm) )
	 FATAL("RooWorkspace::import() failed");
       
     }

     //////////////////////end od adding other Bernstein Pols///////////////
     
     
     /*
     //////Sech * Bernstein

   model = "SechBern";
   cout<<"doing SechBern now "<<endl;
   for(int ideg=2; ideg<=4; ideg++){
     
     poldeg = ideg;
     
     rfit = new RooFitResult();
     
     bgrfit = NULL;
     
     mybkgfit(xmin, xmax,cat);

     	 
     //SJ
     cout<<"now setting bgrfit pdf to name of pdf "<<endl;
     //bgrfit->GetPar(0)->setConstant(true);
     bgrfit->SetName(Form("pdf_SechBern%d",poldeg));
     mypdfs.add(*bgrfit);
	  

     cout<<"===DONE WITH THE FITTING"<<endl;
     bgrfit->Print();
   
     //////////////////////////////////////////////
     OB
     cout<<"in mkdatacards, address of r is "<<rfit<<endl;
     if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
     
     cout<<"Dumping RooFit result ..."<<endl;
     rfit->Dump();
     
     
     
     cout<<"address of bgrfit "<<bgrfit<<endl;
     
     cout<<"set the name to name of pdf "<<endl;
     if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		       RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		       RooFit::RenameVariable("x", "CMS_hzg_mass")))
       
       
       cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
     
     
     cout<<"now saved the bkg pdf"<<endl;
     
     
     //rfit->SetName(("pdf_bgr_"+catString+"_fitresult"));
     rfit->SetName(Form("pdf_bgr_%s_fitresult_SechBern%d",cats,poldeg));
     if ( wspace.import(*rfit) )
       FATAL("RooWorkspace::import() fit result failed");
     
     
     
     RooArgList bkgpar_list;
     bkgpar_list.add(*mean);
     bkgpar_list.add(*sigma);
     bkgpar_list.add(*stepval);
     bkgpar_list.add(*coeflist);
     RooArgSet bkgpar_set(bkgpar_list);
     
     
     cout<<"Values of parameters after fixing ..."<<endl;
     cout<<"Mean : sigma : stepval : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<stepval->getVal()<<endl;
   
     
     // total number of events observed/expected in the region [100, 190]
     int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
     //double expectation = bgrfit->getValV(&bkgpar_set);
     RooArgSet obs(*fX);
     //double expectation = bgrfit->getValV(&obs);
     
     
     RooAbsReal* intbkgfit1 = bgrfit->createIntegral(*fX) ; 
     //double expectation = intbkgfit1->getVal();
     //double expectation = bgrfit1->getVal();
     double expectation = fPar[0]->getVal();
     cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
     
     cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;



     RooRealVar bkgNorm("norm","",expectation);
     //bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
     bkgNorm.SetName(Form("pdf_SechBern%d_bgr_%s_norm",poldeg,cats));
     if ( wspace.import(bkgNorm) )
       FATAL("RooWorkspace::import() failed");
     
     
   }//for(int ideg=2; ideg<5; ideg++)
   //////end of sec * Bern
   */

   /*
   ////Sech * Pow
   cout<<"Doing SechPow now "<<endl;
   model = "SechPow";
   rfit = new RooFitResult();
   bgrfit = NULL;
   mybkgfit(xmin, xmax,cat);
   
   
   cout<<"now setting bgrfit pdf to name of pdf "<<endl;
   bgrfit->SetName("pdf_SechPow");
   mypdfs.add(*bgrfit);
   
   cout<<"===DONE WITH THE FITTING"<<endl;
   bgrfit->Print();
   
   
   cout<<"in mkdatacards, address of r is "<<rfit<<endl;
   if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;
   
   cout<<"Dumping RooFit result ..."<<endl;
   rfit->Dump();
   
   cout<<"address of bgrfit "<<bgrfit<<endl;
   
   cout<<"set the name to name of pdf "<<endl;
   if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		     RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		     RooFit::RenameVariable("x", "CMS_hzg_mass")))
     
     
     cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
   
   
   cout<<"now saved the bkg pdf"<<endl;
   
   
   rfit->SetName(Form("pdf_bgr_%s_fitresult_SechPow",cats));
   if ( wspace.import(*rfit) )
     FATAL("RooWorkspace::import() fit result failed");
   
   
   
   RooArgList bkgpar_list_sechPow;
   bkgpar_list_sechPow.add(*mean);
   bkgpar_list_sechPow.add(*sigma);
   bkgpar_list_sechPow.add(*alphabkg);
   bkgpar_list_sechPow.add(*betabkg);
   RooArgSet bkgpar_set_sechPow(bkgpar_list_sechPow);
   
   
   //cout<<"Values of parameters after fixing ..."<<endl;
   //cout<<"Mean : sigma : lambda : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<lambda->getVal()<<endl;
   
     
   // total number of events observed/expected in the region [100, 190]
   observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
   //obs = *fX;
   
   intbkgfit1 = bgrfit->createIntegral(*fX) ; 
   //double expectation = intbkgfit1->getVal();
   //double expectation = bgrfit1->getVal();
   expectation = fPar[0]->getVal();
   cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
   
   cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;
   
   
   
   RooRealVar bkgNorm_sechPow("norm","",expectation);
   bkgNorm_sechPow.SetName(Form("pdf_SechPow_bgr_%s_norm",cats));
   if ( wspace.import(bkgNorm_sechPow) )
     FATAL("RooWorkspace::import() failed");
   ///end of Sech * Pow
   */
   

   //////write the rootfile for bias study
   RooMultiPdf multipdf("roomultipdf","All Pdfs",cat_bias,mypdfs);
   // As usual make an extended term for the background with _norm for freely floating yield
   //RooRealVar norm("roomultipdf_norm","Number of background events",0,10000);
   RooRealVar norm("roomultipdf_norm","Number of background events",hDataObs.Integral(),0,10000);
   // Save to a new workspace
   
   //TFile *fout_bias = new TFile(Form("output/datacards/background_pdfs_%s.root",cats),"RECREATE");
   TFile *fout_bias = new TFile(Form("output/datacards/background_pdfs_%s_%s.root",cats,schannel.c_str()),"RECREATE");
   RooWorkspace wout("backgrounds","backgrounds");
   wout.import(cat_bias);
   wout.import(norm);
   //wout.import(multipdf);
   wout.import(multipdf, RooFit::RenameVariable("x", "CMS_hzg_mass"));
   wout.Print();
   wout.Write();
   fout_bias->Write();
   fout_bias->Close();
   //////END OF write the rootfile for bias study

   // save TH1 and TF1 objects into root file for later visualization
   wd->cd();
   // rename the X axis variable and add observed mass points into the workspace
   /*dataObs->SetName( ("data_obs_"+catString) );
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass")))
     FATAL("RooWorkspace::import() failed");
   */
   
   /////start with the signal now /////////
   
   cout<<"DOING SIGNAL FIT NOW"<<endl;
   double    expected[nprocess][nmass];
   //const char* proc[nprocess] = {"ggH", "VBF"};
   const char* proc[nprocess] = {"ggH", "VBF", "WplusH", "WminusH", "ZH", "ttH"};
   int mass_exist[nsig] = {120, 125, 130};
   double sigeff_All[nprocess][nmass];
   double sigNraw[nprocess][nmass];
   double sigeff_ggf_vbf[nmass];
   
   ///systematcs ID
   double sigSysonExp_lep[nprocess][nmass];
   double sigSysonExp_pho[nprocess][nmass];

   double sigSysonEff_lep[nprocess][nmass];
   double sigSysonEff_pho[nprocess][nmass];


   for (int p = 0; p < nprocess; p++){
     for (int m=0; m<nsig; m++) {
       
       cout<<"========================At teh very beginning, m is taking value now================"<<m<<endl;
       
       if((int)m >= (int)3){
	 cout<<"----inside >=3 loop, Will continue"<<endl;
	 continue;
	 cout<<"No, wait!!!! not continued"<<endl;
       }
       

       rooMassSig[p][m*5] = NULL;
       sigfit[p][m*5] = NULL;
       rsig[p][m*5] = NULL;
       
       
       int mass = mass_exist[m];
       
       cout<<"At the beginning, m, mass_exist is and m is "<<m<<" " <<mass_exist[m]<< " " <<mass<<endl;
       TString hname = TString::Format("hMass_%s_%d", proc[p], mass);
       
       TH1D hMass(hname, "", 120, xmin,xmax);
       hMass.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
       hMass.SetYTitle("Counts");
       hMass.Sumw2();
       
       
       TString filename;

       RooAbsData* rooMass;
       
       if(cat!=6789){
	 if(channel=="eeg") {
	   filename = TString::Format("minitree_ele_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 }
	 
	 
	 if(channel=="mmg") {
	   filename = TString::Format("minitree_mu_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 }
	 
	 
	 cout<<"proc and mass : Requesting to open "<<proc[p]<<" "<<mass<<" "<<filename<<endl;
	 
	 double trigEff = 1;
	 if( schannel == "eeg" ) trigEff = trigEff_ele;
	 if( schannel == "mmg" ) trigEff = trigEff_mu;
	 
	 cout<<"Calling for fill_events for isgn, trigEff is "<<trigEff<<endl;
	 rooMass = fill_events(&hMass, filename, trigEff, cat);
       }


       TString filename1;
       TString filename2;
       
       if(cat==6789){

	 filename1 = TString::Format("minitree_ele_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 filename2 = TString::Format("minitree_mu_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 
	 cout<<"proc and mass : Requesting to open "<<proc[p]<<" "<<mass<<" "<<filename1<<" "<<filename2<<endl;
	 
	 cout<<"Calling for fill_events for isgn, trigEff1 trigEff2 "<<trigEff_ele<<" "<<trigEff_mu<<endl;
	 rooMass = fill_events(&hMass, filename1, filename2, trigEff_ele, trigEff_mu, cat);
       

       }
       
       
       rooMassSig[p][5*m] = rooMass;
       
       cout<<"=================Channel, cat "<<channel<<" "<<cat<<endl;
       cout<<"===========fitting for mass =================="<<mass<<endl;
       mysigfit(sigmodel,xmin, xmax, _nuisance, p, 5*m, cat, schannel);
       
       setSigParConst(sigmodel, _nuisance);
     
       ///set the parameters const at this point for these mass points
       
    
       // save the normalization factor which adapts the MC scale to real
       // data; evaluate expected signal yield
       double norm = normfactor_SMHiggs_8TeV(channel, p, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double norm = 1;
       
     
       ///SJ
       double sigEff = hMass.Integral()/ngen(channel,p,mass);
       
       if(cat!=6789){
	 sigEff = hMass.Integral()/ngen(channel,p,mass);
       }

       if(cat==6789){
	 sigEff = hMass.Integral()/(ngen("eeg",p,mass)+ngen("mmg",p,mass)); 
       }

       //double sigExp = sigEff * lumi; /// use this for limit on xsec*BR
       double sigExp = hMass.Integral()*norm;///use this for limit on xsec/SMxsec
       
       ///systematics on ID
       double sys_lep = 0;
       double sys_pho = 0;
       
       double trigEff = 1;
       //getIDsys(filename,  trigEff, sys_lep, sys_pho, cat, true);
       //sigSysonExp_lep[p][m*5] = sqrt(sys_lep)*norm/sigExp;
       //sigSysonExp_pho[p][m*5] = sqrt(sys_pho)*norm/sigExp;
       

       sigSysonExp_lep[p][m*5] = sqrt(sys_lep)/hMass.Integral();
       sigSysonExp_pho[p][m*5] = sqrt(sys_pho)/hMass.Integral();
       
       ///sys on Eff is not always right. Only when sigExp=sigEff*lumi, then it is right. 
       ///Only the above sigSysonExp_lep numbers are right. These are on Eff   as well. So sigSysonExp_lep can be used on Exp and Eff
       sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen(channel,p,mass)*sigExp);
       sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen(channel,p,mass)*sigExp);
       
       
       cout<<"sqrt(lepsys) and sqrt(phosys) ; sigExp : "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<hMass.Integral()<<endl;
       
       cout<<" p : mass  : lepID sys EXP : EFF "<<p<<" "<<mass<<" "<<sigSysonExp_lep[p][m*5]<<" "<<sigSysonEff_lep[p][m*5] <<endl;
       cout<<" p : mass  : phoID sys EXP : EFF "<<p<<" "<<mass<<" "<<sigSysonExp_pho[p][m*5]<<" "<<sigSysonEff_pho[p][m*5] <<endl;
       

       expected[p][m*5] = sigExp;
       sigeff_All[p][m*5] = sigEff;
       
       sigNraw[p][m*5] = hMass.Integral();
       
       cout<<"Mass is, Signal eff "<<mass<<" " <<sigEff<<endl;

       cout<<"===p is, channel, cat, expectation=== "<<p<<" "<<channel<<" "<<cat <<" "<<expected[p][m*5]<<endl;
       ///SJ
       
       // unique suffix
       //TString sfx = TString::Format("sig_%s_%d", proc[p], mass[m]);
       TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);

       // - rename parameters of the PDF;
       // - declare parameters of the PDF to be constants (otherwise the
       //   parameters will be considered as freely floating, and the combine
       //   tool will produce weird results)
       

       /*for (int i = 0; i < sfit->GetNPar(); i++) {
	 const char* name = sfit->GetPar(i)->GetName();
	 sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
	 channel, cats, name, sfx.Data()));
	 sfit->GetPar(i)->setConstant(true);
         }
       */
       

       ///CHECK - set signal parameters to const
       
       rooMassSig[p][m*5]->SetName(TString::Format("signaldata_%s_%d_cat%d",proc[p], mass, cat));
       wspace.import(*rooMassSig[p][m*5],RooFit::RenameVariable("x", "CMS_hzg_mass"));
       
       
       cout<<"Inside the signal mass loop for filling the data set"<<endl;
       cout<<"=========Printing the dataset information========="<<endl;
       rooMassSig[p][m*5]->Print();
       cout<<"=========Printed the dataset information========="<<endl;

       // set names for the nuisance parameters of energy scale and resolution
       /*sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
       sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));
       */

       
         // - add signal PDF into the workspace under the name "pdf_sig_...";
         // - add suffix "_sig_..." to all subPDFs;
         // - connect the PDF's X axis variable to the X axis variable of data_obs
         // - connect the nuisance parameters of different signal PDFs together

       cout<<"setting the signal pdf name "<<endl;
       
       sigfit[p][5*m]->SetName("pdf");
       //sfit->SetName(Form("pdf_sig_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*sigfit[p][5*m], RooFit::RenameAllNodes(sfx),
			 RooFit::RenameVariable("x", "CMS_hzg_mass")))
	 FATAL("RooWorkspace::import() failed");
       
       
       cout<<"saved the signal pdf"<<endl;

       rsig[p][m*5]->SetName(Form("fitresult_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*rsig[p][m*5]))
	 FATAL("RooWorkspace::import() failed");
       
       cout<<"saved the signal fit result"<<endl;
       



       //delete sfit;
       
     } // process and mass loops
   }//for (int p = 0; p < nprocess; p++)

   cout<<"nprocess : nsig : "<<nprocess<<" "<<nsig<<endl;


   for (int m = 0; m < 3; m++) {
     
     int mass = mass_exist[m];   
     //double sigEff = (sigNraw[0][m * 5] + sigNraw[1][m * 5])/( ngen(channel,0,mass) + ngen(channel,1,mass) );   
     
     double ntotgen = 0;
     double numtot =  0;
     double systot = 0;
     
     for(int ip=0; ip<nprocess; ip++){
       
       double norm = normfactor_SMHiggs_8TeV(channel, ip, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       ntotgen += ngen(channel,ip,mass)*norm;
       numtot  += sigNraw[ip][m * 5]*norm;
       systot += pow(sigNraw[ip][m * 5]*norm,2); ///check
     }
     double sigEff = numtot/ntotgen;
     sigeff_ggf_vbf[m*5] = sigEff;

     //double sigTotsys_ggf_vbf = sqrt( pow(sigTotSys[0][m * 5],2) + pow(sigTotSys[1][m * 5],2) ); 
     //double sigTotsys_ggf_vbf = sqrt( systot );
     //sigsys_ggf_vbf[m*5] = sigTotsys_ggf_vbf;
   }

  
   ///interpolation
   // use 1 GeV steps
   for (int p = 0; p < nprocess; p++){ // 2 processes
   //for (int p = nprocess-1; p >=0; p--){ // 2 processes
     for (int mm = 0; mm< nsig-1; mm++){ 
       for (int k = 1; k <= 4; k++) {
	 
	 int mass = 120 + mm * 5 + k;
	 
	 cout<<""<<endl;
	 cout<<"p : mm : k : mass "<<p<<" "<<mm<<" "<<k<<" "<<mass<<endl;
	 
	 ///SJ changed
	 double m1 = mass_exist[mm];
	 double m2 = mass_exist[(mm + 1)];
	 
	 cout<<"m1, m2 and their eff are "<<m1<<" "<<m2 <<" Eff "<<sigeff_All[p][5*mm]<<" " <<sigeff_All[p][5*(mm+1)]<<endl;
	 // proportions of first/second neighbour
	 double a = (m2 - mass)/(m2 - m1);
	 double b = 1 - a;
	 
	 //RooRealVar* par1_fitresult = (RooRealVar*) fitresult->floatParsFinal()->find("par1"); 
	 
	 
	 // set values of parameters from linear extrapolation
	 ////y2-y/y2-y1 = x2-x/x2-x1; lets call RHS as a
	 ///gives y = y2(1-a) + ay1
	 //for (int i = 0; i < sfit->GetNPar() - 2; i++) {
	 //double val1 = sigfit[p][m]    ->GetTF1()->GetParameter(i);
	 //double val2 = sigfit[p][m + 1]->GetTF1()->GetParameter(i);
	 	    
	 
	 // evaluate expected signal yield
	 double norm = normfactor_SMHiggs_8TeV(channel, p, mass);
	 //expected[p][m * 5 + k] = norm * sfit->GetPar(0)->getVal();
	 
	 ///extrapolate the eff - SJ
	 expected[p][mm * 5 + k] = expected[p][5*mm] * a + b * expected[p][5*(mm+1)];
	 
	 //cout<<"Setting the element ("<<p<<","<<(m * 5 + k)<<")"
	 //sigeff_All[p][mm * 5 + k] = expected[p][mm * 5 + k]/lumi;
	 sigeff_All[p][mm * 5 + k] = sigeff_All[p][mm*5] * a + b * sigeff_All[p][5*(mm+1)];

	 sigeff_ggf_vbf[mm * 5 + k] = sigeff_ggf_vbf[mm*5] * a + b * sigeff_ggf_vbf[5*(mm+1)];

	 //systematics
	 sigSysonExp_lep[p][mm * 5 + k] = sigSysonExp_lep[p][mm*5] * a + b * sigSysonExp_lep[p][5*(mm+1)];
	 sigSysonExp_pho[p][mm * 5 + k] = sigSysonExp_pho[p][mm*5] * a + b * sigSysonExp_pho[p][5*(mm+1)];

	 sigSysonEff_lep[p][mm * 5 + k] = sigSysonEff_lep[p][mm*5] * a + b * sigSysonEff_lep[p][5*(mm+1)];
	 sigSysonEff_pho[p][mm * 5 + k] = sigSysonEff_pho[p][mm*5] * a + b * sigSysonEff_pho[p][5*(mm+1)];

	 


	 cout<<"========starting to interpolate now for mass "<<mass<<"================="<<endl;
	 cout<<"====masses taken are "<<m1<<" and "<<m2<<endl;
	 siginterpolate(sigmodel, p, mm*5, k, a,b, _nuisance, cat);
	 sigfit[p][5*mm+k] = sigfit1; //filled in mysigfunc
	 
	 setSigParConst(sigmodel, _nuisance);

	 cout<<"========END OF interpolate ================="<<endl;

	 // - add signal PDF into the workspace under the name "pdf_sig_...";
	 // - add suffix "_sig_..." to all subPDFs;
	 // - connect the PDF's X axis variable to the X axis variable of data_obs
	 // - connect the nuisance parameters of different signal PDFs together
	 
	 TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);
	 sigfit[p][5*mm+k]->SetName("pdf");
	 if (wspace.import(*sigfit[p][5*mm+k], RooFit::RenameAllNodes(sfx),
			   RooFit::RenameVariable("x", "CMS_hzg_mass")))
	   FATAL("RooWorkspace::import() failed");
	 

	 /*
	 rsig[p][mm*5+k]->SetName(Form("fitresult_%s_%d_cat%d",proc[p], mass, cat));
	 if (wspace.import(*rsig[p][mm*5+k]))
	   FATAL("RooWorkspace::import() failed");
	 */

	 // unique suffix
	 
	 
	 // - rename parameters of the PDF;
	 // - declare parameters of the PDF to be constants (otherwise the
	 //   parameters will be considered as freely floating, and the combine
	 //   tool will produce weird results)
	 //for (int i = 0; i < sfit->GetNPar(); i++) {
	 //const char* name = sfit->GetPar(i)->GetName();
	 //sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
	 //channel, cats, name, sfx.Data()));
	 //sfit->GetPar(i)->setConstant(true);
	 //}
	 
	 // set names for the nuisance parameters of energy scale and resolution
	 //sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
	 //sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));


	 //deleteSigPar(_nuisance);
	 
       } // extrapolation
     }//for (int mm = 0; mm< nsig; mm++)
   }//for (int p = 0; p < nprocess; p++)


            


   

  //////end of the signal now//////////////





   
   

   cout<<"Writing hte work space now"<<endl;
   wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_8TeV.root",
                                      channel, cats));
   
   




   // close output/for_visualizations.root
   delete fo;
   

      cout<<"WRITING CARD NOW "<<endl;
   //
   // produce datacards
   //
   
   // cuttent UTC time
   TString timestamp(TTimeStamp().AsString());
   timestamp.Remove(timestamp.Last(':'));

   // name of final state
   TString binString = TString::Format("%s_%s_8TeV", channel, cats);
   const char* bin = binString.Data();

   //for (int m = 0; m < 31; m++) {
   for (int m = 0; m < 11; m++) {
      int mass = 120 + m;
      


      sigSysonExp_lep[0][m] = 1 + sigSysonExp_lep[0][m];
      sigSysonExp_pho[0][m] = 1 + sigSysonExp_pho[0][m];

      sigSysonExp_lep[1][m] = 1 + sigSysonExp_lep[1][m];
      sigSysonExp_pho[1][m] = 1 + sigSysonExp_pho[1][m];
      
      cout<<"mass is "<<mass<<endl;
      TString datacard;
      datacard += "# Datacard for the HZg analysis for limit setting\n";
      datacard += "# National Central University, Taiwan\n";
      datacard += "# " + timestamp + " (UTC)\n";
      datacard += TString::Format("# Usage: combine -U -M Asymptotic -m %d datacard.txt\n", mass);
      datacard += "#\n";
      datacard += "imax 1  # number of final states\n";
      datacard += "jmax *  # number of yields given below minus one\n";
      datacard += "kmax *  # number of sources of systematical uncertainties (nuisance parameters)\n";
      datacard += "------------------------------------------------------------------------------------------------------------\n";
     

      datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      
      
      datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ZH_%d_%s\n", channel, cats, mass,cats);

      /*
      datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ZH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_WminusH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_WplusH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ttH_%d_%s\n", channel, cats, mass,cats);
      */
   
      /*
      //21 jan, 2017 - use ggF shape for other processes for cats 1-4 since these events yield is quite small
      if(cat>=1&&cat<=4){
	datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      }


      //21 jan, 2017 - use VBF shape for other processes for cats 5 since these events yield is quite small
      if(cat==5){
	datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      }
      */
      
      
      //datacard += TString::Format("shapes  bgr       * background_pdfs_%s.root backgrounds:roomultipdf\n", cats);
      datacard += TString::Format("shapes  bgr       * background_pdfs_%s_%s.root backgrounds:roomultipdf\n", cats,schannel.c_str());

      datacard += TString::Format("shapes  data_obs  * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:data_obs_%s\n", channel, cats,cats);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %s\n", bin);
      datacard += TString::Format("bin            %s\n", cats);
      //datacard += TString::Format("observation    %d\n", observation);
      datacard += TString::Format("observation    %d\n", -1);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      
      /*
      datacard += TString::Format("bin            %-15s %-15s %-15s\n", cats, cats, cats);
      datacard +=                 "process        VBFH            ggH             bgr\n";
      datacard +=                 "process        -1              0               1\n";
      datacard += TString::Format("rate           %-15f %-15f %-15f\n",
				  expected[1][m], expected[0][m], 1.0);
      
      */


      //datacard += TString::Format("bin            %-15s %-15s %-15s\n", cats, cats, cats);
      datacard += TString::Format("bin            %-15s  %-15s\n", cats, cats);
      if( (cat>=1 && cat<=4) || cat==10 ){ datacard +=                 "process            ggH             bgr\n";
      datacard +=                 "process               0               1\n";
      datacard += TString::Format("rate           %-15f %-15f\n",
				  expected[0][m]+expected[1][m]+expected[2][m]+expected[3][m]+expected[4][m]+expected[5][m], 1.0);
      }


      if(cat==5){ datacard +=                 "process            VBFH             bgr\n";
      datacard +=                 "process               0               1\n";
      datacard += TString::Format("rate           %-15f %-15f\n",
				  //expected[1][m], 1.0);
				  expected[0][m]+expected[1][m]+expected[2][m]+expected[3][m]+expected[4][m]+expected[5][m], 1.0);
      }




      if(cat==6789){ datacard +=                 "process            ZH             bgr\n";
      datacard +=                 "process               0               1\n";
      datacard += TString::Format("rate           %-15f %-15f\n",
				  //expected[1][m], 1.0);
				  expected[0][m]+expected[1][m]+expected[2][m]+expected[3][m]+expected[4][m]+expected[5][m], 1.0);
      }



    


      /*datacard += TString::Format("bin            %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n", cats, cats, cats,cats, cats,cats, cats);
      datacard +=                 "process        ttH       ZH      WminusH        WplusH           VBFH            ggH             bgr\n";
      datacard +=                 "process        -5        -4       -3             -2                -1             0               1\n";

      cout<<"all the expectations ... "<<expected[0][m]<<" "<<expected[1][m]<<" "<<expected[2][m]<<" "<<expected[3][m]<<" "<<expected[4][m]<<" "<<expected[5][m]<<endl;
      
      datacard += TString::Format("rate           %-15f %-15f %-15f %-15f %-15f %-15f %-15f\n",
				  expected[5][m], expected[4][m], expected[3][m], expected[2][m], expected[1][m], expected[0][m], 1.0);

      datacard += "------------------------------------------------------------------------------------------------------------\n";
      */


      ifstream in("theoretical_uncertainties_SM_Higgs.list");
      int count = 0; // simple error protection
      char line[10000];

      /*
      while (in.good()) {
	in.getline(line, 10000);
	if (!in.eof() && in.fail())
	   FATAL("ifstream::getline() failed");

         TString s = line;
	 cout<<""<<endl;
	 cout<<"line is "<<line<<endl;
	 
	 cout<<TString::Format("%.2f ", (double)mass)<<endl;
         //if(s.BeginsWith(TString::Format("%.1f ", (double)mass))) {
	 if(s.BeginsWith(TString::Format("%.2f ", (double)mass))) {  ///needs to be changed when we have intervals of 0.5 GeV
	   cout<<"added this line "<<endl;
            s.Remove(0, s.First(':') + 1); // remove e.g. "155.0   :"
            datacard += s + "\n";
            count++;
         }
      }
      cout<<"count is "<<count<<endl;
      //if (count != 5) FATAL("theoretical uncertainties: line count != 7"); ///when only ggH and VBF are there
      if (count != 8) FATAL("theoretical uncertainties: line count != 8"); /// 8 when other unceratinties for ZH WH etc are ehre
      */

   
      // CMS uncertainties.
      //
      // NOTE: naming conventions are from
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
      // TODO: updated numbers?

      //      datacard += "lumi_8TeV          lnN        1.026          1.026            -\n";
      //datacard += "lumi_8TeV          lnN     1.062          1.062   1.062    1.062   1.062          1.062            -\n";
      //datacard += "eff_8TeV           lnN        1.100          1.100            -\n";
   
      datacard += "lumi_8TeV          lnN                1.026            -\n";
      /*
      // electron vs muon efficiency uncertainties, including uncertainties of
      // triggering efficiencies (it is assumed that the two are not correlated)
      if (TString(channel).Contains("ee"))
         datacard += "CMS_eff_e_8TeV     lnN        1.0001          1.0001          1.0001         -\n";
      else
         datacard += "CMS_eff_m_8TeV     lnN        1.0001          1.0001          1.0001         -\n";


      datacard += TString::Format("lepID    lnN       %-15f %-15f -\n",  sigSysonExp_lep[1][m], sigSysonExp_lep[0][m] );
      datacard += TString::Format("phoID    lnN       %-15f %-15f -\n",  sigSysonExp_pho[1][m], sigSysonExp_pho[0][m] );

      */
      ///we advertise the pdf changing parameter as a discrete nuisance
      datacard += "pdf_index discrete \n";

      cout<<"Writing data card finally"<<endl;
      // make datacard file
      FILE* out = fopen(TString::Format("output/datacards/datacard_hzg_%s_%s_8TeV_%d.txt",
                                        //channel, cats, mass[m]).Data(), "w");
					channel, cats, mass).Data(), "w");
      if (!out) FATAL("fopen() failed");
      if (fputs(datacard.Data(), out) < 0) FATAL("fputs() failed");
      if (fclose(out) != 0) FATAL("fclose() failed");

      cout<<"Wrote finally"<<endl;
   } // mass loop
   

   cout<<"=========================For cat========================="<<cat<<endl;
   for (int p = 0; p < 2; p++) {// 2 processes
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 

	 //cout<<sigeff_All[p][m * 5 + k]<<", ";
	 cout<<"mass is "<<mass<<" and eff is "<<sigeff_All[p][m * 5 + k]<<endl;
	 
       }//for (int k = 1; k <= 4; k++)
     }
     cout<<"}"<<endl;
   }

   ////////////////////////////////////////////////////////////////
   ////total ggF + VBF eff here///////
   cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 
	 cout<<sigeff_ggf_vbf[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;
     

   ///Draw the data and the fit to it with error bands
   int W = 800;
   int H = 600;
   
   int H_ref = 600;
   int W_ref = 800;
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;
   
   TCanvas* c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   cout<<"Fetching result now and plotting......"<<endl;
   //RooRealVar *mass = ws->var("CMS_hzg_mass");
   RooPlot *plot = fX->frame();
          
   cout<<"PLotting the data function now"<<endl;
   dataObs->plotOn(plot);

   cout<<"PLotting the bgrfit function now"<<endl;
   cout<<"address of bgrfit "<<bgrfit<<endl;


   dataObs->plotOn(plot);
   bgrfit->plotOn(plot);
   //bgrfit_ext->plotOn(plot,RooFit::LineColor(2));


   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   rfit->Print();

   cout<<"==================PRINTED r BEFORE PLOTTING===================== "<<endl;


   ///error bands
   ///1sigma
   /*bgrfit->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*rfit,2),RooFit::FillColor(kGreen-4));
   bgrfit->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*rfit,1),RooFit::FillColor(kYellow-4));
   */

   dataObs->plotOn(plot);
   bgrfit->plotOn(plot);
   //bgrfit_ext->plotOn(plot,RooFit::LineColor(2));
   
      
 
   char *outfilename = new char[100];
   char dirName[100] = "plots";
   

   string scat(cats);
   plot->Draw();
   sprintf(outfilename,"%s/%s.gif",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   sprintf(outfilename,"%s/%s.C",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   //GGF
   RooPlot *plotsig_ggf[3];
   int isig = 0;
   for(int imass=120; imass<=130; imass=imass+5){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     plotsig_ggf[isig] = fX->frame();
     cout<<"signal is "<<Form("signaldata_ggf_%d_cat%d",imass,cat) <<endl;
     rooMassSig[0][5*isig]->plotOn(plotsig_ggf[isig]); 
     cout<<"===========Printing the dataset information============="<<endl;
     rooMassSig[0][5*isig]->Print();
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[0][5*isig]->plotOn(plotsig_ggf[isig]);
     cout<<"got the signal"<<endl;
     ///save now
     
     TCanvas *c = new TCanvas("c","c",50,50,W,H);
     
     c->SetLeftMargin( L/W );
     c->SetRightMargin( R/W );
     c->SetTopMargin( T/H );
     c->SetBottomMargin( B/H );
     
     gStyle->SetOptStat(0);
     
     c->SetFillColor(0);
     c->SetBorderMode(0);
     c->SetBorderSize(2);
     c->SetTickx(1);
     c->SetTicky(1);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetLogy();

     //plot->SetMinimum(0);
     plotsig_ggf[isig]->Draw();
     c->Print(Form("plots/sig_ggf_mass%d_cat%d_%s.gif",imass,cat,channel));
     
     c->Print(Form("plots/sig_ggf_mass%d_cat%d.C",imass,cat));
     
     isig++;
   }


   ///VBF
   RooPlot *plotsig_vbf[3];
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+5){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     plotsig_vbf[isig] = fX->frame();
     cout<<"signal is "<<Form("signaldata_vbf_%d_cat%d",imass,cat) <<endl;
     rooMassSig[1][5*isig]->plotOn(plotsig_vbf[isig]); 
     cout<<"===========Printing the dataset information============="<<endl;
     rooMassSig[1][5*isig]->Print();
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[1][5*isig]->plotOn(plotsig_vbf[isig]);
     cout<<"got the signal"<<endl;
     ///save now
     
     TCanvas *c = new TCanvas("c","c",50,50,W,H);
     
     c->SetLeftMargin( L/W );
     c->SetRightMargin( R/W );
     c->SetTopMargin( T/H );
     c->SetBottomMargin( B/H );
     
     gStyle->SetOptStat(0);
     
     c->SetFillColor(0);
     c->SetBorderMode(0);
     c->SetBorderSize(2);
     c->SetTickx(1);
     c->SetTicky(1);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetLogy();

     //plot->SetMinimum(0);
     plotsig_vbf[isig]->Draw();
     c->Print(Form("plots/sig_vbf_mass%d_cat%d_%s.gif",imass,cat,channel));
     
     c->Print(Form("plots/sig_vbf_mass%d_cat%d.C",imass,cat));
     
     isig++;
   }

   
   ///all the signals on one


   int icol[] = {1,2,3,4,5,6,7,8,9,40,41};
   //ggF
   RooPlot *plotsig_ggf_all = fX->frame();;
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+1){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     
     cout<<"signal is "<<Form("signaldata_ggf_%d_cat%d",imass,cat) <<endl;
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[0][isig]->plotOn(plotsig_ggf_all,RooFit::LineColor(icol[isig]));
     cout<<"got the signal"<<endl;
     ///save now
     
     
     isig++;
   }

   c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetLogy();
   
   //plot->SetMinimum(0);
   plotsig_ggf_all->Draw();
   c->Print(Form("plots/sig_ggf_massAll_%s_cat%d.gif",channel,cat));
   c->Print(Form("plots/sig_ggf_massAll_%s_cat%d.C",channel,cat));



   ///VBF
   RooPlot *plotsig_vbf_all = fX->frame();;
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+1){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     
     cout<<"signal is "<<Form("signaldata_vbf_%d_cat%d",imass,cat) <<endl;
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[1][isig]->plotOn(plotsig_vbf_all,RooFit::LineColor(icol[isig]));
     cout<<"got the signal"<<endl;
     ///save now
     
     
     isig++;
   }

   c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetLogy();
   
   //plot->SetMinimum(0);
   plotsig_vbf_all->Draw();
   c->Print(Form("plots/sig_vbf_massAll_%s_cat%d.gif",channel,cat));
   c->Print(Form("plots/sig_vbf_massAll_%s_cat%d.C",channel,cat));
   
   
}

////end of simple mkdatacards for chk











//////////////////////////////////////////////END OF BIAS STUDY///////////////////////////////////////////


////small copy of mkdatacard
//______________________________________________________________________________
void statAn::smallcopymkdatacard(const char* channel, int cat)
{
   /* Does everything for particular channel and particular event category.
    *
    * channel = either of "eeg" (for electrons) or "mmg" (for muons);
    * cat = 1, 2, 3 or 4.
    *
    * NOTE: we follow naming conventions of
    * https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
    *
    * NOTE: in this function, memory is leaked (via "new") for coding brevity.
    */

  //double lumi = 2.7*1000;///pb-1
  double lumi = 2689.8; ///pb-1

   // cat1, cat2, cat3 or cat4
   TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);

   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("hzg_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   //TFile* fo = new TFile("output/for_visualizations.root", "RECREATE");
   /*if (!fo || fo->IsZombie())
      FATAL("TFile::Open() failed");
   */
   
   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_8TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");

   //
   // evaluate PDF of background from real data
   //

   TH1D hDataObs("hDataObs", "", 60, 110, 170);
   hDataObs.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");

   // false = do not add per-event weights
   //TString filename = TString::Format("%s", "input/output_job_2photon_2012abcd_Jan22rereco.root");
   TString filename;
   if(schannel=="eeg") filename = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
   if(schannel=="mmg") filename = TString::Format("%s", "minitree_mu_data_doubleEG_SJ_out.root");
   RooAbsData* dataObs = fill_events(&hDataObs, filename, 1, cat, false);

   //TF1RooGaussStepBernstein bgrfit(110, 170, 5); // 5 = polynomial degree
   TF1RooPolynomial bgrfit(110, 170, 4);

   ///SJ
   //TF1RooPolynomial bgrfit(110, 170, 5);
   //TF1RooPolynomial bgrfit(110, 170, 3);
   //TF1RooPolynomial bgrfit(110, 170, 4);
   //TF1RooGaussStepBernstein bgrfit(110, 170, 4); // 5 = polynomial degree
   //TF1RooGaussStepBernstein bgrfit(110, 170, 5); // 5 = polynomial degree

   //TF1RooPolynomial bgrfit(110, 170, 5);
   RooFitResult *r = bgrfit.FitTo(&hDataObs, dataObs);
   cout<<"Dumping RooFit result ..."<<endl;
   r->Dump();
   // save TH1 and TF1 objects into root file for later visualization
   wd->cd();
   hDataObs.Write("", TObject::kOverwrite);
   bgrfit.GetTF1()->Write("fit_hDataObs", TObject::kOverwrite);
   
   // rename the X axis variable and add observed mass points into the workspace
   dataObs->SetName( ("data_obs_"+catString) );
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass")))
      FATAL("RooWorkspace::import() failed");

   // - rename parameters of the background PDF;
   // - declare parameters of the background PDF to be constants (otherwise the
   //   parameters will be considered as freely floating even without using the
   //   "flatParam" directives below)
   for (int i = 0; i < bgrfit.GetNPar(); i++) {
      const char* name = bgrfit.GetPar(i)->GetName();
      bgrfit.GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s",
                                                channel, cats, name));
      // bgrfit.GetPar(i)->setConstant(true);
   }

   // - make normalization utilized by the combine tool to be equal to expected
   //   number of events written into datacards;
   // - add (extended version of) background PDF into the workspace under the
   //   name "pdf_bgr";
   // - add suffix "_bgr" to all PDF's parameters (and to subPDFs, if any);
   // - connect the PDF's X axis variable to the X axis variable of data_obs
   //bgrfit.GetPar(0)->setVal(1);
   bgrfit.GetExtPdf()->SetName("pdf");

   /*   if (wspace.import(*bgrfit.GetExtPdf(), RooFit::RenameAllNodes("bgr"),
                                          RooFit::RenameAllVariablesExcept("bgr", "x"),
                                          RooFit::RenameVariable("x", "CMS_hzg_mass")))
   */
   //if (wspace.import(*bgrfit.GetExtPdf(), RooFit::RenameAllNodes("bgr"),

   //SJ
   bgrfit.GetPdf()->SetName("pdf");
   if (wspace.import(*bgrfit.GetPdf(), RooFit::RenameAllNodes(("bgr_"+catString) ),
		     RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		     RooFit::RenameVariable("x", "CMS_hzg_mass")))

      FATAL("RooWorkspace::import() failed");

   // total number of events observed/expected in the region [100, 190]
   int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
   double expectation = bgrfit.GetTF1()->GetParameter(0);

   cout<<"EXPECTATION FROM FITTING IS "<<expectation<<endl;

   ///SJ
   RooRealVar bkgNorm("norm","",expectation);
   bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
   if ( wspace.import(bkgNorm) )
     FATAL("RooWorkspace::import() failed");


   r->SetName(("pdf_bgr_"+catString+"_fitresult"));
   if ( wspace.import(*r) )
     FATAL("RooWorkspace::import() fit result failed");

   //RooDouble(expectation).Write("bgr_pdf_norm", TObject::kOverwrite);

   // important consistency check (inconsistency may happen e.g. if by a mistake
   // the mass region in mkminitrees.cc is different from [100, 190]
   if (observation != dataObs->sumEntries())
      FATAL("observation != dataObs->sumEntries()");

   //
   // evaluate signal PDFs for existing MC productions
   //
   
   // - references to signal PDFs for existing MC productions
   // - expected signal yields for all mass points

   wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_8TeV.root",
                                      channel, cats));

   // close output/for_visualizations.root
   delete fo;



   ///Draw the data and the fit to it with error bands
   int W = 800;
   int H = 600;
   
   int H_ref = 600;
   int W_ref = 800;
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;
   
   TCanvas* c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   TFile *fout = TFile::Open(TString::Format("output/datacards/for_datacards_hzg_%s_%s_8TeV.root",channel, cats));
   
   if(fout==NULL) cout<<"Root file is null "<<endl;
   
   RooWorkspace *ws = (RooWorkspace *)fout->Get("hzg_workspace"); 
   //RooRealVar *mass = ws->var("CMS_hzg_mass");
   RooRealVar x("x", "", 125, 110, 170);
   //RooRealVar x("x", "", 125, 130, 170);
   RooPlot *plot = x.frame();
   //RooAbsData *data = ws->data(Form("data_obs_cat%s",cats)); 

   //RooAbsPdf *extbkg = ws->pdf(Form("pdf_bgr_cat%s",cats));

   if(ws==NULL)
     cout<<"WARNING!!!ws "<<endl;

   if(&x==NULL)
     cout<<"WARNING!!!mass "<<endl;

   if(plot==NULL)
     cout<<"WARNING!!!plot "<<endl;

   if(dataObs==NULL)
     cout<<"WARNING!!!data "<<endl;

   if(bgrfit.GetPdf()==NULL)
     cout<<"WARNING!!!extbkg "<<endl;
       
   dataObs->plotOn(plot);
   bgrfit.GetPdf()->plotOn(plot);
   bgrfit.GetExtPdf()->plotOn(plot,RooFit::LineColor(kRed));


   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   r->Print();

   ///error bands
   ///1sigma
   /*bgrfit.GetPdf()->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*r,1),RooFit::FillColor(kOrange-6));
   bgrfit.GetPdf()->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*r,2),RooFit::FillColor(kOrange-10));
   */

   bgrfit.GetExtPdf()->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*r,1),RooFit::FillColor(kYellow-4));
   bgrfit.GetExtPdf()->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*r,2),RooFit::FillColor(kGreen-4));
      
   char *outfilename = new char[100];
   char dirName[100] = "plots";
   

   string scat(cats);
   plot->Draw();
   sprintf(outfilename,"%s/%s.gif",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   sprintf(outfilename,"%s/%s.C",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
}

///end of small copy of mkdatacards



//______________________________________________________________________________
void mkdatacards()
{
   /* Steering function.
    */

   // FIXME: this line is required on lxplus and pcncu1X machines
   gSystem->SetIncludePath("-I$ROOFITSYS/include");

   gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
   gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");

   // suppress useless messages from RooFit
//    RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   ///SJ commented it
   //RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   // prevent segfaults due to various TH1 objects not accurately attributed to
   // correct TDirectory
   TH1::AddDirectory(0);

   // prepare empty directory for datacards
   if (gSystem->Exec("rm -fr output/datacards") != 0 || gSystem->mkdir("output/datacards") != 0)
      FATAL("failed to remove or create directory \"./output/datacards\"");

   // remove previous version of for_visualizations.root, if any
   if (!gSystem->AccessPathName("output/for_visualizations.root"))
      if (gSystem->Unlink("output/for_visualizations.root") != 0)
         FATAL("TSystem::Unlink() failed");



   int _poldeg = 3;
   statAn s("RooGaussStepBernstein", _poldeg); ;
   
   // actual work
   //for (int cat = 1; cat <=4; cat++) {
   for (int cat = 1; cat <=5; cat++) {
   //for (int cat = 1; cat <=1; cat++) {
     //mkdatacard("eeg", cat);
     //mkdatacard("mmg", cat);
     cout<<""<<endl;
     cout<<"Electron cat "<<cat<<endl;
     s.simplemkdatacard("eeg", cat);
     cout<<""<<endl;
     cout<<"Muon cat "<<cat<<endl;
     s.simplemkdatacard("mmg", cat);
     
   }
   

   ///if inclusive category which is either 1 till 5
   cout<<"=========Electron============"<<endl;
   //mkdatacard("eeg", 0);
   s.simplemkdatacard("eeg", 0);
   cout<<"=========Muon============"<<endl;
   s.simplemkdatacard("mmg", 0);
   //mkdatacard("mmg", 0);
   

   const char* proc[2] = {"ggH", "VBF"};
   int mass_exist[3] = {120, 125, 130};
   
   double sigeff_All[10][100];
   double sigeff_All_mu[10][100];

   for (int p = 0; p < 2; p++){
     //cout<<"{";
     for (int m = 0; m < 3; m++) {
       int mass = mass_exist[m];
       
       TString filename = TString::Format("minitree_ele_sig%s_%d_out.root",
					  proc[p], mass);
       sigeff_All[p][m] = getSigEff(filename,true)/ngen("eeg",p,mass);
       
       filename = TString::Format("minitree_mu_sig%s_%d_out.root",
					  proc[p], mass);
       sigeff_All[p][m] = getSigEff(filename,true)/ngen("mmg",p,mass);
       
       
     }
     cout<<"}"<<endl;
   }

   for (int p = 0; p < 2; p++){
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       cout<<sigeff_All[p][m]<<", ";
     }
     cout<<"}"<<endl;
   }
   

}



/////////////call bias studyn function here////
void biasStudy()
{
   /* Steering function.
    */

   // FIXME: this line is required on lxplus and pcncu1X machines
   gSystem->SetIncludePath("-I$ROOFITSYS/include");

   gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
   gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");

   // suppress useless messages from RooFit
//    RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   ///SJ commented it
   //RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   // prevent segfaults due to various TH1 objects not accurately attributed to
   // correct TDirectory
   TH1::AddDirectory(0);

   // prepare empty directory for datacards
   if (gSystem->Exec("rm -fr output/datacards") != 0 || gSystem->mkdir("output/datacards") != 0)
     FATAL("failed to remove or create directory \"./output/datacards\"");
   
   // remove previous version of for_visualizations.root, if any
   if (!gSystem->AccessPathName("output/for_visualizations.root"))
     if (gSystem->Unlink("output/for_visualizations.root") != 0)
       FATAL("TSystem::Unlink() failed");
   
   
   
   int _poldeg = 3;
   statAn s("RooGaussStepBernstein", _poldeg); ;
   
   // actual work
   //for (int cat = 1; cat <=4; cat++) {
   for (int cat = 1; cat <=5; cat++) {
   //for (int cat = 1; cat <=1; cat++) {
     cout<<""<<endl;
     cout<<"Electron cat "<<cat<<endl;
     s.biasStudymkdatacards("eeg", cat);
     cout<<""<<endl;
     cout<<"Muon cat "<<cat<<endl;
     s.biasStudymkdatacards("mmg", cat);
     
   }
   
   s.biasStudymkdatacards("eeg_mmg", 6789);

   cout<<"===doing for cat 10==="<<endl;
   s.biasStudymkdatacards("eeg", 10);
   s.biasStudymkdatacards("mmg", 10);



   /*
   ///if inclusive category which is either 1 till 5
   cout<<"=========Electron============"<<endl;
   s.biasStudymkdatacards("eeg", 0);
   cout<<"=========Muon============"<<endl;
   s.biasStudymkdatacards("mmg", 0);
   //mkdatacard("mmg", 0);
   */

   
}


/////////////////////end of call to bias function//////////////////////////

///get expectation within some range of the Mlg mass
///////
///scat tells if we want untagged or the tagged one 
///scat is either "untagged" OR "all"

void getExpectedEvents(int cat, string scat, string channel, double xmin, double xmax, double sigma_away, bool useRange) ///if useRange is set to true, then the events are calculated within that range else within the sig_eff
{
  //cat 0 is the combined category: 1,2,3,4,5

  std::cout<<"Inside getExpectedEvents; cat is "<<cat <<std::endl;

  CalcSigma t;

  double Estbkg = 0;
  double data = 0;
  double sig125_ggf = 0;
  double sig125_vbf = 0;

 
  double sigN[5] = {0,0,0,0,0};
  double sigTotN = 0;
  string infilename;

  double dataN[5] = {0,0,0,0,0};
  double dataTotN = 0;

  if(channel == "eeg") infilename = "files_ele.list";
  if(channel == "mmg") {
    infilename = "files_mu.list";
  }


  string infile = "";

  double sigma_ggF[5];
  double sigma_VBF[5];

  //ggF
  if(channel=="eeg")
    infile = Form("minitree_ele_sigggH_%d_out.root",125);

  else if(channel=="mmg")
    infile = Form("minitree_mu_sigggH_%d_out.root",125);
  
  else cout<<"Channel WRONGLY given, give either eeg or mmg"<<endl;

  for(int ii=1; ii<=4; ii++)
    {
      double xmin_s = 0;
      double xmax_s = 0;
      t.Loop(infile,xmin, xmax, ii, xmin_s, xmax_s);
      sigma_ggF[ii-1] = (xmax_s - xmin_s)/2.;
    }
  
  //VBF
  if(channel=="eeg")
    infile = Form("minitree_ele_sigVBF_%d_out.root",125);

  else if(channel=="mmg")
    infile = Form("minitree_mu_sigVBF_%d_out.root",125);
  
  else cout<<"Channel WRONGLY given, give either eeg or mmg"<<endl;

  for(int ii=1; ii<=4; ii++)
    {
      double xmin_s = 0;
      double xmax_s = 0;

      t.Loop(infile,xmin, xmax, ii, xmin_s, xmax_s);
      sigma_VBF[ii-1] = (xmax_s - xmin_s)/2.;
    }

  cout<<"For cat "<<cat<<" and mass 125 GeV, sigma_eff for ggF and VBf "<<sigma_ggF<<" "<<sigma_VBF<<endl;
  
  cout<<"using the same sigma for ggF and VBF ... ggF"<<endl;

  double low[5];
  double high[5];
  
  for(int ii=1; ii<=4; ii++){
    low[ii-1] = (125-sigma_away*sigma_ggF[ii-1]);
    high[ii-1] = (125+sigma_away*sigma_ggF[ii-1]);
  }

  ///for now - VBF category - take the same as the untagged sigma 
  low[4] = low[3];
  high[4] = high[3];
  

  cout<<" low an high "<<low <<" "<<high<<endl;
  
   // open file and get requested minitree
  ifstream file;
  file.open(infilename.c_str(), ifstream::in );
  char filename[200];
  TChain *tree = (TChain *)new TChain("minitree");
  while(!file.eof()){
    file >> filename;
    cout<<"filename = "<<filename<<endl;
    if(strncmp(filename,"#",1)==0)
      {
	continue;
      }
    tree->Add(filename);
  }//while(!datafile.eof())
                           
  if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float mcwei;
   Int_t itype;
   
   double br = 3.36/100.;
   double xsec_ggf[] = {47.38*1.1e-3*br,44.14*1.533e-3*br,41.23*1.941e-3*br};
   double xsec_vbf[] = {3.935*1.1e-3*br, 3.782*1.533e-3*br,3.637*1.941e-3*br};

   //double lumi = 2569.10 *(1+0.023);
   double lumi = 2689.8; //pb-1

   // associate tree branches with variables
   /*
   if (tree->SetBranchAddress("category", &category) != 0)
     FATAL("TTree::SetBranchAddress() category failed");
   
   if (tree->SetBranchAddress("itype", &itype) != 0)
     FATAL("TTree::SetBranchAddress() itype failed");
   
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
     FATAL("TTree::SetBranchAddress() hzg_mass failed");
   
   if (tree->SetBranchAddress("mcwei", &mcwei) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() mcwei failed");
   */

   tree->SetBranchAddress("category", &category);
   tree->SetBranchAddress("itype", &itype);
   tree->SetBranchAddress("hzg_mass", &hzg_mass);
   tree->SetBranchAddress("mcwei", &mcwei);
    
   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if (cat > 0 && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=5) && scat=="all")
       continue;

     if (cat==0 && !(category>=1 && category<=4) && scat=="untagged")
       continue;


     if (cat==5 && !(category==5))
       continue;
     
          
     if(useRange){
       if(hzg_mass<xmin || hzg_mass>xmax) continue;
     }
     
     else{
       if(hzg_mass<low[category-1] || hzg_mass>high[category-1]){
	 cout<<"cat is : low : high : hzg "<<category<<" "<<low[category-1]<<" "<<high[category-1]<<" "<<hzg_mass<<endl;
	 continue;  
       }
     }
     

     if(hzg_mass<low[category-1] || hzg_mass>high[category-1]) {

       cout<<"CHECK THE CODE ... Range not rejected "<<endl;
     }
     
     
     //cout<<"came here "<<endl;
 
     if(itype==100 || itype==101) ///100 for electron and 101 for double muon
       {
	 data += 1;
	 dataN[category-1]++;
	 dataTotN++;
       }


     double tmp_xsec_vbf = 1;
     if(itype==2) tmp_xsec_vbf = xsec_vbf[1];

     double tmp_xsec_ggf = 1;
     if(itype==12) tmp_xsec_ggf = xsec_ggf[1];


     double trigEff = 1;
     if( channel == "eeg" ) trigEff = trigEff_ele;
     if( channel == "mmg" ) trigEff = trigEff_mu;
     
     
     if(itype<0) Estbkg += mcwei*lumi*trigEff;
     if(itype==2) {
       //cout<<"xsec VBF "<<tmp_xsec_vbf<<endl;
       //sig125_vbf += mcwei*tmp_xsec_vbf*lumi/ngen(channel.c_str(),1,125); //vbF
       sig125_vbf += tmp_xsec_vbf*lumi/ngen(channel.c_str(),1,125); //vbF

       sigN[category-1] = sigN[category-1]+tmp_xsec_vbf*lumi/ngen(channel.c_str(),1,125);
       sigTotN += tmp_xsec_vbf*lumi/ngen(channel.c_str(),1,125);
     }
     
     if(itype==12){
       //cout<<"xsec ggF "<<tmp_xsec_ggf<<endl;
       //sig125_ggf += mcwei*tmp_xsec_ggf*lumi/ngen(channel.c_str(),0,125); //ggF
       sig125_ggf += tmp_xsec_ggf*lumi/ngen(channel.c_str(),0,125); //ggF

       sigN[category-1] = sigN[category-1]+tmp_xsec_ggf*lumi/ngen(channel.c_str(),0,125);
       sigTotN += tmp_xsec_ggf*lumi/ngen(channel.c_str(),0,125);
     }


   }

   cout<<"Data observed \t SM expectation \t Signal expectation (ggF) \t Signal expectaion(VBF) \t"<<data<<" \t "<<Estbkg<<" \t "<<sig125_ggf<<" \t "<<sig125_vbf<<endl;

   if ( (cat==0 && scat=="all") || (cat==0 && scat=="untagged") )
     for(int icat=1; icat<=5; icat++){
       cout<<""<<endl;
       double perc = sigN[icat-1]/sigTotN;
       cout<<"SIGNAL : frac for icat is "<<icat <<" "<<perc<<endl;


       perc = dataN[icat-1]/dataTotN;
       cout<<"DATA : frac for icat is "<<icat <<" "<<perc<<endl;
     }

}



double ngen(const char* channel, int p, double mass){

  /* p = production process number (0=ggH, 1=qqH, 2=WH, 3=ZH, 4=ttH);*/
  /*
  cout<<"Inside ngen"<<endl;
  cout<<"Mass is "<<mass<<endl;
  cout<<"channel is "<<channel<<endl;
  */

  if(mass!=120 && mass!=125 && mass!=130) {
    cout<<"Give either 120, 125 or 130 as mass!! exiting "<<endl;
    FATAL("wrong mass given");
  }


  ///120, 125, 130
  double nev_eeg[6][3] = {
    /*
    { 33463,33267, 66968},
    {33315,33075,32695}
    */

    {33463+66519, 33267+129713, 33116+66968}, //ggF
    {33315, 33075+66792, 32524}, //VBF
    {3289, 10014, 3173}, //W+H
    {3221, 10039, 3272}, //W-H
    {6675, 18918, 6408}, //ZH
    //{4462, 8675, 4387} //ttH
    {4462, 8675, 3925} //ttH

  };

  double nev_mmg[6][3] = {

    /*
    {33246,33492,66435 },
    {33398,32695,32948}
    */

    {33245+66353, 33492+129925, 32934+66434}, //ggF
    {33399, 32695+66482, 32948}, //VBF
    {3376, 9962, 3201}, //W+H
    {3237, 9912, 3361}, //W-H
    {6554, 18856, 6338}, //ZH
    //{4500, 8621, 4456} //ttH
    {4500, 8621, 3980} //ttH


};


  int index = (mass - 120)/5;

  string schan(channel);
  
  if(schan=="eeg"){
    //cout<<"events for eeg, Mass "<<mass<<" and channel "<<p<<"are "<<nev_eeg[p][index]<<endl;
    return nev_eeg[p][index];
  }
  else if(schan=="mmg"){
    //cout<<"events for mmg, Mass "<<mass<<" and channel "<<p<<"are "<<nev_mmg[p][index]<<endl;
    return nev_mmg[p][index];  
  }
  
  //else  FATAL("wrong channel given");   
  else return nev_eeg[p][index] +  nev_mmg[p][index];

}


