/* Performs unbinned fitting of invariant mass shapes and creates text datacards
 1;2c* for the "combine" tool from the HiggsAnalysis/CombinedLimit package.
 *
 1;2c1;2c* Minitrees are taken from output/minitrees.root. The following branches are
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

//====================================VISUALIZE ERROR===================
//https://root.cern.ch/root/html/tutorials/roofit/rf610_visualerror.C.html
 // Visualize 1-sigma error encoded in fit result 'r' as orange band using linear error propagation
  // This results in an error band that is by construction symmetric
  //
  // The linear error is calculated as
  // error(x) = Z* F_a(x) * Corr(a,a') F_a'(x)
  //
  // where     F_a(x) = [ f(x,a+da) - f(x,a-da) ] / 2, 
  // 
  //         with f(x) = the plotted curve 
  //              'da' = error taken from the fit result
  //        Corr(a,a') = the correlation matrix from the fit result
  //                Z = requested significance 'Z sigma band'
  //
  // The linear method is fast (required 2*N evaluations of the curve, where N is the number of parameters), 
  // but may not be accurate in the presence of strong correlations (~>0.9) and at Z>2 due to linear and 
  // Gaussian approximations made
  //

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

#include "RooAbsArg.h"

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

#include <RooGaussModel.h>
#include <RooTruthModel.h>
#include <RooDecay.h>

#include "RooAddModel.h"

#include <RooNumConvPdf.h>
#include <RooGenericPdf.h>

#include <RooExponential.h>
#include <RooFFTConvPdf.h>

#include <fitting_functions/RooStepBernstein.h>

#include "RooHist.h"
#include "CMS_lumi.C"
#include "TLegend.h"

#include "TLegendEntry.h"

#include "TMathText.h"

#include "RooCurve.h"

#include "RooMinimizer.h"

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

//bool verbose_ = true;
bool verbose_ = false;


//double trigEff_ele = 0.98*0.98*0.997; //leg1 * leg2 * DZ
//double trigEff_mu  = 0.95*0.95*0.97; //leg1 * leg2 * DZ


//double trigEff_ele = 1.0; //leg1 * leg2 * DZ
//double trigEff_mu  = 1.0; //leg1 * leg2 * DZ

//AN2016-408
//double trigEff_ele = 0.966*0.966*0.99; //leg1 * leg2 * DZ
//double trigEff_mu  = 0.937*0.954*0.98; //leg1 * leg2 * DZ

//double trigEff_ele = 0.99*0.99; //leg1 * leg2 * DZ
//double trigEff_mu  = 0.976*0.986; //leg1 * leg2 * DZ

//double trigEff_ele = 0.997686; //dZ
double trigEff_ele = 1; //dZ
double trigEff_mu  = 1; //dZ


int nBinsForMass = 55;


///fit range
//double xmin = 100; 
double xmin = 115; 
double xmax = 170;

double blind_min = 120;
double blind_max = 130;


//double xmin = 100; 
//double xmax = 170;

//double mylumi =12900; //pb-1
double mylumi =36460; //pb-1

double ngen(const char* channel = "eeg", int p=0, double mass=125);
//void getExpectedEvents(int cat, string channel, double xmin, double xmax, double sigma_away=2,bool useRange=false);
void getExpectedEvents(int cat, string scat, string channel, double xmin, double xmax, double sigma_away=2,bool useRange=false); ///if useRange is set to true, then the events are calculated within that range else within the sig_eff
double getSigEff(const char* filename, bool usewei = true);



TCanvas* setTCanvasNicev1(string name){
  int W = 800;
  int H = 600;
  
  int H_ref = 600;
  int W_ref = 800;
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;
  
  TCanvas* c = new TCanvas(name.c_str(),name.c_str(),50,50,W,H);
  
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
  return c;
}

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
     if(hzg_mass<110 || hzg_mass>170) continue;
      
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
  cout<<"file name is "<<filename<<endl;
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
   //RooRealVar x("x", "", 125, 100, 170);
   //RooRealVar x("x", "", 125, xmin, xmax);
   /////do it 
   //x.setBins(500);
   RooRealVar w("w", "", 1);  // NOTE: weights are ignored without this variable
   RooDataSet* rooMass;
   if (usewei)
      rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX, w), RooFit::WeightVar(w));
   else
      rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX));

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if ( (cat > 0 && cat!=6789) && category != cat)
       continue;
     
     //if (cat==0 && !(category>=1 && category<=5) )
     //if (cat==0 && !(category>=1 && category<=9) ) ///include the eff of lepton tagged events when i qupte the total eff
     if (cat==0 && !(category>=1 && category<=10) ) ///include the eff of lepton tagged events when i qupte the total eff
       continue;


     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;
     
     *fX = hzg_mass;
     
     
     ///SJ
     if(hzg_mass<xmin || hzg_mass>xmax) continue;
     
     if (usewei) {
	hMass->Fill(hzg_mass, puwei*trigEff);
	//rooMass->add(RooArgSet(x), puwei*trigEff);
	rooMass->add(RooArgSet(*fX),puwei);
	//rooMass->add(RooArgSet(x),1);
	//cout<<"For signal trigEff = "<<trigEff<<endl;
	//hMass->Fill(hzg_mass, trigEff);
	//rooMass->add(RooArgSet(x), trigEff);
	
     } else {
       hMass->Fill(hzg_mass);
       rooMass->add(RooArgSet(*fX));
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
   //RooRealVar x("x", "", 125, xmin, xmax);
   RooRealVar w("w", "", 1);  // NOTE: weights are ignored without this variable
   RooDataSet* rooMass;
   if (usewei)
      rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX, w), RooFit::WeightVar(w));
   else
      rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX));

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
     //if (cat==0 && !(category>=1 && category<=9) ) ///include the eff of lepton tagged events when i qupte the total eff
     if (cat==0 && !(category>=1 && category<=10) ) ///include the eff of lepton tagged events when i qupte the total eff
       continue;


     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;
     
     *fX = hzg_mass;
     

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

       hMass->Fill(hzg_mass, puwei*trigEff);
       
       if(puwei>0){
	 //rooMass->add(RooArgSet(x), puwei*trigEff);
	 rooMass->add(RooArgSet(*fX),puwei);
	 
       }//if(puwei>0)
     } else {
       hMass->Fill(hzg_mass);
	rooMass->add(RooArgSet(*fX));
     }
   }
   
   delete tree;
   //delete fi;

   return rooMass;
}

///////end of fill_events for LT///////////////////


/////////////////////////////////fill events for LT wplus + wminus - 4 rootfiles to combine//////////////////////

RooAbsData* statAn::fill_events(TH1* hMass, const char* filename1, const char* filename2,const char* filename3,  const char* filename4 , double trigEff1, double trigEff2, int cat, bool usewei)
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
  tree->Add(filename3);
  tree->Add(filename4);
  
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
   //RooRealVar x("x", "", 125, xmin, xmax);
   RooRealVar w("w", "", 1);  // NOTE: weights are ignored without this variable
   RooDataSet* rooMass;
   if (usewei)
      rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX, w), RooFit::WeightVar(w));
   else
      rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX));

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
     //if (cat==0 && !(category>=1 && category<=9) ) ///include the eff of lepton tagged events when i qupte the total eff
     if (cat==0 && !(category>=1 && category<=10) ) ///include the eff of lepton tagged events when i qupte the total eff
       continue;


     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;
     
     *fX = hzg_mass;
     

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

       hMass->Fill(hzg_mass, puwei*trigEff);
       
       if(puwei>0){
	 //rooMass->add(RooArgSet(x), puwei*trigEff);
	 rooMass->add(RooArgSet(*fX),puwei);
	 
       }//if(puwei>0)
     } else {
       hMass->Fill(hzg_mass);
	rooMass->add(RooArgSet(*fX));
     }
   }
   
   delete tree;
   //delete fi;

   return rooMass;
}


///////////////////////////////// END of fill events for LT wplus + wminus - 4 rootfiles to combine//////////////////////


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

   double nev = 0;
   double nev_raw = 0;
   

   double tmpsys_lep = 0;
   double tmpsys_pho = 0;
   
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
     
     if (cat==0 && !(category>=1 && category<=10) )
       continue;
     
     if(hzg_mass<xmin || hzg_mass>xmax) continue;

     nev += puwei*trigEff;

     nev_raw++;

     //cout<<"idSys_lep idSys_pho "<<idSys_lep << " "<<idSys_pho<<endl;
     //sys_lep += idSys_lep*idSys_lep*trigEff*trigEff;
     //sys_pho += idSys_pho*idSys_pho*trigEff*trigEff;


     sys_lep += idSys_lep*trigEff;
     sys_pho += idSys_pho*trigEff;
     
     tmpsys_lep += idSys_lep*trigEff;
     tmpsys_pho += idSys_pho*trigEff;

     /*
     cout<<"nev_raw nev , sqrt(sys_lep) and sqrt(sys_pho) till now "<<nev_raw<<" "<<nev<<" "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<endl;
     cout<<"Just adding gives "<<tmpsys_lep<<" "<<tmpsys_pho<<endl;
     cout<<"====Till now .... % Ele - quad : just add "<<sqrt(sys_lep)/nev<<" "<<tmpsys_lep/nev<<endl;
     cout<<"====Till now .... % Pho - quad : just add "<<sqrt(sys_pho)/nev<<" "<<tmpsys_pho/nev<<endl;
     */

     //sys_lep += idSys_lep*idSys_lep;
     //sys_pho += idSys_pho*idSys_pho;


   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   ///bec in the end it is sqrt
   sys_lep = pow(sys_lep,2);
   sys_pho = pow(sys_pho,2);
       
   //cout<<"final sys_lep sys_pho "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<endl; 
   
}/////end of ID sys


////////////////////PU systematics

/////////////for ID sys
void statAn::getPUsys(const char* filename, double trigEff , double &sys_pu, int cat, bool usewei)
{
  std::cout<<"Inside getIDsys; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  sys_pu = 0;

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
   float puwei_up;
   float puwei_down;
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


   if (usewei && tree->SetBranchAddress("mcwei_puup", &puwei_up) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");


   if (usewei && tree->SetBranchAddress("mcwei_pudown", &puwei_down) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");
   
   if (tree->SetBranchAddress("idSys_lep", &idSys_lep) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");

   if (tree->SetBranchAddress("idSys_pho", &idSys_pho) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");

   double sys_up = 0;
   double sys_down = 0;

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if (cat > 0 && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=10) )
       continue;
     

     if(hzg_mass<xmin || hzg_mass>xmax) continue;
     
     //sys_up += sqrt( pow(puwei_up*trigEff,2) );
     //sys_down += sqrt( pow(puwei_down*trigEff,2) );

     sys_up += puwei_up*trigEff;
     sys_down += puwei_down*trigEff;


     //sys_up += ( pow(puwei_up*trigEff,2) );
     //sys_down += ( pow(puwei_down*trigEff,2) );
     
   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   //sys_pu = sqrt(TMath::Max( fabs(sys_up), fabs(sys_down) ));

   sys_pu = (TMath::Max( fabs(sys_up), fabs(sys_down) ));
   

}/////end of PUsys



/////////////////////PU systematics


/////////////////////ID and PU sys for LT lepton tagged//////////////////////////

void statAn::getIDsys(const char* filename1, const char* filename2, double trigEff1 , double trigEff2, double &sys_lep, double &sys_pho, int cat, bool usewei)
{
  std::cout<<"Inside getIDsys; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  sys_lep = 0;
  sys_pho = 0;

  cout<<"trigEff1 trigEff2 "<<trigEff1<<" "<<trigEff2<<endl;
   // open file and get requested minitree

  TChain *tree = new TChain("minitree");
  tree->Add(filename1);
  tree->Add(filename2);
  
  if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;
   float idSys_lep;
   float idSys_pho;

   double nev = 0;
   double nev_raw = 0;
   

   double tmpsys_lep = 0;
   double tmpsys_pho = 0;
   
   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }
   
   if (tree->SetBranchAddress("idSys_lep", &idSys_lep) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   if (tree->SetBranchAddress("idSys_pho", &idSys_pho) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if ( (cat > 0  && cat!=6789) && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=10) )
       continue;

     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;


     if(hzg_mass<xmin || hzg_mass>xmax) continue;
     
     double trigEff = 1;
     string fname(tree->GetCurrentFile()->GetName());
     
     int fstr;
     fstr = fname.find("mu",0);
     if(fstr!=string::npos)
       trigEff = trigEff2;
     
     fstr = fname.find("ele",0);
     if(fstr!=string::npos)
       trigEff = trigEff1;



     nev += puwei*trigEff;
     
     nev_raw++;

     //cout<<"trigEff "<<trigEff<<endl;
     //cout<<"idSys_lep idSys_pho "<<idSys_lep << " "<<idSys_pho<<endl;
     //sys_lep += idSys_lep*idSys_lep*trigEff*trigEff;
     //sys_pho += idSys_pho*idSys_pho*trigEff*trigEff;


     sys_lep += idSys_lep*trigEff;
     sys_pho += idSys_pho*trigEff;
     
     tmpsys_lep += idSys_lep*trigEff;
     tmpsys_pho += idSys_pho*trigEff;

     /*
     cout<<"nev_raw nev , sqrt(sys_lep) and sqrt(sys_pho) till now "<<nev_raw<<" "<<nev<<" "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<endl;
     cout<<"Just adding gives "<<tmpsys_lep<<" "<<tmpsys_pho<<endl;
     cout<<"====Till now .... % Ele - quad : just add "<<sqrt(sys_lep)/nev<<" "<<tmpsys_lep/nev<<endl;
     cout<<"====Till now .... % Pho - quad : just add "<<sqrt(sys_pho)/nev<<" "<<tmpsys_pho/nev<<endl;
     */

     //sys_lep += idSys_lep*idSys_lep;
     //sys_pho += idSys_pho*idSys_pho;


   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   ///bec in the end it is sqrt
   sys_lep = pow(sys_lep,2);
   sys_pho = pow(sys_pho,2);
       
   //cout<<"final sys_lep sys_pho "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<endl; 
   
}/////end of ID sys


////////////////////PU systematics

/////////////for ID sys
void statAn::getPUsys(const char* filename1, const char* filename2, double trigEff1 , double trigEff2, double &sys_pu, int cat, bool usewei)
{
  std::cout<<"Inside getIDsys; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  sys_pu = 0;

  cout<<"trigEff are "<<trigEff1<<" "<<trigEff2<<endl;
   // open file and get requested minitree
  TChain *tree = new TChain("minitree");
  tree->Add(filename1);
  tree->Add(filename2);
  
  if (!tree) FATAL("TFile::Get() failed");


   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;
   float puwei_up;
   float puwei_down;
   float idSys_lep;
   float idSys_pho;

   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
     {
       // FATAL("TTree::SetBranchAddress() failed");
     }
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }


   if (usewei && tree->SetBranchAddress("mcwei_puup", &puwei_up) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }


   if (usewei && tree->SetBranchAddress("mcwei_pudown", &puwei_down) != 0) ///SJ
     {
       // FATAL("TTree::SetBranchAddress() failed");
     }
   
   if (tree->SetBranchAddress("idSys_lep", &idSys_lep) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   if (tree->SetBranchAddress("idSys_pho", &idSys_pho) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   double sys_up = 0;
   double sys_down = 0;

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if ( (cat > 0  && cat!=6789) && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=10) )
       continue;
     
     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;
     

     if(hzg_mass<xmin || hzg_mass>xmax) continue;

     double trigEff = 1;
     string fname(tree->GetCurrentFile()->GetName());
     
     int fstr;
     fstr = fname.find("mu",0);
     if(fstr!=string::npos)
       trigEff = trigEff2;
     
     fstr = fname.find("ele",0);
     if(fstr!=string::npos)
       trigEff = trigEff1;




     
     
     //sys_up += sqrt( pow(puwei_up*trigEff,2) );
     //sys_down += sqrt( pow(puwei_down*trigEff,2) );

     sys_up += puwei_up*trigEff;
     sys_down += puwei_down*trigEff;

     //sys_up += ( pow(puwei_up*trigEff,2) );
     //sys_down += ( pow(puwei_down*trigEff,2) );
     
   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   //sys_pu = sqrt(TMath::Max( fabs(sys_up), fabs(sys_down) ));

   sys_pu = (TMath::Max( fabs(sys_up), fabs(sys_down) ));
   

}/////end of PUsys


///////////////////end of ID and PU sys for LT lepton tagged////////////////////

////////////////////////////HLT systematics/////////////////////
void statAn::getHLTsys(const char* filename, double trigEff , double &sys_lep, int cat, bool usewei)
{
  std::cout<<"Inside getHLTsys; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  sys_lep = 0;


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
   float hltSys_lep;


   double nev = 0;
   double nev_raw = 0;
   

   double tmpsys_lep = 0;

   
   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
      FATAL("TTree::SetBranchAddress() failed");
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
      FATAL("TTree::SetBranchAddress() failed");

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");
   
   if (tree->SetBranchAddress("hltSys_lep", &hltSys_lep) != 0) ///SJ
      FATAL("TTree::SetBranchAddress() failed");

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if (cat > 0 && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=10) )
       continue;
     
     if(hzg_mass<xmin || hzg_mass>xmax) continue;

     nev += puwei*trigEff;

     nev_raw++;
     
     string fname(tree->GetCurrentFile()->GetName());
     int fstr = fname.find("mu",0);
     if(fstr!=string::npos)
       {
	 hltSys_lep = 1.3/100;
       }


     sys_lep += hltSys_lep*trigEff;

     
     tmpsys_lep += hltSys_lep*trigEff;



   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   ///bec in the end it is sqrt
   sys_lep = pow(sys_lep,2);

       
   //cout<<"final sys_lep sys_pho "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<endl; 
   
}/////end of ID sys


/////////get HLT sys for LT
/////////////////////ID and PU sys for LT lepton tagged//////////////////////////

void statAn::getHLTsys(const char* filename1, const char* filename2, double trigEff1 , double trigEff2, double &sys_lep, int cat, bool usewei)
{
  std::cout<<"Inside getHLTsys; cat is "<<cat <<std::endl;
  cout<<"cat is and useweigh is "<<cat<<" "<<usewei<<endl;

  sys_lep = 0;

  cout<<"trigEff1 trigEff2 "<<trigEff1<<" "<<trigEff2<<endl;
   // open file and get requested minitree

  TChain *tree = new TChain("minitree");
  tree->Add(filename1);
  tree->Add(filename2);
  
  if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category;
   float hzg_mass;
   float puwei;
   float hltSys_lep;
  

   double nev = 0;
   double nev_raw = 0;
   

   double tmpsys_lep = 0;

   
   // associate tree branches with variables
   if (tree->SetBranchAddress("category", &category) != 0)
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }
   if (tree->SetBranchAddress("hzg_mass", &hzg_mass) != 0)
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   //if (usewei && tree->SetBranchAddress("puwei", &puwei) != 0)
   if (usewei && tree->SetBranchAddress("mcwei", &puwei) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }
   
   if (tree->SetBranchAddress("hltSys_lep", &hltSys_lep) != 0) ///SJ
     {
       //FATAL("TTree::SetBranchAddress() failed");
     }

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if ( (cat > 0  && cat!=6789) && category != cat)
       continue;
     
     if (cat==0 && !(category>=1 && category<=10) )
       continue;

     if (cat==6789 && !(category>=6 && category<=9) ) ///lepton tagged
       continue;


     if(hzg_mass<xmin || hzg_mass>xmax) continue;
     
     double trigEff = 1;
     string fname(tree->GetCurrentFile()->GetName());
     
     int fstr;
     fstr = fname.find("mu",0);
     if(fstr!=string::npos)
       {
	 trigEff = trigEff2;
	 hltSys_lep = 1.3/100;
       }
     
     fstr = fname.find("ele",0);
     if(fstr!=string::npos)
       trigEff = trigEff1;



     nev += puwei*trigEff;
     
     nev_raw++;

     sys_lep += hltSys_lep*trigEff;

     
     tmpsys_lep += hltSys_lep*trigEff;
     
   }//for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) 

   ///bec in the end it is sqrt
   sys_lep = pow(sys_lep,2);

       
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



  ///with 2.66 fb-1 lumi
  /*
  double norm_eeg[2][11] = {
    {6.99921e-05, 6.25682e-05, 5.70901e-05, 5.28857e-05, 4.95557e-05, 4.67455e-05, 5.38892e-05, 6.26585e-05, 7.37110e-05, 8.79976e-05, 1.07473e-04}, 
    {1.17197e-05, 1.04330e-05, 9.50575e-06, 8.80897e-06, 8.26751e-06, 7.81746e-06, 9.09090e-06, 1.06691e-05, 1.26792e-05, 1.53112e-05, 1.89609e-05}};
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

  //{"ggH", "VBF", "WplusH", "WminusH", "ZH","ttH"}
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

  /*
  double norm_eeg[6][11] = {
    {2.11478e-03, 2.24872e-03, 2.37852e-03, 2.50662e-03, 2.63312e-03, 2.75275e-03, 2.86643e-03, 2.97609e-03, 3.08076e-03, 3.17547e-03, 3.26395e-03},
    {1.59354e-04, 1.71007e-04, 1.82536e-04, 1.94085e-04, 2.05704e-04, 2.16918e-04, 2.27130e-04, 2.37076e-04, 2.46701e-04, 2.55598e-04, 2.64125e-04},
    {3.87179e-05, 2.91123e-05, 2.36827e-05, 2.01498e-05, 1.76907e-05, 1.58447e-05, 1.88806e-05, 2.30382e-05, 2.90243e-05, 3.85125e-05, 5.61572e-05},
    {2.52410e-05, 1.88637e-05, 1.52930e-05, 1.30202e-05, 1.14052e-05, 1.02038e-05, 1.20848e-05, 1.46049e-05, 1.82022e-05, 2.37422e-05, 3.34815e-05},
    {4.02512e-05, 3.07220e-05, 2.52567e-05, 2.16839e-05, 1.91539e-05, 1.72218e-05, 2.05429e-05, 2.49733e-05, 3.13497e-05, 4.12918e-05, 5.92335e-05},
    {2.31706e-04, 2.04980e-04, 1.85018e-04, 1.69730e-04, 1.57444e-04, 1.46963e-04, 1.70590e-04, 2.00069e-04, 2.38749e-04, 2.91748e-04, 3.69797e-04}};
  */


  ///preapproval
  /*
  double norm_eeg[6][11] = {
    {7.04918e-04, 6.65316e-04, 6.32514e-04, 6.05239e-04, 5.82129e-04, 5.61145e-04, 6.32915e-04, 7.16845e-04, 8.16361e-04, 9.35264e-04, 1.08215e-03}, 
    {1.59354e-04, 1.22036e-04, 1.01120e-04, 8.77570e-05, 7.84911e-05, 7.15310e-05, 8.64925e-05, 1.06814e-04, 1.36073e-04, 1.81726e-04, 2.64125e-04}, 
    {3.87179e-05, 2.91123e-05, 2.36827e-05, 2.01498e-05, 1.76907e-05, 1.58447e-05, 1.88806e-05, 2.30382e-05, 2.90243e-05, 3.85125e-05, 5.61572e-05}, 
    {2.52410e-05, 1.88637e-05, 1.52930e-05, 1.30202e-05, 1.14052e-05, 1.02038e-05, 1.20848e-05, 1.46049e-05, 1.82022e-05, 2.37422e-05, 3.34815e-05}, 
    {4.02512e-05, 3.07220e-05, 2.52567e-05, 2.16839e-05, 1.91539e-05, 1.72218e-05, 2.05429e-05, 2.49733e-05, 3.13497e-05, 4.12918e-05, 5.92335e-05}, 
    {2.31706e-04, 2.04980e-04, 1.85018e-04, 1.69730e-04, 1.57444e-04, 1.46963e-04, 1.70590e-04, 2.00069e-04, 2.38749e-04, 2.91748e-04, 3.69797e-04}};
  */

  ///21st feb - wrong nentries gor VH. assumed its only llg
  /*
  double norm_eeg[6][11] = {
    {1.05740e-03, 9.41565e-04, 8.56445e-04, 7.91543e-04, 7.40280e-04, 6.97301e-04, 8.04577e-04, 9.36773e-04, 1.10396e-03, 1.32108e-03, 1.61892e-03}, 
    {1.59357e-04, 1.42053e-04, 1.29552e-04, 1.20142e-04, 1.12821e-04, 1.06727e-04, 1.24168e-04, 1.45805e-04, 1.73399e-04, 2.09592e-04, 2.59898e-04}, 
    {3.87179e-05, 2.91122e-05, 2.36826e-05, 2.01497e-05, 1.76905e-05, 1.58446e-05, 1.88513e-05, 2.29535e-05, 2.88279e-05, 3.80623e-05, 5.49912e-05}, 
    {2.62923e-05, 1.94191e-05, 1.56419e-05, 1.32629e-05, 1.15852e-05, 1.03435e-05, 1.22372e-05, 1.47677e-05, 1.83668e-05, 2.38794e-05, 3.34821e-05}, 
    {4.02512e-05, 3.04737e-05, 2.49380e-05, 2.13473e-05, 1.88179e-05, 1.68941e-05, 2.01812e-05, 2.45823e-05, 3.09480e-05, 4.09493e-05, 5.92347e-05}, 
    {2.31706e-04, 2.03399e-04, 1.82572e-04, 1.66785e-04, 1.54207e-04, 1.43564e-04, 1.65090e-04, 1.91331e-04, 2.24795e-04, 2.68912e-04, 3.30496e-04}};
  */

  ///21st feb - corrected for correct llg
  /*
  double norm_eeg[6][11] = {
    {1.05740e-03, 9.41565e-04, 8.56445e-04, 7.91543e-04, 7.40280e-04, 6.97301e-04, 8.04577e-04, 9.36773e-04, 1.10396e-03, 1.32108e-03, 1.61892e-03}, 
    {1.59357e-04, 1.42053e-04, 1.29552e-04, 1.20142e-04, 1.12821e-04, 1.06727e-04, 1.24168e-04, 1.45805e-04, 1.73399e-04, 2.09592e-04, 2.59898e-04}, 
    {1.94195e-04, 1.46132e-04, 1.18929e-04, 1.01216e-04, 8.88807e-05, 7.96176e-05, 9.46913e-05, 1.15239e-04, 1.44624e-04, 1.90728e-04, 2.74965e-04}, 
    {1.31450e-04, 9.65920e-05, 7.75886e-05, 6.56737e-05, 5.72980e-05, 5.11124e-05, 6.04983e-05, 7.30544e-05, 9.09409e-05, 1.18402e-04, 1.66428e-04}, 
    {9.97957e-05, 7.59082e-05, 6.22829e-05, 5.34053e-05, 4.71327e-05, 4.23509e-05, 5.05117e-05, 6.13945e-05, 7.70510e-05, 1.01446e-04, 1.45420e-04}, 
    {2.31706e-04, 2.03399e-04, 1.82572e-04, 1.66785e-04, 1.54207e-04, 1.43564e-04, 1.65090e-04, 1.91331e-04, 2.24795e-04, 2.68912e-04, 3.30496e-04}};
  */

  ///21st feb - corrected for genweights
  /*double norm_eeg[6][11] = {
    {1.05740e-03, 9.41565e-04, 8.56445e-04, 7.91543e-04, 7.40280e-04, 6.97301e-04, 8.04577e-04, 9.36773e-04, 1.10396e-03, 1.32108e-03, 1.61892e-03}, 
    {1.59606e-04, 1.42289e-04, 1.29776e-04, 1.20356e-04, 1.13026e-04, 1.06925e-04, 1.24396e-04, 1.46070e-04, 1.73710e-04, 2.09963e-04, 2.60346e-04}, 
    {4.14283e-04, 3.11390e-04, 2.53263e-04, 2.15454e-04, 1.89143e-04, 1.69395e-04, 2.01428e-04, 2.45071e-04, 3.07447e-04, 4.05208e-04, 5.83520e-04}, 
    {2.79648e-04, 2.04952e-04, 1.64398e-04, 1.39029e-04, 1.21225e-04, 1.08091e-04, 1.27921e-04, 1.54441e-04, 1.92200e-04, 2.50130e-04, 3.51318e-04}, 
    {2.22915e-04, 1.70011e-04, 1.39706e-04, 1.19910e-04, 1.05899e-04, 9.52025e-05, 1.13497e-04, 1.37867e-04, 1.72875e-04, 2.27296e-04, 3.25012e-04}, 
    {2.40899e-04, 2.11106e-04, 1.89259e-04, 1.72735e-04, 1.59595e-04, 1.48495e-04, 1.70852e-04, 1.98141e-04, 2.32998e-04, 2.79048e-04, 3.43517e-04}};
  */




  //{"ggH", "VBF", "WplusH", "WminusH", "ZH","ttH"}                         
  ///22st feb - corrected for ZH norm                         
  /*double norm_eeg[6][11] = {
    {1.05740e-03, 9.41565e-04, 8.56445e-04, 7.91543e-04, 7.40280e-04, 6.97301e-04, 8.04577e-04, 9.36773e-04, 1.10396e-03, 1.32108e-03, 1.61892e-03},
    {1.59606e-04, 1.42289e-04, 1.29776e-04, 1.20356e-04, 1.13026e-04, 1.06925e-04, 1.24396e-04, 1.46070e-04, 1.73710e-04, 2.09963e-04, 2.60346e-04},
    {4.14283e-04, 3.11390e-04, 2.53263e-04, 2.15454e-04, 1.89143e-04, 1.69395e-04, 2.01428e-04, 2.45071e-04, 3.07447e-04, 4.05208e-04, 5.83520e-04},
    {2.79648e-04, 2.04952e-04, 1.64398e-04, 1.39029e-04, 1.21225e-04, 1.08091e-04, 1.27921e-04, 1.54441e-04, 1.92200e-04, 2.50130e-04, 3.51318e-04},
    {4.25516e-04, 3.23708e-04, 2.65624e-04, 2.27775e-04, 2.01029e-04, 1.80639e-04, 2.15330e-04, 2.61530e-04, 3.27872e-04, 4.30951e-04, 6.15873e-04},
    {2.40899e-04, 2.11106e-04, 1.89259e-04, 1.72735e-04, 1.59595e-04, 1.48495e-04, 1.70852e-04, 1.98141e-04, 2.32998e-04, 2.79048e-04, 3.43517e-04}}; 

  */


  //{"ggH", "VBF", "WH", "ZH","ttH"}                         
  ////W+H W-H combined
  /*  double norm_eeg[5][11] = {
    {1.05740e-03, 9.41565e-04, 8.56445e-04, 7.91543e-04, 7.40280e-04, 6.97301e-04, 8.04577e-04, 9.36773e-04, 1.10396e-03, 1.32108e-03, 1.61892e-03}, 
    {1.59606e-04, 1.42289e-04, 1.29776e-04, 1.20356e-04, 1.13026e-04, 1.06925e-04, 1.24396e-04, 1.46070e-04, 1.73710e-04, 2.09963e-04, 2.60346e-04}, 
    {3.48859e-04, 2.59159e-04, 2.09340e-04, 1.77558e-04, 1.55386e-04, 1.38876e-04, 1.64699e-04, 1.99707e-04, 2.49532e-04, 3.26632e-04, 4.64753e-04}, 
    {4.25516e-04, 3.23708e-04, 2.65624e-04, 2.27775e-04, 2.01029e-04, 1.80639e-04, 2.15330e-04, 2.61530e-04, 3.27872e-04, 4.30951e-04, 6.15873e-04}, 
    {2.40899e-04, 2.11106e-04, 1.89259e-04, 1.72735e-04, 1.59595e-04, 1.48495e-04, 1.70852e-04, 1.98141e-04, 2.32998e-04, 2.79048e-04, 3.43517e-04}};

  */

  /*
  ////7th march
  //{"ggH", "VBF", "WH", "ZH","ttH"}  - luminosity updated to 35.86 fb-1
  double norm_eeg[5][11] = {
    {1.04000e-03, 9.26070e-04, 8.42351e-04, 7.78517e-04, 7.28098e-04, 6.85826e-04, 7.91337e-04, 9.21358e-04, 1.08579e-03, 1.29934e-03, 1.59228e-03}, 
    {1.56980e-04, 1.39947e-04, 1.27641e-04, 1.18376e-04, 1.11166e-04, 1.05165e-04, 1.22349e-04, 1.43667e-04, 1.70852e-04, 2.06508e-04, 2.56062e-04}, 
    {3.43118e-04, 2.54894e-04, 2.05896e-04, 1.74636e-04, 1.52829e-04, 1.36591e-04, 1.61989e-04, 1.96421e-04, 2.45425e-04, 3.21257e-04, 4.57105e-04}, 
    {4.18514e-04, 3.18381e-04, 2.61253e-04, 2.24026e-04, 1.97721e-04, 1.77666e-04, 2.11786e-04, 2.57226e-04, 3.22476e-04, 4.23859e-04, 6.05738e-04}, 
    {2.36935e-04, 2.07632e-04, 1.86144e-04, 1.69892e-04, 1.56969e-04, 1.46052e-04, 1.68040e-04, 1.94880e-04, 2.29163e-04, 2.74456e-04, 3.37864e-04}};
  */

  /*
  ///16th march - update to reminiAOD samples - its the same as before
  double norm_eeg[5][11] = {
    {1.04000e-03, 9.26070e-04, 8.42351e-04, 7.78517e-04, 7.28098e-04, 6.85826e-04, 7.91337e-04, 9.21358e-04, 1.08579e-03, 1.29934e-03, 1.59228e-03}, 
    {1.56980e-04, 1.39947e-04, 1.27641e-04, 1.18376e-04, 1.11166e-04, 1.05165e-04, 1.22349e-04, 1.43667e-04, 1.70852e-04, 2.06508e-04, 2.56062e-04}, 
    {3.43118e-04, 2.54894e-04, 2.05896e-04, 1.74636e-04, 1.52829e-04, 1.36591e-04, 1.61989e-04, 1.96421e-04, 2.45425e-04, 3.21257e-04, 4.57105e-04}, 
    {4.18514e-04, 3.18381e-04, 2.61253e-04, 2.24026e-04, 1.97721e-04, 1.77666e-04, 2.11786e-04, 2.57226e-04, 3.22476e-04, 4.23859e-04, 6.05738e-04}, 
    {2.36935e-04, 2.07632e-04, 1.86144e-04, 1.69892e-04, 1.56969e-04, 1.46052e-04, 1.68040e-04, 1.94880e-04, 2.29163e-04, 2.74456e-04, 3.37864e-04}};
  */


  ////6th april, 2017 - luminosity 35.86--->35.9
  /*double norm_eeg[5][11] = {
    {1.04116e-03, 9.27103e-04, 8.43290e-04, 7.79385e-04, 7.28910e-04, 6.86591e-04, 7.92219e-04, 9.22385e-04, 1.08700e-03, 1.30079e-03, 1.59405e-03}, 
    {1.57155e-04, 1.40103e-04, 1.27783e-04, 1.18508e-04, 1.11290e-04, 1.05283e-04, 1.22486e-04, 1.43827e-04, 1.71042e-04, 2.06738e-04, 2.56348e-04}, 
    {3.43501e-04, 2.55178e-04, 2.06125e-04, 1.74831e-04, 1.52999e-04, 1.36743e-04, 1.62169e-04, 1.96640e-04, 2.45699e-04, 3.21616e-04, 4.57615e-04}, 
    {4.18981e-04, 3.18736e-04, 2.61544e-04, 2.24276e-04, 1.97941e-04, 1.77864e-04, 2.12023e-04, 2.57513e-04, 3.22836e-04, 4.24332e-04, 6.06414e-04}, 
    {2.37199e-04, 2.07864e-04, 1.86352e-04, 1.70082e-04, 1.57144e-04, 1.46215e-04, 1.68228e-04, 1.95098e-04, 2.29419e-04, 2.74762e-04, 3.38241e-04}};
  */


  //23rd april - new samples in v2 which was run to produce OR of triggers - this is actually the same as before
  double norm_eeg[5][11] = {
    {1.04116e-03, 9.27103e-04, 8.43290e-04, 7.79385e-04, 7.28910e-04, 6.86591e-04, 7.92219e-04, 9.22385e-04, 1.08700e-03, 1.30079e-03, 1.59405e-03}, 
    {1.57155e-04, 1.40103e-04, 1.27783e-04, 1.18508e-04, 1.11290e-04, 1.05283e-04, 1.22486e-04, 1.43827e-04, 1.71042e-04, 2.06738e-04, 2.56348e-04}, 
    {3.43501e-04, 2.55178e-04, 2.06125e-04, 1.74831e-04, 1.52999e-04, 1.36743e-04, 1.62169e-04, 1.96640e-04, 2.45699e-04, 3.21616e-04, 4.57615e-04}, 
    {4.18981e-04, 3.18736e-04, 2.61544e-04, 2.24276e-04, 1.97941e-04, 1.77864e-04, 2.12023e-04, 2.57513e-04, 3.22836e-04, 4.24332e-04, 6.06414e-04}, 
    {2.37199e-04, 2.07864e-04, 1.86352e-04, 1.70082e-04, 1.57144e-04, 1.46215e-04, 1.68228e-04, 1.95098e-04, 2.29419e-04, 2.74762e-04, 3.38241e-04}};
  


  double fac_gstar = 0.969; ///from Ming-Yan's study

  int index = mass - 120;
  
   if (index < 0 || index > 11)
      FATAL("wrong mass given");
   /*
   if (TString(channel).Contains("ee"))
      return norm_eeg[p][index]/fac_gstar;
   */
   
   return norm_eeg[p][index]/fac_gstar;
}


/////not needed for now
/*
///added on 21st jan
double statAn::ngen_Tot(const char* channel, int p, int mass)
{
  //* Values were calculated by normfactors_and_theoretical_uncertainties.py.
  //*
  //* channel = either of "eeg" (for electrons) or "mmg" (for muons);
  //* p = production process number (0=ggH, 1=qqH, 2=WH, 3=ZH, 4=ttH);
  //* mass = value of SM Higgs mass (from the region 120-150GeV with step 1GeV).
  //*

  //{"ggH", "VBF", "WplusH", "WminusH", "ZH","ttH"}
  double ngen_Tot[6][11] = {
    {9.99980e+04, 9.99184e+04, 9.98388e+04, 9.97592e+04, 9.96796e+04, 9.96000e+04, 9.95196e+04, 9.94392e+04, 9.93588e+04, 9.92784e+04, 9.91980e+04}, 
    {1.00000e+05, 9.96800e+04, 9.93600e+04, 9.90400e+04, 9.87200e+04, 9.84000e+04, 9.83996e+04, 9.83992e+04, 9.83988e+04, 9.83984e+04, 9.83980e+04}, 
    {9.99920e+04, 1.39827e+05, 1.79662e+05, 2.19496e+05, 2.59331e+05, 2.99166e+05, 2.41262e+05, 1.83358e+05, 1.25454e+05, 6.75500e+04, 9.64600e+03}, 
    {9.77400e+04, 1.37120e+05, 1.76500e+05, 2.15879e+05, 2.55259e+05, 2.94639e+05, 2.55671e+05, 2.16703e+05, 1.77734e+05, 1.38766e+05, 9.97980e+04}, 
    {9.99960e+04, 1.37929e+05, 1.75863e+05, 2.13796e+05, 2.51730e+05, 2.89663e+05, 2.50789e+05, 2.11915e+05, 1.73040e+05, 1.34166e+05, 9.52920e+04}, 
    {9.95700e+03, 8.16540e+03, 6.37380e+03, 4.58220e+03, 2.79060e+03, 9.99000e+02, 9.90000e+02, 9.81000e+02, 9.72000e+02, 9.63000e+02, 9.54000e+02}};

  int index = mass - 120;
   // if (index < 0 || index > 30)
   if (index < 0 || index > 11)
      FATAL("wrong mass given");

   if (TString(channel).Contains("ee"))
      return ngen_Tot[p][index];

   return ngen_Tot[p][index];
}
*/


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
  //double lumi =2689.8; //pb-1
  double lumi =mylumi; //pb-1

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
void statAn::mybkgfit(double xmin, double xmax, int icat, string schannel){

  rfit = NULL;

  //RooRealVar  fX("x", "", xmin, xmax); 


   PdfModelBuilder pdfsModel;
   pdfsModel.setObsVar(fX);

  RooExtendPdf *bgrfit_ext = new RooExtendPdf();

  cout<<"Model is "<<model<<endl;
  
  
  else if(model=="RooPolynomial"){
    

    coeflist = new RooArgList;
    cout<<"inside RooPolynomial"<<endl;

    RooAbsPdf *bgrfit1 = pdfsModel.getBernstein(Form("Bernstein_%d_model%s",poldeg,model.c_str()),poldeg);
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+6);    
    //fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 200, 0, 1e+10);    
    fPar[0] = new RooRealVar(TString::Format("norm_model%s",model.c_str()), "", 500, 0, 1e+10);    
    
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


  
  
  ////Exp
  else if(model=="Exp"){
    

    cout<<"inside Exp"<<endl;

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
    

    name[0]= "mean";
    name[1] = "sigma";
    name[2] = "stepval";
    
    for(int ipol=1; ipol<=poldeg; ipol++){  
      name[ipol+2] = TString::Format("coef%d",ipol);
    }
  }//if(model=="RooGaussStepBernstein")



  if(model=="RooBernstein"){   
    cout<<"inside RooBernstein"<<endl;
    

    name[0]= "mean";
    name[1] = "sigma";
    name[2] = "stepval";
    
    for(int ipol=1; ipol<=poldeg; ipol++){  
      name[ipol+2] = TString::Format("coef%d",ipol);
    }
  }//if(model=="RooBernstein")

  
  return name[ipar];
}

int statAn::getNbkgPar(){

  if(model=="RooGaussStepBernstein"){
    
    return poldeg+3;
  }
  
  
  if(model=="RooBernstein"){
    
    return poldeg;
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
    //pdf2 = new RooGaussian("subpdf2", "", *fX, *sigmean, *sigma2);
    sigfit1 = new RooAddPdf("CBGaussian", "", *pdf1, *pdf2, *frac);


    cout<<"----values inside mysigfunc-----"<<endl;
    cout<<"meana : sigmaa : alpha : power : mean1 : sigma2 : frac "<<sigmean->getVal()<<" "<<sigsigma->getVal()<<" "<<alpha->getVal()<<" "<<power->getVal()<<" "<<sigmean->getVal()<<" "<<sigma2->getVal()<<" "<<frac->getVal()<<endl;
    
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

    if( (icat==1 && np==1  && nm==10 && schannel=="mmg") || (icat==4&&np==0&&nm==10&&schannel=="eeg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.1, 20);

    if( (icat==4&&np==1&&nm==10&&schannel=="eeg") || (icat==4&&np==0&&nm==10&&schannel=="eeg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 8);

    if( (icat==4&&np==0&&nm==5) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 5);

    if( (icat==4&&np==0&&nm==5 && schannel=="eeg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 50);

    if( (icat==1&&np==0&&nm==5 && schannel=="mmg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 50);

    if( (icat==3&&np==0&&nm==0) &&  (schannel=="mmg"|| schannel=="eeg")) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 5);


    if( (icat==3&&np==1&&nm==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 5);
    
    if( (icat==3&&np==1&&nm==5) && schannel=="eeg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 0.5, 0.1, 5);


    if( (icat==1 && np==1  && nm==5 && schannel=="mmg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 3, 0.1, 30);

    //if( (icat==3 && np==1  && nm==0 && schannel=="eeg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 3, 0.1, 30);

    if( (icat==3 && np==1  && nm==0 && schannel=="eeg") ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 3, 0.1, 5);

    //if( (icat==1 && np==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 20);

    //if( (icat==1 && np==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 20);

    ///try
    if( (icat==2||icat==5) &&(np==0||np==1) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 100);
    //if( (icat==5) &&(np==0||np==1) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 100);

    //if( (icat==1 && np==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 100);

    
    if( (icat==5 && np==0) && schannel=="mmg" && nm==0) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30);

    
    if( icat==2 &&np==1 && schannel=="mmg" && nm==10 ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30);
    if( icat==3 &&np==0 && schannel=="mmg" && nm==0 ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30);
    
    if( (np==0&&schannel=="mmg"&&icat==10 && nm==10) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30); 



    if( (np==1&&schannel=="mmg"&&icat==1 && nm==5) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30); 

    if( icat==3 && np==1 &&schannel=="mmg" && nm==0 ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30); 
    
    if( icat==2 &&np==0 && schannel=="mmg" && (nm==0||nm==5) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30);

    if( icat==1 && np==1 && schannel=="eeg" && nm==10)  sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30); 


    ////23rd april
    if( icat==3 && np==1  && nm==10 && schannel=="eeg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 5, 0.1, 100);

    if( icat==2 && np==0  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 5, 0.1, 100);

    //if( icat==1 && np==0  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 5, 0.1, 100);



    

    ///22jan  
    //if( (icat==1 && np==0 && nm==0) && schannel=="mmg") sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 20);
    
    //if( (icat==3 && np==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 1, 20);
    

    if( (icat==5 && np==0) && schannel=="mmg" && nm==5) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 100);

    ///9th may
    if( (icat==5 && np==1) && schannel=="mmg" && nm==10) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30);

    //if( icat==1 && np==0  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.01, 100);


    if( icat==5 && np==1  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 10);

    if( icat==3 && np==1  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 50);

    if( icat==5 && np==0  && nm==0 && schannel=="eeg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 30);


    if( icat==5 && np==0  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.001, 70);

    if( (np==0&&schannel=="mmg"&&icat==10 && nm==0) ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 10); 
    //if( icat==1 && np==0  && nm==0 && schannel=="mmg" ) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 10, 0.1, 100);

    //if( (icat==5 && np==1) && schannel=="mmg" && nm==10) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 10);

    ///17th july
    if( (icat==5 && np==1) && schannel=="mmg" && nm==10) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 2, 0.5, 30);


    if( (icat==3 && np==0) && schannel=="eeg" && nm==0) sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 1, 0.01, 10);


    if( icat==1 && np==0 &&schannel=="mmg" && nm==10 )  sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 1, 0.01, 30);

    //17th sept
    if( icat==5 && np==0 &&schannel=="eeg" && nm==0 )  sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 1, 0.001, 30);
    if( icat==10 && np==0 &&schannel=="mmg" && nm==5 )  sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 1, 0.001, 30);
    if( icat==2 && np==0 &&schannel=="mmg" && nm==5 )  sigma1  = new RooRealVar(Form("sig_sigma1_chan%d_m%d_cat%d",np,nm,icat),"", 1, 0.001, 30);


    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.5, 10);
    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.01, 20);
    alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.01, 50);

    //if( (icat==1&&np==0)  && nm == 0 && schannel=="mmg") alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",0.1, 0.01, 6);
    if( (icat==1&&np==0)  && nm == 0 && schannel=="mmg" ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",0.1, 0.01, 6);

    if( (np==1&&schannel=="eeg"&&nm==5&&icat==3) || (np==0&&schannel=="eeg"&&nm==5&&icat==2) ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",0.1, 0.01, 100);


    if((np==0&&schannel=="eeg"&&nm==0&&icat==2) ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",0.1, 0.01, 100);
    
    if( icat==1 && np==0 &&schannel=="mmg" && nm==10 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.01, 100);



    if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",4, 0.01, 100);


    if( icat==3 && np==0 &&schannel=="mmg" && nm==0 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",20, 0.01, 100);

    if( icat==2 && np==0 &&schannel=="mmg" && nm==0 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.01, 300);

    if( icat==2 && np==0 &&schannel=="eeg" && nm==0 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.01, 300);

    if( icat==3 && np==0 &&schannel=="eeg" && nm==0 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",30, 0.01, 300);

    if( icat==5 && np==1 &&schannel=="mmg" && nm==0 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.01, 300);

    //if( icat==10 && np==0 &&schannel=="mmg" && nm==5 ) alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "", 4,0.01, 300);
    
    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.001, 100);

    //if( (icat==2&&np==0)  && nm == 5 )  alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",0.1, 0.01, 50);
    //alpha   = new RooRealVar(Form("alpha_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.01, 50);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.5, 10);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.01, 20);
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.001, 30);

    ///26th jan
    //power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.001, 70);
    power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",30, 29,31);

    
    if( icat==1 && np==0 &&schannel=="mmg" && nm==10 )  power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.36794,1,2);
    
    if( icat==1 && np==0 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",0.9, 0.8,100);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==0 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.31812, 1,2);
    

    if( icat==2 && np==0 &&schannel=="mmg" && nm==10 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",0.99, 0.7,100);

    if( icat==2 && np==0 &&schannel=="mmg" && nm==0 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.88191, 1,2);
    if( icat==2 && np==0 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.88191, 1,2);

    if( icat==3 && np==0 &&schannel=="mmg" && nm==10 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",0.99,0.5,100);

    if( icat==3 && np==0 &&schannel=="mmg" && nm==0 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.6,1,2);
    if( icat==3 && np==0 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.60032,1,2);

    if( icat==1 && np==1 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.72837,1,2);

    if( icat==1 && np==1 &&schannel=="mmg" && nm==0 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.72837,1,2);
    if( icat==1 && np==1 &&schannel=="mmg" && nm==10 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1.72837,1,2);
    
    
    if( icat==5 && np==1 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",2.14829,2,3);
    if( icat==5 && np==1 &&schannel=="mmg" && nm==0 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",2.14829,2,3);
    if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",2.14829,2,3);

    if( icat==10 && np==1 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",0.99,0.7,100);
    //if( icat==10 && np==1 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3.11633,3,4);
    if( icat==10 && np==1 &&schannel=="mmg" && nm==0 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3.11633,3,4);
    if( icat==10 && np==1 &&schannel=="mmg" && nm==10 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",3.11633,3,4);


    if( icat==10 && np==0 &&schannel=="mmg" && nm==5 ) power   = new RooRealVar(Form("power_chan%d_m%d_cat%d",np,nm,icat), "",1,0.7,100);


    
    
    s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 0.1, 30);       // sigma2/sigma1

    
    ///25th jan
    //if( (icat==3 && np==0 && nm==10) || icat==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 2, 30);       // sigma2/sigma1
    //if( (icat==3 && np==0 && nm==10) || icat==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 2, 30);
    if( (icat==3 && np==0 && nm==10) ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 2, 30);

    //if( (icat==4 && np==0 && nm==0 && schannel=="mmg") ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.6, 30);

    if( (icat==4 && np==0 && nm==0) ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.6, 30);
  

    if( (icat==1 && np==0 && nm==10 &&schannel=="eeg") ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 0.1, 30);
  
    if( (icat==3 && np==0 && nm==10 &&schannel=="mmg") ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",6, 0.1, 30);

    if( ((icat==2||icat==1||icat==5) && (np==0||np==1) &&schannel=="mmg") ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 50);
    //if( ((icat==2||icat==1||icat==5 || icat==10) && (np==0||np==1) &&schannel=="mmg") ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 100);

    if( ((icat==5&&np==1) || (icat==1&&np==0) ) &&schannel=="mmg")  s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 100);


    if( icat==5 && np==0 &&schannel=="mmg" ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 100);

    if( icat==2 && np==1 &&schannel=="mmg") s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 100);
    if( icat==3 && np==0 &&schannel=="mmg") s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 100);

    //if( icat==10 && (np==0||np==1) &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);
    //if( icat==1 && (np==0||np==1) &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    if( icat==2 && np==0 &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);


    if( icat==10 && np==0 &&schannel=="mmg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==5 && np==0 &&schannel=="mmg") s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==5 && np==0 &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 100);


    if( icat==3 && np==0 &&schannel=="mmg" && nm==0) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    //if( icat==2 && np==1 &&schannel=="mmg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",30, 0.1, 100);


    ////16th march//

    if( icat==1 && np==1 &&schannel=="eeg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==3 && np==1 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==1 && np==1 &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==2 && np==0 &&schannel=="eeg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    //if( icat==5 && np==0 &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    
    if( icat==10 && np==0 &&schannel=="mmg" && (nm==5||nm==0 || nm==10) ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    
    if( icat==10 && np==0 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==2 && np==0 &&schannel=="mmg" && (nm==5||nm==0)) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==1 && np==1 && schannel=="eeg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==1 && np==1 && schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    
    if( icat==1 && np==0 && schannel=="mmg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    
    //if( icat==3 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==5 && np==1 &&schannel=="mmg" && nm==0) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    

    if( icat==1 && np==0 &&schannel=="mmg" && nm==0) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    if( icat==3 && np==0 &&schannel=="mmg" && nm==10 ) s21     =  new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);


    if( icat==5 && np==0 &&schannel=="mmg" && nm==5) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    /////23th april
    if( icat==1 && np==1 &&schannel=="eeg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);
    if( icat==3 && np==0 &&schannel=="eeg" && nm==0) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==10 && np==0 &&schannel=="eeg" && nm==0) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==3 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);
    if( icat==5 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    if( icat==2 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==10) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);
    

    /////9th may
    if( icat==3 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);
    if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);
 
    if( icat==2 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==2 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",1, 0.001, 50);

    if( icat==1 && np==0 &&schannel=="eeg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);

    if( icat==5 && np==0 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);


    if( icat==3 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==5 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);


    if( icat==5 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==10 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==2 && np==0 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    //if( icat==5 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 200);

    if( icat==5 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==3 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==2 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==5 && np==0 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 30);


    if( icat==5 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 70);


    if( icat==2 && np==1 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==5 && np==1 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==3 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==10 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==5 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    //if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 10);

    if( icat==2 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 10);

    if( icat==2 && np==0 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 10);

    ///17th july


    if( icat==3 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==3 && np==0 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 100);

    if( icat==10 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 100);

    if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);

    ///27th july
    if( icat==4 && np==0 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);
    if( icat==4 && np==0 &&schannel=="eeg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",3, 0.1, 50);


    ///6th aug
    if( icat==5 && np==1 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);
    if( icat==3 && np==1 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);

    if( icat==2 && np==0 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);
    if( icat==3 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);
    if( icat==10 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);


    if( icat==5 && np==0 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);
    if( icat==3 && np==0 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);
    if( icat==1 && np==0 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);
    if( icat==2 && np==0 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);

    //###7th aug
    if( icat==5 && np==1 &&schannel=="eeg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);
    if( icat==4 && np==0 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);
    if( icat==5 && np==0 &&schannel=="eeg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200); ///try again

    if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);

    if( icat==1 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);
    
    if( icat==1 && np==0 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);

    if( icat==2 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 70);

    if( icat==3 && np==0 &&schannel=="mmg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);

    if( icat==5 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);


    if( icat==10 && np==1 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);

    
    if( icat==2 && np==0 &&schannel=="mmg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",2, 0.1, 10);

    if( icat==4 && np==0 &&schannel=="eeg" && nm==10 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",10, 0.1, 200);


    if( icat==5 && np==1 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",7, 0.1, 200);


    ///17th sept
    if( icat==5 && np==0 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.1, 50);

    if( icat==10 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);
    if( icat==2 && np==0 &&schannel=="mmg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",15, 0.1, 200);


    ///31 oct
    if( icat==3 && np==1 &&schannel=="eeg" && nm==0 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.1, 200);

    if( icat==4 && np==0 &&schannel=="eeg" && nm==5 ) s21     = new RooRealVar(Form("s21_chan%d_m%d_cat%d",np,nm,icat), "",5, 0.1, 200);

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
    

    fPar_sig[0] = new RooRealVar( Form("norm_chan%d_m%d_cat%d",np,nm,icat), "", 200, 0, 1e+6);    
    RooExtendPdf *sigfit_ext = new  RooExtendPdf("sigfit", "", *sigfit1, *fPar_sig[0], "our_window");
    
    
    cout<<"inside mysigfit np and nm "<<np<<" "<<nm<<endl;
    
    /*
    rsig[np][nm] = sigfit_ext->fitTo(*rooMassSig[np][nm],
				     RooFit::PrintLevel(2),
				     RooFit::Range(xmin, xmax),
				     RooFit::NumCPU(4),
				     RooFit::SumW2Error(kTRUE),
				     RooFit::Save(true)
				     //	     ,RooFit::Extended(true)
				     );
    */
    

    
    rsig[np][nm] = sigfit1->fitTo(*rooMassSig[np][nm],
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

    
    /*RooRealVar* delta21_fitresult1 = (RooRealVar*) rsig[np][nm]->floatParsFinal().find(Form("delta21_chan%d_m%d_cat%d",np,nm,icat));
    //cout<<"Found delta21_1"<<endl;

    RooRealVar* delta21_fitresult2 = (RooRealVar*) rsig[np][nm + 5]->floatParsFinal().find(Form("delta21_chan%d_m%d_cat%d",np,nm+5,icat));
    //cout<<"Found delta21_2"<<endl;
    
    double delta21_fr =  a * delta21_fitresult1->getVal() + b * delta21_fitresult2->getVal();
    cout<<"delta1, delta2 and calculated delta1_fr"<<delta21_fitresult1->getVal() <<" "<<delta21_fitresult2->getVal()<< " "<<delta21_fr<<endl;
    
    //double delta21_fr = 0;
    */

        
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    

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
  //int nBinsForMass = 55;
  //bool blind = true;
  bool blind = false;

  double num_zero_shapesys = 0.00001;

  //bool plotRooFITErrBands = true;
  bool plotRooFITErrBands = false;

  bool drawPulldistribution = false;
  //bool floatbkg = false;
  bool floatbkg = true;

  string sigmodel = "RooCBxGaus";
  bool _nuisance = false;
  //bool _nuisance = true;
  
  //double lumi = 2.7*1000;///pb-1
  //double lumi = 2689.8; //pb-1
  double lumi =mylumi; //pb-1
  // cat1, cat2, cat3 or cat4
  TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);


   ///read files for systematics
   ///ROCHOR
   ifstream fsys_mu; //sysMuon_hzg_cat1_8TeV.txt

   double sysrochor_mean_arr[10][20];
   double sysrochor_reso_arr[10][20];

     
   if(cat>0 && (schannel == "mmg" || schannel == "eeg_mmg") ){
     fsys_mu.open(TString::Format("sysMuon_hzg_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_mu.eof()){
       
       string proc;
       int mass;
       double sysrochor_mean, sysrochor_reso;
       
       fsys_mu >> proc >> mass >> sysrochor_mean >> sysrochor_reso;
       
       
       cout<<"reading sysMuon file "<< TString::Format("sysMuon_hzg_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean sysrochor_reso "<<proc<<" "<<mass<<" "<<sysrochor_mean<<" "<<sysrochor_reso<<endl;
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sysrochor_mean_arr[iproc][5*imass] = sysrochor_mean;
       sysrochor_reso_arr[iproc][5*imass] = sysrochor_reso;
       
     }//while(!datafile.eof())
   }//if(cat>0 && schannel == "mmg")



     
   
   //////////////////////electron, photons energy scale
   
   ///1. electron channel
   ifstream fsys_emscale_e;
   ifstream fsys_emreso_e; 
   
   double sysem_mean_arr_e[10][20];
   double sysem_reso_arr_e[10][20];
   
   
   ///scale
   if(cat>0 && (schannel == "eeg"||schannel=="eeg_mmg")){
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg") fsys_emscale_e.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     //if(schannel == "eeg_mmg") fsys_emscale_e.open(TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objelepho.txt",cats), ifstream::in );

     while(!fsys_emscale_e.eof()){
       
       string proc;
       int mass;
       double sysem_mean_e;
       
       fsys_emscale_e >> proc >> mass >> sysem_mean_e;
	 
       
       if(schannel == "eeg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_e<<endl;

       //if(schannel == "eeg_mmg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_e<<endl;
       
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sysem_mean_arr_e[iproc][5*imass] = sysem_mean_e;
       
     }//while(!datafile.eof())
     
     for(int iproc=0; iproc<3; iproc++){
       sysem_mean_arr_e[iproc][0] = sysem_mean_arr_e[iproc][5];
       sysem_mean_arr_e[iproc][10] = sysem_mean_arr_e[iproc][5];
     }

     ////reso
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg")  fsys_emreso_e.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     //if(schannel == "eeg_mmg")  fsys_emreso_e.open(TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     while(!fsys_emreso_e.eof()){
       
       string proc;
       int mass;
       double sysem_reso_e;
       
       fsys_emreso_e >> proc >> mass >> sysem_reso_e;
       
     
       if(schannel == "eeg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_e<<endl;
       //if(schannel == "eeg_mmg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_e<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_reso_arr_e[iproc][5*imass] = sysem_reso_e;
     
     }//while(!datafile.eof())

     for(int iproc=0; iproc<3; iproc++){
       sysem_reso_arr_e[iproc][0] = sysem_reso_arr_e[iproc][5];
       sysem_reso_arr_e[iproc][10] = sysem_reso_arr_e[iproc][5];
     }
     

   }//if(cat>0 && schannel == "eeg")


   ///2. muon channel
   ifstream fsys_emscale_m;
   ifstream fsys_emreso_m; 

   double sysem_mean_arr_m[10][20];
   double sysem_reso_arr_m[10][20];

   ///scale
   if(cat>0 && schannel == "mmg"){
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emscale_m.open(TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     while(!fsys_emscale_m.eof()){
       
     string proc;
     int mass;
     double sysem_mean_m;
     
     fsys_emscale_m >> proc >> mass >> sysem_mean_m;

     
     cout<<"reading em reso file "<< TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_mean_arr_m[iproc][5*imass] = sysem_mean_m;
     
     }//while(!datafile.eof())

     for(int iproc=0; iproc<3; iproc++){
       sysem_mean_arr_m[iproc][0] = sysem_mean_arr_m[iproc][5];
       sysem_mean_arr_m[iproc][10] = sysem_mean_arr_m[iproc][5];
     }

     ////reso
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emreso_m.open(TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     while(!fsys_emreso_m.eof()){
       
       string proc;
       int mass;
       double sysem_reso_m;
       
       fsys_emreso_m >> proc >> mass >> sysem_reso_m;
       
     
       cout<<"reading em reso file "<< TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_reso_arr_m[iproc][5*imass] = sysem_reso_m;
     
     }//while(!datafile.eof())

     for(int iproc=0; iproc<3; iproc++){
       sysem_reso_arr_m[iproc][0] = sysem_reso_arr_m[iproc][5];
       sysem_reso_arr_m[iproc][10] = sysem_reso_arr_m[iproc][5];
     }
   
   }//if(cat>0 && schannel == "mmg")

   /////////////////////////////////////SYStematics due to e-pho different/////////////////////////////////////////

   ///1. electron channel
   ifstream fsys_phoscale_e;
   ifstream fsys_phoreso_e; 
   
   double syspho_mean_arr_e[10][20];
   double syspho_reso_arr_e[10][20];
   
   
   ///scale
   if(cat>0 && (schannel == "eeg"||schannel=="eeg_mmg")){

     ///syselepho_hzg_ele_scale_cat1_8TeV_Objpho.txt
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg")  fsys_phoscale_e.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     if(schannel == "eeg_mmg")  fsys_phoscale_e.open(TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoscale_e.eof()){
       
       string proc;
       int mass;
       double syspho_mean_e;
       
       fsys_phoscale_e >> proc >> mass >> syspho_mean_e;
	 
       
       if(schannel == "eeg")  cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_e<<endl;

       if(schannel == "eeg_mmg")  cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_e<<endl;

       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       syspho_mean_arr_e[iproc][5*imass] = syspho_mean_e;
       
     }//while(!datafile.eof())
     
   
       ////reso
       //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg")  fsys_phoreso_e.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     if(schannel == "eeg_mmg")  fsys_phoreso_e.open(TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoreso_e.eof()){
       
       string proc;
       int mass;
       double syspho_reso_e;
       
       fsys_phoreso_e >> proc >> mass >> syspho_reso_e;
       
     
       if(schannel == "eeg")   cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_e<<endl;
       if(schannel == "eeg_mmg")   cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_e<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     syspho_reso_arr_e[iproc][5*imass] = syspho_reso_e;
     
     }//while(!datafile.eof())

   }//if(cat>0 && schannel == "eeg")


   ///2. moun channel
   ifstream fsys_phoscale_m;
   ifstream fsys_phoreso_m; 
   
   double syspho_mean_arr_m[10][20];
   double syspho_reso_arr_m[10][20];
   
   
   ///scale
   if(cat>0 && (schannel == "mmg") ){

     ///syselepho_hzg_ele_scale_cat1_8TeV_Objpho.txt
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     //fsys_phoscale_m.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     fsys_phoscale_m.open(TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoscale_m.eof()){
       
       string proc;
       int mass;
       double syspho_mean_m;
       
       fsys_phoscale_m >> proc >> mass >> syspho_mean_m;
	 
       
       //cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_m<<endl;
       cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_m<<endl;
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       syspho_mean_arr_m[iproc][5*imass] = syspho_mean_m;
       
     }//while(!datafile.eof())
     
   
       ////reso
       //syselepho_hzg_ele_scale_cat1_8TeV.txt
     //fsys_phoreso_m.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     fsys_phoreso_m.open(TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoreso_m.eof()){
       
       string proc;
       int mass;
       double syspho_reso_m;
       
       fsys_phoreso_m >> proc >> mass >> syspho_reso_m;
       
     
       //cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_m<<endl;
       cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     syspho_reso_arr_m[iproc][5*imass] = syspho_reso_m;
     
     }//while(!datafile.eof())

   }//if(cat>0 && schannel == "eeg")






   //////////////////////////////////////////END of systematics due to e-pho difference////////////////////////////////

   /////JEC


   double sys_jec[10][20];
   //if(cat>=1&&cat<=5){
   if(cat>0){
     ifstream fsys_jet; //sys_ele_JEC_8TeV.txt

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_jet.open(TString::Format("sys_%s_JEC_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_jet.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_jet >> proc >> mass >> icat >> sys;
       
       cout<<"reading JEC file "<<TString::Format("sys_%s_JEC_8TeV.txt",schannel.c_str())<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 sys_jec[iproc][5*imass] = sys;
     
     }//while(!datafile.eof())
     
   }


   ////JER
   double sys_jer[10][20];
   //if(cat>=1&&cat<=5){
   if(cat>0){
     ifstream fsys_jet; //sys_ele_JEC_8TeV.txt

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_jet.open(TString::Format("sys_%s_JER_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_jet.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_jet >> proc >> mass >> icat >> sys;
       
       cout<<"reading JER file "<<TString::Format("sys_%s_JER_8TeV.txt",schannel.c_str())<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 sys_jer[iproc][5*imass] = sys;
     
     }//while(!datafile.eof())
   }//if(cat>0)

   //////////////////////////////////END OF JER/////////////////////



   ////////////////////////////////UE/////////////////////

   double sys_ue[10][20];
   if(cat>=0){ ///= after approval
     ifstream fsys_ue; //

     for(int ip=0; ip<7; ip++){ 
       if(cat==0){
	 sys_ue[ip][5] = 0;
	 sys_ue[ip][0] = 0;
	 sys_ue[ip][10] = 0;
       }
     }

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_ue.open(TString::Format("sys_hzg_%s_UE_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_ue.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_ue >> proc >> mass >> icat >> sys;
       
       cout<<"reading UE file "<<TString::Format("sys_hzg_%s_UE_8TeV.txt",name.c_str())<<"... proc mass sysUE "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH") iproc = 4;
       if(proc == "WminusH") iproc = 5;
       if(proc == "WplusH") iproc = 6;
       if(proc == "ZH") iproc = 3;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 {
	   sys_ue[iproc][5*imass] = sys;

	   ///since it is estimated only for 125 GeV
	   sys_ue[iproc][0] = sys;
	   sys_ue[iproc][10] = sys;
	 }
     
       ///after approval
       if(cat==0 && (icat==1||icat==2||icat==3||icat==4||icat==5||icat==10)){

	 sys_ue[iproc][5*imass] += sys/6.;
	 ///since it is estimated only for 125 GeV
	 sys_ue[iproc][0] += sys/6.;
	 sys_ue[iproc][10] += sys/6.;
	 
	 cout<<"UE for cat 0 === cat : icat : sys : tot sys till now "<<cat<<" "<<icat<<" "<<sys/6.<<" "<<sys_ue[iproc][5*imass]<<endl;
       }


     }//while(!datafile.eof())


     for(int ii=0; ii<3; ii++){
       ///Wh = avg(W+H and W-H)


       sys_ue[2][5*ii] = (sys_ue[5][5*ii] + sys_ue[6][5*ii])/2.; 
     }

   }


   for(int ip=0; ip<2; ip++){

     cout<<"UE sys for 125 for cat "<<cat<<" is "<<sys_ue[ip][5]<<endl;
   }

   //////////////////////END OF UE////////////////////////////


   /////////////////////////////////////PS ///////////////////////////
   double sys_ps[10][20];
   if(cat>=0){ //= after approval
     ifstream fsys_ps; //

     for(int ip=0; ip<7; ip++){ 
       if(cat==0){
	 sys_ps[ip][5] = 0;
	 sys_ps[ip][0] = 0;
	 sys_ps[ip][10] = 0;
       }
     }

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_ps.open(TString::Format("sys_hzg_%s_PS_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_ps.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_ps >> proc >> mass >> icat >> sys;

       cout<<"reading PS file "<<TString::Format("sys_hzg_%s_PS_8TeV.txt",name.c_str())<<"... proc mass sysPS "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;       

     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH") iproc = 4;
       if(proc == "WminusH") iproc = 5;
       if(proc == "WplusH") iproc = 6;
       if(proc == "ZH") iproc = 3;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 {
	   sys_ps[iproc][5*imass] = sys;

	   ///since it is estimated only for 125 GeV
	   sys_ps[iproc][0] = sys;
	   sys_ps[iproc][10] = sys;
	 }

       if(cat==0 && (icat==1||icat==2||icat==3||icat==4||icat==5||icat==10)){

	 sys_ps[iproc][5*imass] += sys/6.;
	 ///since it is estimated only for 125 GeV
	 sys_ps[iproc][0] += sys/6.;
	 sys_ps[iproc][10] += sys/6.;
	 
	 cout<<"PS for cat 0 === cat : icat : sys : tot sys till now "<<cat<<" "<<icat<<" "<<sys/6.<<" "<<sys_ps[iproc][5*imass]<<endl;
       }

     
     }//while(!datafile.eof())


   for(int ip=0; ip<2; ip++){

     cout<<"PS sys for 125 for cat "<<cat<<" is "<<sys_ps[ip][5]<<endl;
   }


     for(int ii=0; ii<3; ii++){
       ///Wh = avg(W+H and W-H)


       sys_ps[2][5*ii] = (sys_ps[5][5*ii] + sys_ps[6][5*ii])/2.; 
     }

   }
   ////////////////////////////////////////END OF PS/////////////////

   /////////////////////////////Due to R9/////////////////////////////////

   double sys_r9[10][20];
   if(cat>0){
     ifstream fsys_r9; //

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     //else name = "ele_mu";
     else name = "ele";
     
     fsys_r9.open(TString::Format("sys_%s_R9_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_r9.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_r9 >> proc >> mass >> icat >> sys;
       
       cout<<"reading R9 file "<<TString::Format("sys_%s_R9_8TeV.txt",schannel.c_str())<<"... proc mass sys_xx_R9 "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 sys_r9[iproc][5*imass] = sys;
     
       if(cat>=5) {
	 sys_r9[iproc][5*imass] = 0;
       }

     }//while(!datafile.eof())

   }



   ////////////////////////////////END Of DUE TO R9//////////////////



   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("hzg_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   
   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_13TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");

   //
   // evaluate PDF of background from real data
   //

   //double xmin = 110;
   //double xmin = 100;
   //double xmin = 115;
   //double xmax = 170;

   //TH1D hDataObs("hDataObs", "", 60, xmin, xmax);
   TH1D hDataObs("hDataObs", "", 110, xmin, xmax);
   hDataObs.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");

   fX = new RooRealVar("x", "", xmin, xmax); 
   //fX->setBins(1000);
   //fX->setBins(500);
   fX->setBins(700);


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

   ///SJ
  
   // make background function


   rfit = new RooFitResult();
   
   bgrfit = NULL;
   //poldeg = 4;
   //poldeg = 3;


   //model = "Pow";
   //poldeg = 1;


   model = "Exp";
   poldeg = 1;

   if( (cat>=1 && cat<=4) || cat==10 ) {
     model = "RooPolynomial";
   }
   
   if(cat==1 || cat==3 || cat==4) poldeg = 4;
   if(cat==2) poldeg = 5;
   if(cat==10) poldeg = 3;
   
   if(cat==5 || cat==6789){
   //if(cat==5){
     model = "Pow";
     //model = "Laurent";
     poldeg = 1;
   }

   
   cout<<"Channel cat model "<<channel<<" "<<cat<<" "<<model<<endl;
   /*
   if(cat==6789){
     model = "Exp";
     poldeg = 1;
   }
   */
   /*
   if(cat==6789){
     model = "RooPolynomial";
     //poldeg = 2;
     poldeg = 1;
   }
   */

   /*
   if(cat==10){
     model = "Pow";
     poldeg = 1;
   }
   */



   //RooGaussStepBernstein *bgrfit = new RooGaussStepBernstein();
   //bgrfit = mybkgfit("RooGaussStepBernstein", poldeg);
   mybkgfit(xmin, xmax,cat, schannel);
   //////try importing here//////////////////

   //bkgfit1_modelLaurent_l1_bgr_cat1
   
   //bgrfit1_modelPow_p1
   //bgrfit1_model%s,model.c_str()"%s_p%d"

   
   ///set the parameters constant  - is it the right way? - CHECK
   /*if(model=="RooGaussStepBernstein"){
     mean->setConstant(true);
     sigma->setConstant(true);
     stepval->setConstant(true);
   
   
     if(hDataObs.Integral()!=0){
       for(int ipol=1; ipol<=poldeg; ipol++){
	 
	 pol[ipol]->setConstant(true); 
       }
     }
   }
   */

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
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
     FATAL("RooWorkspace::import() failed");
   

   //SJ
   cout<<"now setting bgrfit pdf to name of pdf "<<endl;
   //bgrfit->GetPar(0)->setConstant(true);
   bgrfit->SetName("pdf");

   cout<<"address of bgrfit "<<bgrfit<<endl;
   
   cout<<"set the name to name of pdf "<<endl;


   if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		     RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		     RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
     cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
   

   ///set the bkg parameter constant
   RooAbsPdf *bgrfit_ws = wspace.pdf(Form("pdf_bgr_cat%d",cat)); 
   RooArgSet *list = bgrfit_ws->getParameters(*fX);
   cout<<"PRINTING LIST"<<endl;
   list->Print();
   
   TIterator* iter = list->createIterator();


   ///set the bkg parameter constant
   for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
     RooRealVar *v = dynamic_cast<RooRealVar *>(a);
     cout<<"Printing my var "<<endl;
     v->Print();
     cout<<"Name of bkg parameter is "<<v->GetName()<<endl;
     //wspace.var(v->GetName())->setConstant(true); 
     if( v->GetName() == "CMS_hzg_mass_13TeV" ) continue;
     if(floatbkg) wspace.var(v->GetName())->setConstant(false); 
     else wspace.var(v->GetName())->setConstant(true); 
   }
   ///set the bkg parameter constant

   
   cout<<"now saved the bkg pdf"<<endl;
   
   
   rfit->SetName(("pdf_bgr_"+catString+"_fitresult"));
  if ( wspace.import(*rfit) )
     FATAL("RooWorkspace::import() fit result failed");



  if(model=="RooGaussStepBernstein"){
    RooArgList bkgpar_list;
    bkgpar_list.add(*mean);
    bkgpar_list.add(*sigma);
    bkgpar_list.add(*stepval);
    bkgpar_list.add(*coeflist);
    RooArgSet bkgpar_set(bkgpar_list);
    
    
    cout<<"Values of parameters after fixing ..."<<endl;
    cout<<"Mean : sigma : stepval : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<stepval->getVal()<<endl;
  }
    
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
    
    
    double sigma = 1.;
    double integxmin = 125.0-2*sigma;
    double integxmax = 125.0+2*sigma;
    
    integxmin = 115;

    integxmax = 170;
    
    integxmin = 125.0-2*sigma;
    integxmax = 125.0+2*sigma;
    
    
    //fX->setRange("2sigma",xmin,xmax) ;
    fX->setRange("2sigma",integxmin,integxmax) ;
    RooAbsReal* intbkgfit12s = bgrfit->createIntegral(*fX,NormSet(*fX),Range("2sigma")); 
    
    int ibin = hDataObs.FindBin(integxmin);
    int jbin = hDataObs.FindBin(integxmax)-1;
    int observation2s = TMath::Nint(hDataObs.Integral(ibin, jbin));

    //cout<<"Checks on createIntegral -- Obs : Exp : "<<observation2s<<" "<<intbkgfit12s->getVal()*observation<<endl;
    cout<<"Checks on createIntegral -- Obs : Exp : "<<observation2s<<" "<<intbkgfit12s->getVal()*expectation<<endl;

   RooRealVar bkgNorm("norm","",expectation);
   bkgNorm.setConstant(false);

   bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
   if ( wspace.import(bkgNorm) )
     FATAL("RooWorkspace::import() failed");
   



  /////start with the signal now /////////

   cout<<"DOING SIGNAL FIT NOW"<<endl;
   double    expected[nprocess][nmass];
   //const char* proc[nprocess] = {"ggH", "VBF", "WplusH", "WminusH", "ZH", "ttH"};
   const char* proc[nprocess] = {"ggH", "VBF", "WH", "ZH", "ttH"};
   //const char* proc[nprocess] = {"ggH", "VBF", "ZH", "WminusH", "ttH"};
   int mass_exist[nsig] = {120, 125, 130};
   double sigeff_All[nprocess][nmass];
   double sigNraw[nprocess][nmass];
   double sigeff_ggf_vbf[nmass];
   double sigsys_ggf_vbf[nmass];

   double sigsys_ggf_vbf_theoryalso[nmass];
   double sigsys_ggf_vbf_theoryonly[nmass];

   ///systematcs ID
   double sigSysonExp_lep[nprocess][nmass];
   double sigSysonExp_pho[nprocess][nmass];

   double sigSysonEff_lep[nprocess][nmass];
   double sigSysonEff_pho[nprocess][nmass];

   double sigTotSys[nprocess][nmass];
   double sigSysPU[nprocess][nmass];
   
   double sigTothltSys[nprocess][nmass]; 

   double sigSysEff[nprocess][nmass];

   double sigSysEff_theoryalso[nprocess][nmass];
   double sigSysEff_theoryonly[nprocess][nmass];

   for (int p = 0; p < nprocess; p++){
     if(p>nprocess) break;
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
       
       //TH1D hMass(hname, "", 120, xmin,xmax);
       TH1D hMass(hname, "", nBinsForMass, xmin,xmax);
       hMass.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
       hMass.SetYTitle("Counts");
       hMass.Sumw2();
       
       //hMasssig[p][5*m] = new TH1D(hname, "", 120, xmin,xmax);
       hMasssig[p][5*m] = new TH1D(hname, "", nBinsForMass, xmin,xmax);

       RooAbsData* rooMass;
       
       TString filename;

       double trigEff = 1;
       
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
	 
	 if( schannel == "eeg" ) trigEff = trigEff_ele;
	 if( schannel == "mmg" ) trigEff = trigEff_mu;
	 
	 cout<<"Calling for fill_events for isgn, trigEff is "<<trigEff<<endl;
	 rooMass = fill_events(&hMass, filename, trigEff, cat);
	 ///for WH, combine both W+H and W-H
	
	 
       }
     
     

       TString filename1;
       TString filename2;
       TString filename3;
       TString filename4;
       
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

       *hMasssig[p][5*m] = hMass;
       
       cout<<"=================Channel, cat "<<channel<<" "<<cat<<endl;
       cout<<"===========fitting for mass =================="<<mass<<endl;
       mysigfit(sigmodel,xmin, xmax, _nuisance, p, 5*m, cat, schannel);
       
       setSigParConst(sigmodel, _nuisance);
       
       ////get expectation from ZH and ttH
       /*RooAbsData* rooMassZH = fill_events(&hMassZH, filename_ZH, trigEff, cat);
       RooAbsData* rooMassttH = fill_events(&hMassttH, filename_ttH, trigEff, cat);
       */

       ///set the parameters const at this point for these mass points
       
    
       // save the normalization factor which adapts the MC scale to real
       // data; evaluate expected signal yield
       double norm = normfactor_SMHiggs_8TeV(channel, p, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double norm = 1;
       
     
       ///SJ
       double sigEff = 1;

       if(cat!=6789){
	 sigEff = hMass.Integral()/ngen(channel,p,mass);
       }

       if(cat==6789){
	 sigEff = hMass.Integral()/(ngen("eeg",p,mass)+ngen("mmg",p,mass)); 
       }



       if(p>1){
	 ///UE
	 sys_ue[p][5*m] = sys_ue[0][5*m];
	 if(cat==5)
	   sys_ue[p][5*m] = sys_ue[0][5*m];
	 ///PS
	 sys_ps[p][5*m] = sys_ps[0][5*m];
	 if(cat==5)
	   sys_ps[p][5*m] = sys_ps[0][5*m];
	 
       }
       
       //double sigExp = sigEff * lumi; /// use this for limit on xsec*BR
       double sigExp = hMass.Integral()*norm;///use this for limit on xsec/SMxsec
       
       ///systematics on ID
       double sys_lep = 0;
       double sys_pho = 0;
       
       double sys_hlt = 0;
       
       if(cat!=6789){
	 getIDsys(filename,  trigEff, sys_lep, sys_pho, cat, true);
	 getHLTsys(filename,  trigEff, sys_hlt, cat, true);
       }

       if(cat==6789){
	 getIDsys(filename1, filename2, trigEff_ele, trigEff_mu, sys_lep, sys_pho, cat, true);
	 getHLTsys(filename1, filename2, trigEff_ele, trigEff_mu, sys_hlt, cat, true);
       }
       
       double sys_pu = 0;
       if(cat!=6789) getPUsys(filename,  trigEff, sys_pu, cat, true);
       if(cat==6789) getPUsys(filename1, filename2, trigEff_ele, trigEff_mu, sys_pu, cat, true);
       
       //sigSysonExp_lep[p][m*5] = sqrt(sys_lep)*norm/sigExp;
       //sigSysonExp_pho[p][m*5] = sqrt(sys_pho)*norm/sigExp;
       
       sigSysPU[p][m*5] = sys_pu/hMass.Integral();
       
       sigSysonExp_lep[p][m*5] = sqrt(sys_lep)/hMass.Integral();
       sigSysonExp_pho[p][m*5] = sqrt(sys_pho)/hMass.Integral();
       


       sigTotSys[p][m*5] = sqrt( pow(sigSysonExp_lep[p][m*5],2) + pow(sigSysonExp_pho[p][m*5],2) );
       sigTothltSys[p][m*5] = sqrt(sys_hlt)/hMass.Integral();

       cout<<"==ID sys==="<<endl;
       cout<<"Mass cat p "<<m*5 << " "<<cat <<" "<<p<<" "<<endl;
       cout<<"sys_lep : sys_pho : % of sys_lep : % of sys_pho : % tot : hmass : "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<sigSysonExp_lep[p][m*5]<< " "<<sigSysonExp_pho[p][m*5]<<" "<<sigTotSys[p][m*5]<< " "<<hMass.Integral()<<endl;
       cout<<"===HLT sys - "<<sigTothltSys[p][m*5]<<endl;

       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) );

       //includes lumi
       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) +pow(2.7*hMass.Integral()/100.,2) );

       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) );


       ///idsys gives square  - Till 7th Feb, 2018 - for approval
       //sigSysEff[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) );

       ///include HLT eff and UE/PS - after approval - 7th Feb, 2018
       cout<<"after approval numbers: sys_lep : sys_pho : sus_pu : sys_hlt : UE : PS "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<sys_pu<<" "<<sqrt(sys_hlt)<<" " << sys_ue[p][m*5]*hMass.Integral() << " "<< sys_ps[p][m*5]*hMass.Integral()<<endl;
       cout<<"individual numbers : UE : PS : integ : "<<sys_ue[p][m*5]<<" "<<sys_ps[p][m*5]<<" "<<hMass.Integral()<<endl;

       sigSysEff[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) + sys_hlt );


       ///sys_ue and sys_ps are estimated as: change in yield/total yield so multiply by total yield ot get the correct number
       sigSysEff_theoryalso[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) + sys_hlt + pow(sys_ue[p][m*5]*hMass.Integral(),2) + pow(sys_ps[p][m*5]*hMass.Integral(),2) );

       sigSysEff_theoryonly[p][m*5] = sqrt( pow(sys_ue[p][m*5]*hMass.Integral(),2) + pow(sys_ps[p][m*5]*hMass.Integral(),2) );
     

       ///sys on Eff is not always right. Only when sigExp=sigEff*lumi, then it is right. 
       ///Only the above sigSysonExp_lep numbers are right. These are on Eff   as well. So sigSysonExp_lep can be used on Exp and Eff
       if(cat!=6789){
	 sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen(channel,p,mass)*sigExp);
	 sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen(channel,p,mass)*sigExp);
       }

       if(cat==6789){

	 sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen("eeg",p,mass)*sigExp+ngen("mmg",p,mass)*sigExp);
	 sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen("eeg",p,mass)*sigExp+ngen("mmg",p,mass)*sigExp);
       }
       
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
       wspace.import(*rooMassSig[p][m*5],RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV"));
       
       


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
			 RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
	 FATAL("RooWorkspace::import() failed");
       
       
       cout<<"saved the signal pdf"<<endl;

       rsig[p][m*5]->SetName(Form("fitresult_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*rsig[p][m*5]))
	 FATAL("RooWorkspace::import() failed");
       
       cout<<"saved the signal fit result"<<endl;



       ////setting the n, alpha, siggaus fixed 

       RooAbsPdf *sigfit_ws = wspace.pdf(Form("pdf_sig_%s_%d_cat%d",proc[p],mass,cat) );
       RooArgSet *list_sig = sigfit_ws->getParameters(*fX);
       cout<<"PRINTING SIG LIST"<<endl;
       list_sig->Print();
   
       TIterator* iter_sig = list_sig->createIterator();

       /////Signal parameter constant
       for (RooAbsArg *a = (RooAbsArg *) iter_sig->Next(); a != 0; a = (RooAbsArg *) iter_sig->Next()) {
	 RooRealVar *v = dynamic_cast<RooRealVar *>(a);
	 cout<<"Printing my var "<<endl;
	 v->Print();
	 cout<<"Name of sig parameter is "<<v->GetName()<<endl;

	 string name = v->GetName();

	 int fstr = name.find("alpha");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }


	 fstr = name.find("frac");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }

	 fstr = name.find("power");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }

	 fstr = name.find("s21");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }


       }
       ////END of setting the n, alpha, siggaus fixed 



       /////quality of fit////
       cout<<"process is " << proc[p] <<" for cat "<<cat<<" channel "<<schannel<<" and mass "<<mass<<" quality of fit (NLL) "<<rsig[p][m*5]->minNll()<<endl;
       
       /////for systematics on resolution and scale
     
        if(schannel=="mmg"){
	 
	 double mean = sigmean->getVal();
	 double sigma = sigsigma->getVal();
	
 	 //for now taking all the lepton and photon energy/mom related sys for other samples (WH, Zh and ttH) from ggF
	 if(p>2){
	   /*sysrochor_mean_arr[p][5*m] = sysrochor_mean_arr[0][5*m];
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[0][5*m];
	   sysem_mean_arr_m[p][5*m] = sysem_mean_arr_m[0][5*m];
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[0][5*m];

	   syspho_mean_arr_m[p][5*m] = syspho_mean_arr_m[0][5*m];
	   syspho_reso_arr_m[p][5*m] = syspho_reso_arr_m[0][5*m];
	   */
	   sysrochor_mean_arr[p][5*m] = sysrochor_mean_arr[2][5*m];
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[2][5*m];
	   sysem_mean_arr_m[p][5*m] = sysem_mean_arr_m[2][5*m];
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[2][5*m];

	   syspho_mean_arr_m[p][5*m] = syspho_mean_arr_m[2][5*m];
	   syspho_reso_arr_m[p][5*m] = syspho_reso_arr_m[2][5*m];
	   

	 }

	 
	 cout<<"Mass process cat channel "<<5*m<<" "<<p<<" "<<cat<<" "<<channel<<endl;
	 cout<<"sysrochor_reso_arr sysem_reso_arr_m sigma  "<<sysrochor_reso_arr[p][5*m]<<" "<<sysem_reso_arr_m[p][5*m]<<" "<<sigma<<endl;
	 cout<<"(%) sysrochor_reso_arr sysem_reso_arr_m sigma  "<<sysrochor_reso_arr[p][5*m]/sigma<<" "<<sysem_reso_arr_m[p][5*m]/sigma<<" "<<sigma<<endl;
	 cout<<"(%) sysrochor_mean_arr sysem_mean_arr_m sigma  "<<sysrochor_mean_arr[p][5*m]/mean<<" "<<sysem_mean_arr_m[p][5*m]/mean<<" "<<sigma<<endl;


	 ///because the other samples already have added values
	 if(p<=2){
	   sysrochor_mean_arr[p][5*m] =  sysrochor_mean_arr[p][5*m]/mean;
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[p][5*m]/sigma;
	   
	   ///1 a;ready added above so no need to add here - see mail from Louie: 6th march: multiply the resolution by (1+nuisance1+nuisance2+...).
	   sysem_mean_arr_m[p][5*m] =  sysem_mean_arr_m[p][5*m]/mean;
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[p][5*m]/sigma;

	   syspho_mean_arr_m[p][5*m] =  syspho_mean_arr_m[p][5*m]/mean;
	   syspho_reso_arr_m[p][5*m] = syspho_reso_arr_m[p][5*m]/sigma;
	 
	   sysrochor_mean_arr[p][5*m] =  sqrt( pow(sysrochor_mean_arr[p][5*m],2) + pow(sysem_mean_arr_m[p][5*m],2) + pow(syspho_mean_arr_m[p][5*m],2) );
	   sysrochor_reso_arr[p][5*m] =  sqrt( pow(sysrochor_reso_arr[p][5*m],2) + pow(sysem_reso_arr_m[p][5*m],2) + pow(syspho_reso_arr_m[p][5*m],2) );
	   
	 }


	 
	 cout<<"final sysrochor_mean_arr is "<<sysrochor_mean_arr[p][5*m]<<endl;
	 cout<<"final sysrochor_reso_arr is "<<sysrochor_reso_arr[p][5*m]<<endl;

	 cout<<"final sysem_mean_arr is "<<sysem_mean_arr_m[p][5*m]<<endl;
	 cout<<"final sysem_reso_arr is "<<sysem_reso_arr_m[p][5*m]<<endl;

	 cout<<"final syspho_mean_arr is "<<syspho_mean_arr_m[p][5*m]<<endl;
	 cout<<"final syspho_reso_arr is "<<syspho_reso_arr_m[p][5*m]<<endl;

	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d[0]",p,5*m,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod

	 wspace.factory(Form("CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d[0]",p,5*m,cat));

	 //wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d , CMS_hzg_delta_muonEM_mean_chan%d_m%d_cat%d , CMS_hzg_delta_muonPho_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat, p,5*m,cat,p,5*m,cat));
	 //wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_muonEM_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_muonPho_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat, p,5*m,cat,p,5*m,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat));

	 

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );
	 
       }
     
	if(p>2){
	  sys_jec[p][5*m] = sys_jec[2][5*m];
	}

	/*	////////since UE/PS is not god for VH samples. I use it from ggF/ VBF
	if(p>1){
	  ///UE
	  sys_ue[p][5*m] = sys_ue[0][5*m];
	  if(cat==5)
	    sys_ue[p][5*m] = sys_ue[0][5*m];
	  ///PS
	  sys_ps[p][5*m] = sys_ps[0][5*m];
	  if(cat==5)
	    sys_ps[p][5*m] = sys_ps[0][5*m];

	}
	*/

	///em scale 
       if(schannel=="eeg" || schannel=="eeg_mmg"){
	 
	 double mean = sigmean->getVal();
	 double sigma = sigsigma->getVal();

	 //for now taking all the lepton and photon energy/mom related sys for other samples (WH, Zh and ttH) from ggF
	 if(p>2){
	   /*sysem_mean_arr_e[p][5*m] = sysem_mean_arr_e[0][5*m];
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[0][5*m];

	   syspho_mean_arr_e[p][5*m] = syspho_mean_arr_e[0][5*m];
	   syspho_reso_arr_e[p][5*m] = syspho_reso_arr_e[0][5*m];
	   */
	   sysem_mean_arr_e[p][5*m] = sysem_mean_arr_e[2][5*m];
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[2][5*m];

	   syspho_mean_arr_e[p][5*m] = syspho_mean_arr_e[2][5*m];
	   syspho_reso_arr_e[p][5*m] = syspho_reso_arr_e[2][5*m];
	   

	 }


	 cout<<"Mass process cat channel "<<5*m<<" "<<p<<" "<<cat<<" "<<channel<<endl;
	 cout<<"sysem_reso_arr_e sigma  "<<sysem_reso_arr_e[p][5*m]<<" "<<" "<<sigma<<endl;
	 cout<<"(%) sysem_reso_arr_e sigma  "<<sysem_reso_arr_e[p][5*m]/sigma<<" "<<" "<<sigma<<endl;
	 cout<<"(%) sysem_mean_arr_e sigma  "<<sysem_mean_arr_e[p][5*m]/mean<<" "<<" "<<mean<<endl;
	 
	 if(p<=2){
	   sysem_mean_arr_e[p][5*m] =  sysem_mean_arr_e[p][5*m]/mean;
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[p][5*m]/sigma;

	   syspho_mean_arr_e[p][5*m] =  syspho_mean_arr_e[p][5*m]/mean;
	   syspho_reso_arr_e[p][5*m] = syspho_reso_arr_e[p][5*m]/sigma;


	   sysem_mean_arr_e[p][5*m] =  sqrt( pow(sysem_mean_arr_e[p][5*m] ,2) + pow(syspho_mean_arr_e[p][5*m],2) );
	   sysem_reso_arr_e[p][5*m] =  sqrt( pow(sysem_reso_arr_e[p][5*m] ,2) + pow(syspho_reso_arr_e[p][5*m],2) );

	   
	 }


	 cout<<"final sysem_mean_arr is "<<sysem_mean_arr_e[p][5*m]<<endl;
	 cout<<"final sysem_reso_arr_e is "<<sysem_reso_arr_e[p][5*m]<<endl;

	 cout<<"final syspho_mean_arr is "<<syspho_mean_arr_e[p][5*m]<<endl;
	 cout<<"final syspho_reso_arr is "<<syspho_reso_arr_e[p][5*m]<<endl;


	 

	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d[0]",p,5*m,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod

	 wspace.factory(Form("CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d[0]",p,5*m,cat));

	 /*wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_ele_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d, CMS_hzg_delta_elePho_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat, p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_ele_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_elePho_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat,  p,5*m,cat, p,5*m,cat));
	 
	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );
	 */


	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_ele_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_ele_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat,  p,5*m,cat));
	 
	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );

	 
       }

       /////for systematics on resolution and scale


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
     double systot_theoryalso = 0;
     double systot_theoryonly = 0;
     
     for(int ip=0; ip<nprocess; ip++){
       
       double norm = normfactor_SMHiggs_8TeV(channel, ip, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double ngenTot = ngen_Tot(channel, ip, mass); ///total no of ngen e vents = ee mm tt
       if(cat!=6789) ntotgen += ngen(channel,ip,mass)*norm;
       
       if(cat==6789) ntotgen += ngen("eeg",ip,mass)*norm + ngen("mmg",ip,mass)*norm;

       numtot  += sigNraw[ip][m * 5]*norm;
       //systot += pow(sigNraw[ip][m * 5]*norm,2); ///check
       cout<<"channel "<<schannel<<" category=="<<cat<<"=====for process "<<ip<<" mass "<<mass<<", numtot "<<numtot<<" and ntotgen "<<ntotgen<<endl;

       systot += pow(sigSysEff[ip][m * 5]*norm,2); ///check
       systot_theoryalso += pow(sigSysEff_theoryalso[ip][m * 5]*norm,2); ///check
       systot_theoryonly += pow(sigSysEff_theoryonly[ip][m * 5]*norm,2); ///check

       cout<<" wo theory : theory : all "<<sigSysEff[ip][m * 5]<<" "<<sigSysEff_theoryonly[ip][m * 5]<<endl;
       cout<<"Till now wo theory : theor : all "<<systot<<" "<<systot_theoryonly<<" "<<systot_theoryalso<<endl;
       cout<<"total "<<ntotgen/norm<<endl;
     }
     double sigEff = numtot/ntotgen;

     sigeff_ggf_vbf[m*5] = sigEff;

     //double sigTotsys_ggf_vbf = sqrt( pow(sigTotSys[0][m * 5],2) + pow(sigTotSys[1][m * 5],2) ); 
     double sigTotsys_ggf_vbf = sqrt( systot )/ntotgen;
     sigsys_ggf_vbf[m*5] = sigTotsys_ggf_vbf;

     sigsys_ggf_vbf_theoryalso[m*5] = sqrt( systot_theoryalso )/ntotgen;
     sigsys_ggf_vbf_theoryonly[m*5] = sqrt( systot_theoryonly )/ntotgen;

     cout<<"After everything, eff : wo theory thoery only all  "<<sigEff<<" "<<sigsys_ggf_vbf_theoryonly[m*5]<<" "<<sigsys_ggf_vbf_theoryalso[m*5]<<endl;


   }


  
   ///interpolation
   // use 1 GeV steps
   for (int p = 0; p < nprocess; p++){ // 2 processes
   //for (int p = nprocess-1; p >=0; p--){ // 2 processes
     for (int mm = 0; mm< nsig-1; mm++){ 
       for (int k = 1; k <= 4; k++) {
	 
	 int mass = 120 + mm * 5 + k;

	 /*
	 rooMassSig[p][mm * 5 + k] = NULL;
	 sigfit[p][mm * 5 + k] = NULL;
	 rsig[p][mm * 5 + k] = NULL;
	 */
       
	 cout<<""<<endl;
	 cout<<"cat : p : mm : k : mass : channel : "<<cat<<" "<<p<<" "<<mm<<" "<<k<<" "<<mass<<" "<<schannel<<endl;
	 
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

	 sigsys_ggf_vbf[mm * 5 + k] = sigsys_ggf_vbf[mm*5] * a + b * sigsys_ggf_vbf[5*(mm+1)];

	 sigsys_ggf_vbf_theoryalso[mm * 5 + k] = sigsys_ggf_vbf_theoryalso[mm*5] * a + b * sigsys_ggf_vbf_theoryalso[5*(mm+1)];
	 sigsys_ggf_vbf_theoryonly[mm * 5 + k] = sigsys_ggf_vbf_theoryonly[mm*5] * a + b * sigsys_ggf_vbf_theoryonly[5*(mm+1)];


	 sigTotSys[p][mm * 5 + k] = sigTotSys[p][mm*5]* a + sigTotSys[p][(mm+1)*5] * b;
       
	 ///
	 sigSysonExp_lep[p][mm * 5 + k] = sigSysonExp_lep[p][mm*5]* a + sigSysonExp_lep[p][(mm+1)*5] * b;
	 sigSysonExp_pho[p][mm * 5 + k] = sigSysonExp_pho[p][mm*5]* a + sigSysonExp_pho[p][(mm+1)*5] * b;
	 

	 sigSysPU[p][mm * 5 + k] = sigSysPU[p][mm*5]* a + sigSysPU[p][(mm+1)*5] * b;

	 sigTothltSys[p][mm * 5 + k] = sigTothltSys[p][mm*5]* a + sigTothltSys[p][(mm+1)*5] * b;

	 
	 //sigsys_ggf_vbf[mm * 5 + k] = sigsys_ggf_vbf[mm*5] * a + b * sigsys_ggf_vbf[5*(mm+1)]; 

	 sigSysonEff_lep[p][mm * 5 + k] = sigSysonEff_lep[p][mm*5] * a + b * sigSysonEff_lep[p][5*(mm+1)];
	 sigSysonEff_pho[p][mm * 5 + k] = sigSysonEff_pho[p][mm*5] * a + b * sigSysonEff_pho[p][5*(mm+1)];
	 
	 ///interpolate all the sys
	 sys_jec[p][mm * 5 + k] = sys_jec[p][mm*5]* a + sys_jec[p][(mm+1)*5] * b;

	 sys_ue[p][mm * 5 + k] = sys_ue[p][mm*5]* a + sys_ue[p][(mm+1)*5] * b;
	 sys_ps[p][mm * 5 + k] = sys_ps[p][mm*5]* a + sys_ps[p][(mm+1)*5] * b;
	 
	 sysrochor_mean_arr[p][mm * 5 + k] = sysrochor_mean_arr[p][mm*5]* a + sysrochor_mean_arr[p][(mm+1)*5] * b;
	 sysem_mean_arr_m[p][mm * 5 + k] = sysem_mean_arr_m[p][mm*5]* a + sysem_mean_arr_m[p][(mm+1)*5] * b;
	 syspho_mean_arr_m[p][mm * 5 + k] = syspho_mean_arr_m[p][mm*5]* a + syspho_mean_arr_m[p][(mm+1)*5] * b;

	 sysrochor_reso_arr[p][mm * 5 + k] = sysrochor_reso_arr[p][mm*5]* a + sysrochor_reso_arr[p][(mm+1)*5] * b;
	 sysem_reso_arr_m[p][mm * 5 + k] = sysem_reso_arr_m[p][mm*5]* a + sysem_reso_arr_m[p][(mm+1)*5] * b;
	 syspho_reso_arr_m[p][mm * 5 + k] = syspho_reso_arr_m[p][mm*5]* a + syspho_reso_arr_m[p][(mm+1)*5] * b;

	 sysem_mean_arr_e[p][mm * 5 + k] = sysem_mean_arr_e[p][mm*5]* a + sysem_mean_arr_e[p][(mm+1)*5] * b;
	 syspho_mean_arr_e[p][mm * 5 + k] = syspho_mean_arr_e[p][mm*5]* a + syspho_mean_arr_e[p][(mm+1)*5] * b;

	 sysem_reso_arr_e[p][mm * 5 + k] = sysem_reso_arr_e[p][mm*5]* a + sysem_reso_arr_e[p][(mm+1)*5] * b;
	 syspho_reso_arr_e[p][mm * 5 + k] = syspho_reso_arr_e[p][mm*5]* a + syspho_reso_arr_e[p][(mm+1)*5] * b;
	 

	 


	 cout<<"========starting to interpolate now for mass "<<mass<<"================="<<endl;
	 cout<<"====masses taken are "<<m1<<" and "<<m2<<endl;
	 siginterpolate(sigmodel, p, mm*5, k, a,b, _nuisance, cat);
	 cout<<"BACK after interpolating"<<endl;
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
			   RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
	   FATAL("RooWorkspace::import() failed");
	 

       if(schannel=="mmg"){
	 
	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));
	 
	 wspace.factory(Form("CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));


	 /*wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d, CMS_hzg_delta_muonEM_mean_chan%d_m%d_cat%d, CMS_hzg_delta_muonPho_mean_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat, p,5*mm+k,cat,p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_muonEM_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_muonPho_sigma_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat, p,5*mm+k,cat, p,5*mm+k,cat));

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );

	 */

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat));

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 
       }

	///em scale 
       if(schannel=="eeg"  || schannel=="eeg_mmg"){
	 
	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d[0]",p,5*mm+k,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod

	 wspace.factory(Form("CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));

	 
	 /*wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_ele_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d, CMS_hzg_delta_elePho_mean_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat, p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_ele_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_elePho_sigma_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat, p,5*mm+k,cat));

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 */


	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_ele_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_ele_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat));

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 
       }


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





   

   const char* tmpproc[nprocess] = {"ggH", "qqH", "WH", "ZH", "ttH"};
   for(int ip=0; ip<nprocess; ip++){
     cout<<"Writing hte work space now"<<endl;
     
     wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_%s_13TeV.root",
					tmpproc[ip], channel, cats));
     
   }//for(int ip=0; ip<nprocess; ip++)


   wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_13TeV.root",
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
   TString binString = TString::Format("%s_%s_13TeV", channel, cats);
   const char* bin = binString.Data();

   //for (int m = 0; m < 31; m++) {
   for (int m = 0; m < 11; m++) {
      int mass = 120 + m;
      

      /*
      sigSysonExp_lep[0][m] = 1 + sigSysonExp_lep[0][m];
      sigSysonExp_pho[0][m] = 1 + sigSysonExp_pho[0][m];

      sigSysonExp_lep[1][m] = 1 + sigSysonExp_lep[1][m];
      sigSysonExp_pho[1][m] = 1 + sigSysonExp_pho[1][m];
      */

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
      //datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      //datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      
      //when energy/mom scale sys is applied  
      //datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      //datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);

      /*datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ZH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_WminusH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_WplusH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ttH_%d_%s\n", channel, cats, mass,cats);
      */


      
      //21 jan, 2017 - use ggF shape for other processes for cats 1-4 since these events yield is quite small
      if(cat>=0&&cat<=4){

	datacard += TString::Format("shapes  qqH_hzg         * for_datacards_hzg_qqH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ggH_hzg         * for_datacards_hzg_ggH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);

      
	datacard += TString::Format("shapes  ZH_hzg         * for_datacards_hzg_ZH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WH_hzg         * for_datacards_hzg_WH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH_hzg         * for_datacards_hzg_ttH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      }


      //21 jan, 2017 - use VBF shape for other processes for cats 5 since these events yield is quite small
      if(cat==5){

	datacard += TString::Format("shapes  qqH_hzg         * for_datacards_hzg_qqH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ggH_hzg         * for_datacards_hzg_ggH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);


	datacard += TString::Format("shapes  ZH_hzg         * for_datacards_hzg_ZH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WH_hzg         * for_datacards_hzg_WH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH_hzg         * for_datacards_hzg_ttH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      }
      

      if(cat==6789){
	datacard += TString::Format("shapes  ZH_hzg         * for_datacards_hzg_ZH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ZH_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_WminusH_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_WplusH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WH_hzg         * for_datacards_hzg_WH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_WH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH_hzg         * for_datacards_hzg_ttH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ttH_%d_%s\n", channel, cats, mass,cats);

	datacard += TString::Format("shapes  qqH_hzg         * for_datacards_hzg_qqH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ttH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ggH_hzg         * for_datacards_hzg_ggH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ZH_%d_%s\n", channel, cats, mass,cats);
	
      }


      if(cat==10){
	datacard += TString::Format("shapes  ZH_hzg         * for_datacards_hzg_ZH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	//datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WH_hzg         * for_datacards_hzg_WH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH_hzg         * for_datacards_hzg_ttH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);

	datacard += TString::Format("shapes  qqH_hzg         * for_datacards_hzg_qqH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ggH_hzg         * for_datacards_hzg_ggH_%s_%s_13TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	
      }
      

      /*
      //21 jan, 2017 - use ggF shape for other processes for cats 1-4 since these events yield is quite small
      if(cat>=1&&cat<=4){
	datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      }


      //21 jan, 2017 - use VBF shape for other processes for cats 5 since these events yield is quite small
      if(cat==5){
	datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
	datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      }
      */
      
      
      
      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:pdf_bgr_%s\n", channel, cats,cats);
      datacard += TString::Format("shapes  data_obs  * for_datacards_hzg_%s_%s_13TeV.root hzg_workspace:data_obs_%s\n", channel, cats,cats);
      
      datacard += "------------------------------------------------------------------------------------------------------------\n";


      //datacard += TString::Format("bin            %s\n", bin);
      datacard += TString::Format("bin            %s\n", cats);
      //datacard += TString::Format("observation    %d\n", observation);
      datacard += TString::Format("observation    %d\n", -1);
      cout<<"wrote about teh observation ... "<<endl;
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %-15s %-15s %-15s\n", bin, bin, bin);
      /*datacard += TString::Format("bin            %-15s %-15s %-15s\n", cats, cats, cats);
      datacard +=                 "process        VBFH            ggH             bgr\n";
      datacard +=                 "process        -1              0               1\n";
      datacard += TString::Format("rate           %-15f %-15f %-15f\n",
				  //			  expected[1][m], expected[0][m], expectation);
				  expected[1][m], expected[0][m], 1.0);
      */

      
      cout<<"writing the process expectation"<<endl;
      cout<<"ggF "<<expected[0][m]<<endl;
      cout<<"VBF "<<expected[1][m]<<endl;
      cout<<"WH "<<expected[2][m]<<endl;
      cout<<"ZH "<<expected[3][m]<<endl;
      cout<<"ttH "<<expected[4][m]<<endl;
      //cout<<"W+H "<<expected[2][m]<<endl;
      //cout<<"W-H "<<expected[3][m]<<endl;
      //cout<<"ZH "<<expected[4][m]<<endl;
      //cout<<"ttH "<<expected[5][m]<<endl;
      
      /*
      datacard += TString::Format("bin            %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n", cats, cats, cats,cats, cats,cats, cats);
      datacard +=                 "process        ttH       ZH      WminusH        WplusH           VBFH            ggH             bgr\n";
      datacard +=                 "process        -5        -4       -3             -2                -1             0               1\n";

      cout<<"all the expectations ... "<<expected[0][m]<<" "<<expected[1][m]<<" "<<expected[2][m]<<" "<<expected[3][m]<<" "<<expected[4][m]<<" "<<expected[5][m]<<endl;
      
      datacard += TString::Format("rate           %-15f %-15f %-15f %-15f %-15f %-15f %-15f\n",
				  expected[5][m], expected[4][m], expected[3][m], expected[2][m], expected[1][m], expected[0][m], 1.0);
      */


      datacard += TString::Format("bin            %-15s %-15s %-15s %-15s %-15s %-15s\n", cats, cats, cats,cats, cats,cats);
      datacard +=                 "process        ttH_hzg       ZH_hzg        WH_hzg           qqH_hzg            ggH_hzg             bgr\n";
      datacard +=                 "process       -4       -3             -2                -1             0               1\n";

      //cout<<"all the expectations ... "<<expected[0][m]<<" "<<expected[1][m]<<" "<<expected[2][m]<<" "<<expected[3][m]<<" "<<expected[4][m]<<" "<<expected[5][m]<<endl;
      cout<<"all the expectations ... "<<expected[0][m]<<" "<<expected[1][m]<<" "<<expected[2][m]<<" "<<expected[3][m]<<" "<<expected[4][m]<<endl;
      
      datacard += TString::Format("rate           %-15f %-15f %-15f %-15f %-15f %-15f\n",
				  expected[4][m], expected[3][m], expected[2][m],  expected[1][m], expected[0][m], 1.0);

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

   
      // CMS uncertainties.
      //
      // NOTE: naming conventions are from
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
      datacard += "lumi_13TeV          lnN     1.025          1.025   1.025    1.025   1.025                -\n";
      
      if(schannel=="mmg"){
	datacard += TString::Format("CMS_HLTeff_m_13TeV     lnN        %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1.013,1.013,1.013,1.013,1.013);
      }
      else{
	datacard += TString::Format("CMS_HLTeff_e_13TeV     lnN        %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigTothltSys[4][m], 1+sigTothltSys[3][m], 1+sigTothltSys[2][m], 1+sigTothltSys[1][m], 1+sigTothltSys[0][m]);
      }
      
      //datacard += TString::Format("CMS_IDeff_13TeV     lnN        %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigTotSys[4][m], 1+sigTotSys[3][m], 1+sigTotSys[2][m], 1+sigTotSys[1][m], 1+sigTotSys[0][m]);

      if(schannel=="mmg"){
	datacard += TString::Format("CMS_IDeff_m_13TeV     lnN        %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysonExp_lep[4][m], 1+sigSysonExp_lep[3][m], 1+sigSysonExp_lep[2][m], 1+sigSysonExp_lep[1][m], 1+sigSysonExp_lep[0][m]);
      }
      else{
	datacard += TString::Format("CMS_IDeff_e_13TeV     lnN        %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysonExp_lep[4][m], 1+sigSysonExp_lep[3][m], 1+sigSysonExp_lep[2][m], 1+sigSysonExp_lep[1][m], 1+sigSysonExp_lep[0][m]);
      }

      datacard += TString::Format("CMS_IDeff_g_13TeV     lnN        %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysonExp_pho[4][m], 1+sigSysonExp_pho[3][m], 1+sigSysonExp_pho[2][m], 1+sigSysonExp_pho[1][m], 1+sigSysonExp_pho[0][m]);

   
      //datacard += TString::Format("CMS_PU_13TeV     lnN        %-15f          %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysPU[5][m], 1+sigSysPU[4][m], 1+sigSysPU[3][m], 1+sigSysPU[2][m], 1+sigSysPU[1][m], 1+sigSysPU[0][m]);

      //datacard += TString::Format("CMS_hzg_PU_13TeV     lnN              %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysPU[4][m], 1+sigSysPU[3][m], 1+sigSysPU[2][m], 1+sigSysPU[1][m], 1+sigSysPU[0][m]);

      datacard += TString::Format("CMS_PU_13TeV     lnN              %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysPU[4][m], 1+sigSysPU[3][m], 1+sigSysPU[2][m], 1+sigSysPU[1][m], 1+sigSysPU[0][m]);


      //if(cat!=6789) datacard += TString::Format("CMS_JEC_13TeV     lnN        %-15f          %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_jec[5][m], 1+sys_jec[4][m], 1+sys_jec[3][m], 1+sys_jec[2][m], 1+sys_jec[1][m], 1+sys_jec[0][m]);

      //if(cat!=6789) datacard += TString::Format("CMS_JEC_13TeV     lnN      %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_jec[4][m], 1+sys_jec[3][m], 1+sys_jec[2][m], 1+sys_jec[1][m], 1+sys_jec[0][m]);
      
      ///dont include it for lepton tagged
      if(cat!=6789) {
	
	datacard += TString::Format("CMS_JEC_13TeV     lnN      %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_jec[4][m], 1+sys_jec[3][m], 1+sys_jec[2][m], 1+sys_jec[1][m], 1+sys_jec[0][m]);

	datacard += TString::Format("CMS_JER_13TeV     lnN      %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_jer[4][m], 1+sys_jer[3][m], 1+sys_jer[2][m], 1+sys_jer[1][m], 1+sys_jer[0][m]);

	///UE
	datacard += TString::Format("CMS_hzg_UE_13TeV     lnN      %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_ue[4][m], 1+sys_ue[3][m], 1+sys_ue[2][m], 1+sys_ue[1][m], 1+sys_ue[0][m]);

	///PS
	datacard += TString::Format("CMS_hzg_PS_13TeV     lnN      %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_ps[4][m], 1+sys_ps[3][m], 1+sys_ps[2][m], 1+sys_ps[1][m], 1+sys_ps[0][m]);


      }

      
      ///R9
      if(cat>=1 && cat<=4){
	datacard += TString::Format("CMS_R9_13TeV     lnN      %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_r9[2][m], 1+sys_r9[2][m], 1+sys_r9[2][m], 1+sys_r9[1][m], 1+sys_r9[0][m]);
      }
      
      int imass = m;
      if(mass>=120&&mass<125) imass = 0;
      if(mass>=125 && mass<130) imass = 5;


      
      if(schannel == "mmg"){

	for(int isys=0; isys<=5; isys++){
	  
	  /*if(sysrochor_mean_arr[isys][imass]==0) sysrochor_mean_arr[isys][imass] = 0.001;
	  if(sysem_mean_arr_m[isys][imass]==0) sysem_mean_arr_m[isys][imass] = 0.001;
	  if(syspho_mean_arr_m[isys][imass]==0) syspho_mean_arr_m[isys][imass] = 0.001;


	  if(sysrochor_reso_arr[isys][imass]==0) sysrochor_reso_arr[isys][imass] = 0.001;
	  if(sysem_reso_arr_m[isys][imass]==0) sysem_reso_arr_m[isys][imass] = 0.001;
	  if(syspho_reso_arr_m[isys][imass]==0) syspho_reso_arr_m[isys][imass] = 0.001;
	  */

	  
	  if(sysrochor_mean_arr[isys][imass]==0) sysrochor_mean_arr[isys][imass] = num_zero_shapesys;
	  if(sysem_mean_arr_m[isys][imass]==0) sysem_mean_arr_m[isys][imass] =  num_zero_shapesys;
	  if(syspho_mean_arr_m[isys][imass]==0) syspho_mean_arr_m[isys][imass] = num_zero_shapesys;


	  if(sysrochor_reso_arr[isys][imass]==0) sysrochor_reso_arr[isys][imass] = num_zero_shapesys;
	  if(sysem_reso_arr_m[isys][imass]==0) sysem_reso_arr_m[isys][imass] =  num_zero_shapesys;
	  if(syspho_reso_arr_m[isys][imass]==0) syspho_reso_arr_m[isys][imass] = num_zero_shapesys;
	  

	}

	datacard += TString::Format("CMS_llg_delta_muonRochor_mean_chan0_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_mean_arr[0][imass]);
	datacard += TString::Format("CMS_llg_delta_muonRochor_mean_chan1_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_mean_arr[1][imass] );
	datacard += TString::Format("CMS_llg_delta_muonRochor_mean_chan2_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_mean_arr[2][imass]);
	datacard += TString::Format("CMS_llg_delta_muonRochor_mean_chan3_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_mean_arr[3][imass] );
	datacard += TString::Format("CMS_llg_delta_muonRochor_mean_chan4_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_mean_arr[4][imass]);
	//datacard += TString::Format("CMS_hzg_delta_muonRochor_mean_chan5_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_mean_arr[5][imass] );

	
	datacard += TString::Format("CMS_llg_delta_muonRochor_sigma_chan0_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_reso_arr[0][imass]);
	datacard += TString::Format("CMS_llg_delta_muonRochor_sigma_chan1_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_reso_arr[1][imass]);
	datacard += TString::Format("CMS_llg_delta_muonRochor_sigma_chan2_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_reso_arr[2][imass]);
	datacard += TString::Format("CMS_llg_delta_muonRochor_sigma_chan3_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_reso_arr[3][imass]);
	datacard += TString::Format("CMS_llg_delta_muonRochor_sigma_chan4_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_reso_arr[4][imass]);
	//datacard += TString::Format("CMS_hzg_delta_muonRochor_reso_chan5_m%d_cat%d    param 0 %-15f\n",  m, cat, sysrochor_reso_arr[5][imass]);



	}


      if(schannel == "eeg" || schannel=="eeg_mmg"){

	for(int isys=0; isys<=5; isys++){
	  
	  /*if(sysem_mean_arr_e[isys][imass]==0) sysem_mean_arr_e[isys][imass] = 0.001;
	  if(syspho_mean_arr_e[isys][imass]==0) syspho_mean_arr_e[isys][imass] = 0.001;


	  if(sysem_reso_arr_e[isys][imass]==0) sysem_reso_arr_e[isys][imass] = 0.001;
	  if(syspho_reso_arr_e[isys][imass]==0) syspho_reso_arr_e[isys][imass] = 0.001;
	  */


	  if(sysem_mean_arr_e[isys][imass]==0) sysem_mean_arr_e[isys][imass] =  num_zero_shapesys;
	  if(syspho_mean_arr_e[isys][imass]==0) syspho_mean_arr_e[isys][imass] = num_zero_shapesys;
	  
	  
	  if(sysem_reso_arr_e[isys][imass]==0) sysem_reso_arr_e[isys][imass] =  num_zero_shapesys;
	  if(syspho_reso_arr_e[isys][imass]==0) syspho_reso_arr_e[isys][imass] = num_zero_shapesys;
	  
	}

	  datacard += TString::Format("CMS_llg_delta_eleEM_mean_chan0_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_mean_arr_e[0][imass]);
	  datacard += TString::Format("CMS_llg_delta_eleEM_mean_chan1_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_mean_arr_e[1][imass] );
	  datacard += TString::Format("CMS_llg_delta_eleEM_mean_chan2_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_mean_arr_e[2][imass]);
	  datacard += TString::Format("CMS_llg_delta_eleEM_mean_chan3_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_mean_arr_e[3][imass] );
	  datacard += TString::Format("CMS_llg_delta_eleEM_mean_chan4_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_mean_arr_e[4][imass]);
	  //datacard += TString::Format("CMS_hzg_delta_eleEM_mean_chan5_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[5][imass] );



	  datacard += TString::Format("CMS_llg_delta_eleEM_sigma_chan0_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_reso_arr_e[0][imass]);
	  datacard += TString::Format("CMS_llg_delta_eleEM_sigma_chan1_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_reso_arr_e[1][imass]);
	  datacard += TString::Format("CMS_llg_delta_eleEM_sigma_chan2_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_reso_arr_e[2][imass]);
	  datacard += TString::Format("CMS_llg_delta_eleEM_sigma_chan3_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_reso_arr_e[3][imass]);
	  datacard += TString::Format("CMS_llg_delta_eleEM_sigma_chan4_m%d_cat%d    param 0 %-15f\n",  m, cat, sysem_reso_arr_e[4][imass]);
	  //datacard += TString::Format("CMS_hzg_delta_eleEM_reso_chan5_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[5][imass]);

      }
   

      /////add bkg uncertainties
      //CMS_hzg_mass

      //TIterator* iter = list->createIterator();

      if(floatbkg){
	RooAbsPdf *bgrfit_ws = wspace.pdf(Form("pdf_bgr_cat%d",cat)); 
	RooArgSet *list = bgrfit_ws->getParameters(*fX);
	cout<<"PRINTING LIST"<<endl;
	list->Print();
	
	TIterator* iter = list->createIterator();
	
	for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
	  RooRealVar *v = dynamic_cast<RooRealVar *>(a);
	  
	  string name =  v->GetName();
	  int fstr = name.find("CMS_hzg_mass");
	  if(fstr!=string::npos) continue;
	  
	  datacard += TString::Format("%s flatParam\n", v->GetName());
	  
	}
      }//if(floatbkg)

      datacard += TString::Format("pdf_bgr_cat%d_norm flatParam\n", cat);
      ///set the bkg parameter constant
      
      
      


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
      FILE* out = fopen(TString::Format("output/datacards/datacard_hzg_%s_%s_13TeV_%d.txt",
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
	 
	 int p = 0;
	 if(cat==5) p = 1;
	 
	 double sigeff = sigeff_All[p][m * 5 + k];
	 if(cat==0) sigeff = sigeff_ggf_vbf[m * 5 + k];
	 //cout<<sigeff_ggf_vbf[m * 5 + k]<<",";
	 cout<<sigeff<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;

     cout<<"====related systematics(ID,ISO,PU,HLT)===="<<endl;
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;

	 cout<<sigsys_ggf_vbf[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;


     cout<<"====related systematics(ID,ISO,PU,HLT, Theory - UE/PS)===="<<endl;
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;

	 cout<<sigsys_ggf_vbf_theoryalso[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;


     cout<<"====related systematics(Theory - UE/PS)===="<<endl;
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;

	 cout<<sigsys_ggf_vbf_theoryonly[m * 5 + k]<<",";
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
   /*
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
     */
     

   cout<<"Fetching result now and plotting......"<<endl;
   //RooPlot *plot = fX->frame();


   
   

   
   //fX->setRange("signal",124.0,126.0);
   fX->setRange("unblindR1",115,blind_min);
   fX->setRange("unblindR2",blind_max,170.0);
   
   //RooPlot *plot = fX->frame(RooFit::Title("M_{Z#gamma} distribution"));
   RooPlot *plot = fX->frame(RooFit::Title("   "));

   //dataObs->addColumn(fX);

   cout<<"Printing the original data"<<endl;
   dataObs->Print();

   //dataObs->plotOn(plot,RooFit::CutRange("unblindR1"),RooFit::CutRange("unblindR2"),RooFit::Binning(nBinsForMass));
   cout<<"printing reduced data"<<endl;
   dataObs->Print();
   
   if(blind) dataObs->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
   else dataObs->plotOn(plot,Binning(nBinsForMass));
   


   cout<<"PLotting the data function now"<<endl;
   //dataObs->plotOn(plot);

   cout<<"PLotting the bgrfit function now"<<endl;
   cout<<"address of bgrfit "<<bgrfit<<endl;


   //   dataObs->plotOn(plot);
   bgrfit->plotOn(plot);

   RooHist* hpull = plot->pullHist() ;
   
   if(blind){
     int n = hpull->GetN();
     for(int in=0; in<n; in++){
       
       double x, y;
       hpull->GetPoint(in,x,y);
       cout<<"For point i, x and y : "<<in<<" "<<x <<" "<<y<<endl;
       
       if(x>=blind_min && x<=blind_max) 
	 hpull->SetPoint(in,x,-100);
     }//for(int in=0; in<n; in++)
   }//if(blind)
      

   // Create a new frame to draw the residual distribution and add the distribution to the frame


   //bgrfit_ext->plotOn(plot,RooFit::LineColor(2));

   
   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   rfit->Print();

   cout<<"==================PRINTED r BEFORE PLOTTING===================== "<<endl;

   /*
   TH1F *confInt1s = new TH1F("confInt1s","",nBinsForMass,xmin,xmax);
   confInt1s->SetStats(false);
   confInt1s->SetFillColor(kGreen-4);
   
   TH1F *confInt2s = new TH1F("confInt2s","",nBinsForMass,xmin,xmax);
   confInt2s->SetStats(false);
   confInt2s->SetFillColor(kYellow-4);

   rfit->GetFitter()->GetConfidenceIntervals(confInt1s, 0.6827);
   rfit->GetFitter()->GetConfidenceIntervals(confInt2s, 0.9545);
   */

   //bgrfit->plotOn(plot);
   ///error bands
   ///1sigma

   if(plotRooFITErrBands){
     bgrfit->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*rfit,2),RooFit::FillColor(kYellow-4));
     bgrfit->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*rfit,1),RooFit::FillColor(kGreen-4));
   }

   bgrfit->plotOn(plot,RooFit::Name("central")); 
   
   if(blind) dataObs->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
   else dataObs->plotOn(plot,Binning(nBinsForMass));
   
   TGraphAsymmErrors *oneSigmaBand = new TGraphAsymmErrors();
   TGraphAsymmErrors *twoSigmaBand = new TGraphAsymmErrors();
   
   RooCurve *nomBkgCurve = plot->getCurve("central");
   int p=0;
   
   double mhLow = xmin;
   double mhHigh = xmax;
   double massStep = 0.5;
   double nllTolerance = 0.05;
   
   //cout<<"Error band"
   if(!plotRooFITErrBands){
     for (double mass=double(mhLow); mass<double(mhHigh)+massStep; mass+=massStep) {
       double lowedge = mass-0.5;
       double upedge = mass+0.5;
       double center = mass;
       double nomBkg = nomBkgCurve->interpolate(center);
       
       // sensible range
       double lowRange = TMath::Max(0.,nomBkg - 3*TMath::Sqrt(nomBkg));
       double highRange = nomBkg + 3*TMath::Sqrt(nomBkg);
       
       double nllBest = getNLL(fX,bgrfit,dataObs,nomBkg,lowedge,upedge);
       
       double nllBest_full = getNLL(fX,bgrfit,dataObs,nomBkg,xmin, xmax);
       
       cout<<"mass nllBest nllBest_full "<<mass<<" "<<nllBest<<" "<<nllBest_full<<endl;
       
       double errLow1Value,errHigh1Value,errLow2Value,errHigh2Value;
       // cant handle having 0 events
       if (nomBkg<1.e-5) {
	 errLow1Value = 0.;
	 errLow2Value = 0.;
	 errHigh1Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,1.,nllTolerance);
	 errHigh2Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,4.,nllTolerance);
       }
       
       else {
	 // error calc algo
	 if (verbose_) cout<< "[INFO] " << "errLow1" << endl;
	 errLow1Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,lowRange,lowedge,upedge,1.,nllTolerance);
	 if (verbose_) cout<< "[INFO] " << "errLow2" << endl;
	 errLow2Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,lowRange,lowedge,upedge,4.,nllTolerance);
	 
	 if (verbose_) cout<< "[INFO] " << "errHigh1" << endl;
	 errHigh1Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,1.,nllTolerance);
	 
	 if (verbose_) cout<< "[INFO] " << "errHigh2" << endl;
	 errHigh2Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,4.,nllTolerance);
	 
	 
     }
       
       double errLow1 = nomBkg - errLow1Value;
       double errHigh1 = errHigh1Value - nomBkg;
       double errLow2 = nomBkg - errLow2Value;
       double errHigh2 = errHigh2Value - nomBkg;
       
       //cout<<"mass errLow1 errHigh1 errLow2 errHigh2nom  "<<mass<<" "<<errLow1<<" "<<errHigh1<<" "<<errLow2<<" "<<errHigh2<<" "<<nomBkg<<endl;
       oneSigmaBand->SetPoint(p,center,nomBkg);
       //oneSigmaBand_r->SetPoint(p,center,0);
       twoSigmaBand->SetPoint(p,center,nomBkg);
       //twoSigmaBand_r->SetPoint(p,center,0);
     oneSigmaBand->SetPointError(p,0.,0.,errLow1,errHigh1);
     //oneSigmaBand_r->SetPointError(p,0.,0.,errLow1,errHigh1);
     twoSigmaBand->SetPointError(p,0.,0.,errLow2,errHigh2);
     //twoSigmaBand_r->SetPointError(p,0.,0.,errLow2,errHigh2);
     p++;
     }
   }//if(!plotRooFITErrBands)

   //dataObs->plotOn(plot,RooFit::CutRange("unblindR1"),RooFit::CutRange("unblindR2"),RooFit::Binning(nBinsForMass));
   /*
   cout<<"=====checking the error bands now ======"<<endl;
   cout<<"channel : cat : "<<schannel<<" "<<cats<<endl;

   RooCurve *central = plot->getCurve("central");
   RooCurve *bounds = plot->getCurve("1sigma");

   TGraph* upBound = new TGraph(central->GetN());
   TGraph* loBound = new TGraph(central->GetN());
   TGraph* centralval = new TGraph(central->GetN());
   
   for( int j = 0; j < bounds->GetN(); ++j ){

     if( j < central->GetN() )
       
       upBound->SetPoint(j, bounds->GetX()[j], bounds->GetY()[j]);

     else

       loBound->SetPoint(j, bounds->GetX()[j], bounds->GetY()[j]);

   }

   for( int j = 0; j < central->GetN(); ++j ){

     centralval->SetPoint(j, central->GetX()[j], central->GetY()[j]);
   }

   cout<<"Central val at 120 : 136 : 137 : 138 : "<<centralval->Eval(120)<<" "<<centralval->Eval(136)<<" "<<centralval->Eval(137)<<" "<<centralval->Eval(138)<<endl;
   cout<<"Err at 120 : 136 : 137 : 138 : "<<upBound->Eval(120)<<" "<<upBound->Eval(136)<<" "<<upBound->Eval(137)<<" "<<upBound->Eval(138)<<endl;
   */

   TCanvas *c = setTCanvasNicev1("ccat");
   //c->Divide(1,2);
   
   

   //TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.4,1.0,1.0,21);
   //TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.398,22);

   //TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
   //TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.198,22);

   TPad *pad1;
   TPad *pad2;

   if(drawPulldistribution){
   pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0,21);
   pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.350,22);
   
   pad1->SetFillColor(0);
   pad2->SetFillColor(0);
   
   
   pad1->Draw();
   pad2->Draw();
   
   
   //c->cd(1); 
   pad1->cd();
   }

   if(!drawPulldistribution)
     c->cd();

   gPad->SetLeftMargin(0.15);
   plot->GetYaxis()->SetTitleOffset(1.6);
   
   double max = plot->GetMaximum();
   plot->SetMaximum(max*1.5);
   //plot->SetMinimum(0);
   if(blind) plot->SetMinimum(0.0001);
   if(!blind) plot->SetMinimum(0);
   
   //plot->GetXaxis()->SetTitle("m_{Z#gamma} [GeV]");

   if(schannel == "eeg") plot->GetXaxis()->SetTitle("m_{ee#gamma} [GeV]");
   if(schannel == "mmg") plot->GetXaxis()->SetTitle("m_{#mu#mu#gamma} [GeV]");
   if(schannel == "eeg_mmg") plot->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");

   plot->GetXaxis()->SetTitleSize(0.05);
   //plot->GetXaxis()->SetLabelFont(42);
   //plot->GetXaxis()->SetTitleFont(42);

   plot->GetXaxis()->SetLabelFont(22);
   plot->GetXaxis()->SetTitleFont(22);

   plot->GetYaxis()->SetTitle("Events / GeV");
   plot->GetYaxis()->SetTitleSize(0.05);

   plot->GetYaxis()->SetTitleOffset(1);

   //plot->GetYaxis()->SetLabelFont(42);
   //plot->GetYaxis()->SetTitleFont(42);

   plot->GetYaxis()->SetLabelFont(22);
   plot->GetYaxis()->SetTitleFont(22);


   plot->GetXaxis()->SetLabelSize(0.05);
   plot->GetYaxis()->SetLabelSize(0.05);

   plot->GetXaxis()->SetTitleSize(0.05);
   plot->GetYaxis()->SetTitleSize(0.05);


   //int status = rbkg->status();
   TLatex *lat = new TLatex();
   lat->SetNDC();
   lat->SetTextFont(42);
   //lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));


   TMathText *latm = new TMathText();
   latm->SetNDC();
   latm->SetTextFont(42);

   string text = "";
   if (schannel == "eeg") text = "ee#gamma";
   if (schannel == "mmg") text = "#mu#mu#gamma";
   //if (schannel == "eeg_mmg") text = "ee#gamma + #mu#mu#gamma";
   if (schannel == "eeg_mmg") text = "ll#gamma";

   
   string tmpcat = Form("%d",cat);
   if(cat==1) 
     tmpcat = "Untagged 1";

   if(cat==2) 
     tmpcat = "Untagged 2";

   if(cat==3) 
     tmpcat = "Untagged 3";

   if(cat==4) 
     tmpcat = "Untagged 4";

   if(cat==5)
     tmpcat = "Dijet tag";
   if(cat==6789)
     tmpcat = "Lepton tag";
   if(cat==10)
     tmpcat = "Boosted";

   //lat->DrawLatex(0.1,0.92,Form("Category  %s",tmpcat.c_str()));

   plot->Draw();
   //lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Fit Status = %d ",chi2,status));
   
   if(!plotRooFITErrBands){ 
     
     /*oneSigmaBand->SetFillColor(kGreen);
     oneSigmaBand->SetLineColor(kGreen);
          
     twoSigmaBand->SetFillColor(kYellow);
     twoSigmaBand->SetLineColor(kYellow);
     */

     ///CWR
     oneSigmaBand->SetFillColor(kGreen+1);
     oneSigmaBand->SetLineColor(kGreen+1);
          
     twoSigmaBand->SetFillColor(kOrange);
     twoSigmaBand->SetLineColor(kOrange);
     ///CWR

     twoSigmaBand->Draw("L3 same");
     oneSigmaBand->Draw("L3 same");
   }

   c->Modified();
   
   c->Update();

   plot->Draw("same");
   ///during PUBCOMM

   //c->GetFrame()->Draw();
   c->Modified();
   c->Update();

   lat->DrawLatex(0.35,0.85,Form("#font[22]{#scale[0.85]{H#rightarrow Z#gamma#rightarrow %s }}",text.c_str()));

   //latm->DrawMathText(0.35,0.85,Form("#font[22]{#scale[0.85]{H\\rightarrow Z\\gamma\\rightarrow %s }}",text.c_str()));
   c->Modified();
   c->Update();
   lat->DrawLatex(0.37,0.80,Form("#font[22]{#scale[0.8]{%s}}",tmpcat.c_str()));
   //latm->DrawMathText(0.37,0.80,Form("#font[22]{#scale[0.8]{%s}}",tmpcat.c_str()));
   c->Modified();
   c->Update();
   
   int iproc = 0;
   if(cat==5)
     iproc = 1;
   if(cat==6789)
     iproc = 4;
   
   double integ = hMasssig[iproc][5]->Integral();
   hMasssig[iproc][5]->Scale(expected[iproc][5]*10/integ);
   cout<<"Signal=====signal integ is "<<hMasssig[iproc][5]->Integral()<<endl;
   hMasssig[iproc][5]->SetLineColor(2);
   hMasssig[iproc][5]->SetLineWidth(2);

   hMasssig[iproc][5]->Draw("sameHIST");   

   //confInt1s->Draw("same");
   //confInt2s->Draw("same");

   int iPeriod = 4;
   int iPos = 11;
   //int iPos = 10;
   //int iPos = 40;
   CMS_lumi( c, iPeriod, iPos );
   c->Update();
   c->Modified();
   c->Update();
   
   c->RedrawAxis();
   c->Update();

   TH1F *h1 = new TH1F("h1","",1,1,2);
   h1->SetFillColor(kGreen-4);
   TH1F *h2 = new TH1F("h2","",1,1,2);
   h2->SetFillColor(kYellow-4);

   
   //TLegend *leg = new TLegend(0.6080402,0.7125436,0.8994975,0.8954704,NULL,"brNDC");
   //TLegend *leg = new TLegend(0.6080402,0.6025436,0.8994975,0.8954704,NULL,"brNDC");
   TLegend *leg = new TLegend(0.6300402,0.6025436,0.9394975,0.8954704,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

   TLegendEntry* entry = leg->AddEntry(dataObs, "Data","lp");
   entry->SetMarkerStyle(20);
   
   entry = leg->AddEntry(bgrfit, "Background model","l");
   entry->SetLineColor(4);
   entry->SetLineWidth(2);

   //leg->AddEntry(h1,"#pm1#sigma");
   //leg->AddEntry(h2,"#pm2#sigma");   

   //leg->AddEntry(h1,"#pm1 st. dev.");
   //leg->AddEntry(h2,"#pm2 st. dev.");   

   leg->AddEntry(oneSigmaBand,"#pm1 st. dev.","f");
   leg->AddEntry(twoSigmaBand,"#pm2 st. dev.","f");   


   leg->AddEntry(hMasssig[iproc][5], "Expected signal #times 10","l");


   /*entry = leg->AddEntry(hMasssig[iproc][5],Form("Category:  %s",tmpcat.c_str()),"" );
   entry->SetLineColor(1);
   entry->SetLineStyle(2);
   */
   //entry->SetMarkerColor(0);

   leg->Draw();
   c->Update();
   c->Modified();
   c->Update();

   if(drawPulldistribution){
     RooPlot* frame2 = fX->frame(Title("Pull distribution")) ;
     frame2->addPlotable(hpull,"P") ;
     
     
     //c->cd(2);
     pad2->cd();
     gPad->SetLeftMargin(0.15);
     frame2->GetYaxis()->SetTitleOffset(1.6);
     frame2->SetMinimum(-5);
     frame2->Draw();
     TLine *l = new TLine(xmin,0,170,0);
     l->SetLineWidth(2);
     l->SetLineColor(2);
     l->Draw("same");
   }



   
   
   char *outfilename = new char[100];
   char dirName[100] = "plots";
   

   string scat(cats);
   sprintf(outfilename,"%s/%s.gif",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   sprintf(outfilename,"%s/%s.pdf",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);

   sprintf(outfilename,"%s/%s.root",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   /*
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
   */
   
}

////end of simple mkdatacards for chk

/////////////////////////////call here specifically for lepton tagged//////
void statAn::simplemkdatacard_LT(const char* channel, int cat)
{

  string sigmodel = "RooCBxGaus";
  bool _nuisance = false;

  double lumi =mylumi; //pb-1
  // cat1, cat2, cat3 or cat4
  TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);


   ///read files for systematics
   ///ROCHOR
   ifstream fsys_mu; //sysMuon_hzg_cat1_8TeV.txt

   double sysrochor_mean_arr[10][20];
   double sysrochor_reso_arr[10][20];

     
   if(cat>=1&&cat<=5 && schannel == "mmg"){
     fsys_mu.open(TString::Format("sysMuon_hzg_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_mu.eof()){
       
       string proc;
       int mass;
       double sysrochor_mean, sysrochor_reso;
       
       fsys_mu >> proc >> mass >> sysrochor_mean >> sysrochor_reso;
       
       
       cout<<"reading sysMuon file "<< TString::Format("sysMuon_hzg_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean sysrochor_reso "<<proc<<" "<<mass<<" "<<sysrochor_mean<<" "<<sysrochor_reso<<endl;
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sysrochor_mean_arr[iproc][5*imass] = sysrochor_mean;
       sysrochor_reso_arr[iproc][5*imass] = sysrochor_reso;
       
     }//while(!datafile.eof())
   }//if(cat>=1&&cat<=5 && schannel == "mmg")
   
     
   
   //////////////////////electron, photons energy scale
   
   ///1. electron channel
   ifstream fsys_emscale_e;
   ifstream fsys_emreso_e; 
   
   double sysem_mean_arr_e[10][20];
   double sysem_reso_arr_e[10][20];
   
   
   ///scale
   if(cat>=1&&cat<=5 && schannel == "eeg"){
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emscale_e.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_emscale_e.eof()){
       
       string proc;
       int mass;
       double sysem_mean_e;
       
       fsys_emscale_e >> proc >> mass >> sysem_mean_e;
	 
       
       cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_e<<endl;
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sysem_mean_arr_e[iproc][5*imass] = sysem_mean_e;
       
     }//while(!datafile.eof())
     
     
       ////reso
       //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emreso_e.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_emreso_e.eof()){
       
       string proc;
       int mass;
       double sysem_reso_e;
       
       fsys_emreso_e >> proc >> mass >> sysem_reso_e;
       
     
     cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_e<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_reso_arr_e[iproc][5*imass] = sysem_reso_e;
     
     }//while(!datafile.eof())

   }//if(cat>=1&&cat<=5 && schannel == "eeg")


   ///2. muon channel
   ifstream fsys_emscale_m;
   ifstream fsys_emreso_m; 

   double sysem_mean_arr_m[10][20];
   double sysem_reso_arr_m[10][20];

   ///scale
   if(cat>=1&&cat<=5 && schannel == "mmg"){
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emscale_m.open(TString::Format("syselepho_hzg_mu_scale_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_emscale_m.eof()){
       
     string proc;
     int mass;
     double sysem_mean_m;
     
     fsys_emscale_m >> proc >> mass >> sysem_mean_m;

     
     cout<<"reading em reso file "<< TString::Format("syselepho_hzg_mu_scale_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_mean_arr_m[iproc][5*imass] = sysem_mean_m;
     
     }//while(!datafile.eof())


     ////reso
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emreso_m.open(TString::Format("syselepho_hzg_mu_smear_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_emreso_m.eof()){
       
       string proc;
       int mass;
       double sysem_reso_m;
       
       fsys_emreso_m >> proc >> mass >> sysem_reso_m;
       
     
       cout<<"reading em reso file "<< TString::Format("syselepho_hzg_mu_smear_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_reso_arr_m[iproc][5*imass] = sysem_reso_m;
     
     }//while(!datafile.eof())
   
   }//if(cat>=1&&cat<=5 && schannel == "mmg")


   double sys_jec[10][20];
   if(cat>=1&&cat<=5){
     ifstream fsys_jet; //sys_ele_JEC_8TeV.txt

     string name;
     if(schannel=="eeg")
       name="ele";
     
     if(schannel=="mmg")
       name="mu";
     
     fsys_jet.open(TString::Format("sys_%s_JEC_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_jet.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_jet >> proc >> mass >> icat >> sys;
       
       cout<<"reading JEC file "<<TString::Format("sys_%s_JEC_8TeV.txt",schannel.c_str())<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sys_jec[iproc][5*imass] = sys;
     
     }//while(!datafile.eof())
   }



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
   //double xmin = 100;
   //double xmin = 115;
   //double xmax = 170;

   TH1D hDataObs("hDataObs", "", 60, xmin, xmax);
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
   //poldeg = 3;


   //model = "Pow";
   //poldeg = 1;


   model = "Exp";
   poldeg = 1;

   if(cat>=1 && cat<=4) {
     model = "RooPolynomial";
   }
   
   if(cat==1 || cat==3 || cat==4) poldeg = 4;
   if(cat==2) poldeg = 5;
   
   if(cat==5){
     model = "Pow";
     poldeg = 1;
   }
   
   if(cat>=6){
     model = "RooPolynomial";
     poldeg = 2;
   }


   //RooGaussStepBernstein *bgrfit = new RooGaussStepBernstein();
   fX = new RooRealVar("x", "", xmin, xmax); 
   //bgrfit = mybkgfit("RooGaussStepBernstein", poldeg);
   mybkgfit(xmin, xmax,cat, schannel);
   //////try importing here//////////////////


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
   

   ///set the bkg parameter constant
   RooAbsPdf *bgrfit_ws = wspace.pdf(Form("pdf_bgr_cat%d",cat)); 
   RooArgSet *list = bgrfit_ws->getParameters(*fX);
   cout<<"PRINTING LIST"<<endl;
   list->Print();
   
   TIterator* iter = list->createIterator();
   
   for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
     RooRealVar *v = dynamic_cast<RooRealVar *>(a);
     cout<<"Printing my var "<<endl;
     v->Print();
     wspace.var(v->GetName())->setConstant(true); 
   }
   ///set the bkg parameter constant

   
   cout<<"now saved the bkg pdf"<<endl;
   
   
   rfit->SetName(("pdf_bgr_"+catString+"_fitresult"));
  if ( wspace.import(*rfit) )
     FATAL("RooWorkspace::import() fit result failed");



  if(model=="RooGaussStepBernstein"){
    RooArgList bkgpar_list;
    bkgpar_list.add(*mean);
    bkgpar_list.add(*sigma);
    bkgpar_list.add(*stepval);
    bkgpar_list.add(*coeflist);
    RooArgSet bkgpar_set(bkgpar_list);
    
    
    cout<<"Values of parameters after fixing ..."<<endl;
    cout<<"Mean : sigma : stepval : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<stepval->getVal()<<endl;
  }
    
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
   bkgNorm.setConstant(false);

   bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
   if ( wspace.import(bkgNorm) )
     FATAL("RooWorkspace::import() failed");
   



  /////start with the signal now /////////

   cout<<"DOING SIGNAL FIT NOW"<<endl;
   double    expected[nprocess][nmass];
   //const char* proc[nprocess] = {"ggH", "VBF", "WplusH", "WminusH", "ZH", "ttH"};
   const char* proc[nprocess] = {"ggH", "VBF", "WH", "ZH", "ttH"};
   //const char* proc[nprocess] = {"ggH", "VBF", "ZH", "WminusH", "ttH"};
   int mass_exist[nsig] = {120, 125, 130};
   double sigeff_All[nprocess][nmass];
   double sigNraw[nprocess][nmass];
   double sigeff_ggf_vbf[nmass];
   double sigsys_ggf_vbf[nmass];

   ///systematcs ID
   double sigSysonExp_lep[nprocess][nmass];
   double sigSysonExp_pho[nprocess][nmass];

   double sigSysonEff_lep[nprocess][nmass];
   double sigSysonEff_pho[nprocess][nmass];

   double sigTotSys[nprocess][nmass];
   double sigSysPU[nprocess][nmass];
   
   double sigSysEff[nprocess][nmass];

   for (int m=0; m<nsig; m++) {


     int mass = mass_exist[m];

     rooMassSig[0][m*5] = NULL;
     sigfit[0][m*5] = NULL;
     rsig[0][m*5] = NULL;
     

     RooDataSet *tmprooMassSig;
     for (int p = 0; p < nprocess; p++){
       
       cout<<"========================At teh very beginning, m is taking value now================"<<m<<endl;
       
       if((int)m >= (int)3){
	 cout<<"----inside >=3 loop, Will continue"<<endl;
	 continue;
	 cout<<"No, wait!!!! not continued"<<endl;
       }
       
       
     
       
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
       
       if(p==0) tmprooMassSig = (RooDataSet*) rooMass;
       tmprooMassSig->append(*(RooDataSet*)rooMass);
       
       
       ///set the parameters const at this point for these mass points
       
    
       // save the normalization factor which adapts the MC scale to real
       // data; evaluate expected signal yield
       double norm = normfactor_SMHiggs_8TeV(channel, p, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double norm = 1;
       
     
       ///SJ
       double sigEff = hMass.Integral()/ngen(channel,p,mass);
       
       //double sigExp = sigEff * lumi; /// use this for limit on xsec*BR
       double sigExp = hMass.Integral()*norm;///use this for limit on xsec/SMxsec
       
       cout<<"LEPTON TAGGED ---- Cat is sigExp "<<cat<<" "<<sigExp<<endl;
       ///systematics on ID
       double sys_lep = 0;
       double sys_pho = 0;
       getIDsys(filename,  trigEff, sys_lep, sys_pho, cat, true);

       double sys_pu = 0;
       getPUsys(filename,  trigEff, sys_pu, cat, true);
       //sigSysonExp_lep[p][m*5] = sqrt(sys_lep)*norm/sigExp;
       //sigSysonExp_pho[p][m*5] = sqrt(sys_pho)*norm/sigExp;
       
       sigSysPU[0][m*5] += pow(sys_pu/hMass.Integral(),2);
       
       sigSysonExp_lep[p][m*5] = sqrt(sys_lep)/hMass.Integral();
       sigSysonExp_pho[p][m*5] = sqrt(sys_pho)/hMass.Integral();
       

       sigTotSys[0][m*5] += pow(sigSysonExp_lep[p][m*5],2) + pow(sigSysonExp_pho[p][m*5],2);

       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) );

       //includes lumi
       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) +pow(2.7*hMass.Integral()/100.,2) );
       sigSysEff[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) +pow(2.7*hMass.Integral()/100.,2) );

       ///sys on Eff is not always right. Only when sigExp=sigEff*lumi, then it is right. 
       ///Only the above sigSysonExp_lep numbers are right. These are on Eff   as well. So sigSysonExp_lep can be used on Exp and Eff
       sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen(channel,p,mass)*sigExp);
       sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen(channel,p,mass)*sigExp);
       
       expected[p][m*5] = sigExp;
       sigeff_All[p][m*5] = sigEff;
       
       sigNraw[p][m*5] = hMass.Integral();
       
       cout<<"Mass is, Signal eff "<<mass<<" " <<sigEff<<endl;

       cout<<"===p is, channel, cat, expectation=== "<<p<<" "<<channel<<" "<<cat <<" "<<expected[p][m*5]<<endl;
       ///SJ
       
       
       /*
       /////for systematics on resolution and scale
     
        if(schannel=="mmg"){
	 
	 double mean = sigmean->getVal();
	 double sigma = sigsigma->getVal();
	
 	 //for now taking all the lepton and photon energy/mom related sys for other samples (WH, Zh and ttH) from ggF
	 if(p>1){
	   sysrochor_mean_arr[p][5*m] = sysrochor_mean_arr[0][5*m];
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[0][5*m];
	   sysem_mean_arr_m[p][5*m] = sysem_mean_arr_m[0][5*m];
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[0][5*m];
	 }

       

	 cout<<"Mass process cat "<<5*m<<" "<<p<<" "<<cat<<endl;
	 cout<<"sysrochor_reso_arr sysem_reso_arr_m sigma  "<<sysrochor_reso_arr[p][5*m]<<" "<<sysem_reso_arr_m[p][5*m]<<" "<<sigma<<endl;
      
	 sysrochor_mean_arr[p][5*m] =  1+sqrt(pow(sysrochor_mean_arr[p][5*m],2) + pow(sysem_mean_arr_m[p][5*m],2))/mean;
	 sysrochor_reso_arr[p][5*m] = 1+sqrt(pow(sysrochor_reso_arr[p][5*m],2) + pow(sysem_reso_arr_m[p][5*m],2))/sigma;
	 

	 cout<<"final sysrochor_reso_arr is "<<sysrochor_reso_arr[p][5*m]<<endl;

	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("delta_muonRochor_mean_chan%d_m%d_cat%d[1]",p,5*m,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod
	 wspace.factory(Form("delta_muonRochor_sigma_chan%d_m%d_cat%d[1]",p,5*m,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,delta_muonRochor_mean_chan%d_m%d_cat%d)",p,5*m,cat,p,5*m,cat,p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,delta_muonRochor_sigma_chan%d_m%d_cat%d)",p,5*m,cat,p,5*m,cat,p,5*m,cat));
	 
	 //wspace.factory(Form("sum::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,delta_muonRochor_mean_chan%d_m%d_cat%d)",p,5*m,cat,p,5*m,cat,p,5*m,cat));
	 //wspace.factory(Form("sum::sigma_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,delta_muonRochor_sigma_chan%d_m%d_cat%d)",p,5*m,cat,p,5*m,cat,p,5*m,cat));
	 
	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );
	 
       }


	sys_jec[p][5*m] = sys_jec[0][5*m];
	
	///em scale 
       if(schannel=="eeg"){
	 
	 double mean = sigmean->getVal();
	 double sigma = sigsigma->getVal();

	 //for now taking all the lepton and photon energy/mom related sys for other samples (WH, Zh and ttH) from ggF
	 if(p>1){
	   sysem_mean_arr_e[p][5*m] = sysem_mean_arr_e[0][5*m];
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[0][5*m];
	   

	 }

	 sysem_mean_arr_e[p][5*m] =  1+sysem_mean_arr_e[p][5*m]/mean;
	 sysem_reso_arr_e[p][5*m] = 1+sysem_reso_arr_e[p][5*m]/sigma;
	 


	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("delta_eleEM_mean_chan%d_m%d_cat%d[1]",p,5*m,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod
	 wspace.factory(Form("delta_eleEM_sigma_chan%d_m%d_cat%d[1]",p,5*m,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,delta_eleEM_mean_chan%d_m%d_cat%d)",p,5*m,cat,p,5*m,cat,p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,delta_eleEM_sigma_chan%d_m%d_cat%d)",p,5*m,cat,p,5*m,cat,p,5*m,cat));
	 
	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );
	 
       }

       /////for systematics on resolution and scale
       */

       //delete sfit;

     }//for (int p = 0; p < nprocess; p++)       

     rooMassSig[0][5*m] = tmprooMassSig;
     
     sigTotSys[0][5*m] = sqrt(sigTotSys[0][5*m]);
     sigSysPU[0][m*5]  = sqrt(sigSysPU[0][m*5]);

     cout<<"=================Channel, cat "<<channel<<" "<<cat<<endl;
     cout<<"===========fitting for mass =================="<<mass<<endl;
     mysigfit(sigmodel,xmin, xmax, _nuisance, 0, 5*m, cat, schannel);
     
     setSigParConst(sigmodel, _nuisance);
       
     // unique suffix
     TString sfx = TString::Format("sig_%s_%d_cat%d", proc[0], mass,cat);
     
     ///CHECK - set signal parameters to const
     
     rooMassSig[0][m*5]->SetName(TString::Format("signaldata_%s_%d_cat%d",proc[0], mass, cat));
     wspace.import(*rooMassSig[0][m*5],RooFit::RenameVariable("x", "CMS_hzg_mass"));
     
     
     cout<<"Inside the signal mass loop for filling the data set"<<endl;
     cout<<"=========Printing the dataset information========="<<endl;
     rooMassSig[0][m*5]->Print();
     cout<<"=========Printed the dataset information========="<<endl;
     
     
     // - add signal PDF into the workspace under the name "pdf_sig_...";
     // - add suffix "_sig_..." to all subPDFs;
     // - connect the PDF's X axis variable to the X axis variable of data_obs
     // - connect the nuisance parameters of different signal PDFs together
     
     cout<<"setting the signal pdf name "<<endl;
     
     sigfit[0][5*m]->SetName("pdf");
     //sfit->SetName(Form("pdf_sig_%s_%d_cat%d",proc[p], mass, cat));
     if (wspace.import(*sigfit[0][5*m], RooFit::RenameAllNodes(sfx),
			 RooFit::RenameVariable("x", "CMS_hzg_mass")))
	 FATAL("RooWorkspace::import() failed");
       
       
       cout<<"saved the signal pdf"<<endl;

       rsig[0][m*5]->SetName(Form("fitresult_%s_%d_cat%d",proc[0], mass, cat));
       if (wspace.import(*rsig[0][m*5]))
	 FATAL("RooWorkspace::import() failed");
       
       cout<<"saved the signal fit result"<<endl;
       

   } // process and mass loops


   cout<<"nprocess : nsig : "<<nprocess<<" "<<nsig<<endl;

   
   for (int m = 0; m < 3; m++) {
     
     int mass = mass_exist[m];   
     //double sigEff = (sigNraw[0][m * 5] + sigNraw[1][m * 5])/( ngen(channel,0,mass) + ngen(channel,1,mass) );   
     
     double ntotgen = 0;
     double numtot =  0;
     double systot = 0;
     double sumexp = 0;
     for(int ip=0; ip<nprocess; ip++){
       
       double norm = normfactor_SMHiggs_8TeV(channel, ip, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double ngenTot = ngen_Tot(channel, ip, mass); ///total no of ngen e vents = ee mm tt
       ntotgen += ngen(channel,ip,mass)*norm;
       numtot  += sigNraw[ip][m * 5]*norm;
       //systot += pow(sigNraw[ip][m * 5]*norm,2); ///check

       systot += pow(sigSysEff[ip][m * 5]*norm,2); ///check

       ////put the expected = total of all processes
       sumexp += expected[ip][m*5];
       
     }
     double sigEff = numtot/ntotgen;
     sigeff_ggf_vbf[m*5] = sigEff;

     //double sigTotsys_ggf_vbf = sqrt( pow(sigTotSys[0][m * 5],2) + pow(sigTotSys[1][m * 5],2) ); 
     double sigTotsys_ggf_vbf = sqrt( systot )/ntotgen;
     sigsys_ggf_vbf[m*5] = sigTotsys_ggf_vbf;

     ///put the 0th element = total
     expected[0][m*5] = sumexp;
     
     cout<<"Finally for lepton tagged, for mass "<<m*5 <<"expectation is "<<expected[0][m*5]<<endl;
     
   }




  
   ///interpolation
   // use 1 GeV steps
   //for (int p = 0; p < nprocess; p++){ 
   for (int p = 0; p < 1; p++){ 
     
     //for (int p = nprocess-1; p >=0; p--){ // 2 processes
     for (int mm = 0; mm< nsig-1; mm++){ 
       for (int k = 1; k <= 4; k++) {
	 
	 int mass = 120 + mm * 5 + k;
	 
	   
       
	 cout<<""<<endl;
	 cout<<"cat : p : mm : k : mass : channel : "<<cat<<" "<<p<<" "<<mm<<" "<<k<<" "<<mass<<" "<<schannel<<endl;
	 
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

	 sigsys_ggf_vbf[mm * 5 + k] = sigsys_ggf_vbf[mm*5] * a + b * sigsys_ggf_vbf[5*(mm+1)];

	 //systematics
	 sigSysonExp_lep[p][mm * 5 + k] = sigSysonExp_lep[p][mm*5] * a + b * sigSysonExp_lep[p][5*(mm+1)];
	 sigSysonExp_pho[p][mm * 5 + k] = sigSysonExp_pho[p][mm*5] * a + b * sigSysonExp_pho[p][5*(mm+1)];


	 sigTotSys[p][mm * 5 + k] = sigTotSys[p][mm*5]* a + sigTotSys[p][(mm+1)*5] * b;

	 sigSysPU[p][mm * 5 + k] = sigSysPU[p][mm*5]* a + sigSysPU[p][(mm+1)*5] * b;


	 
	 sigsys_ggf_vbf[mm * 5 + k] = sigsys_ggf_vbf[mm*5] * a + b * sigsys_ggf_vbf[5*(mm+1)]; 

	 sigSysonEff_lep[p][mm * 5 + k] = sigSysonEff_lep[p][mm*5] * a + b * sigSysonEff_lep[p][5*(mm+1)];
	 sigSysonEff_pho[p][mm * 5 + k] = sigSysonEff_pho[p][mm*5] * a + b * sigSysonEff_pho[p][5*(mm+1)];
	 
	 ///interpolate all the sys
	 sys_jec[p][mm * 5 + k] = sys_jec[p][mm*5]* a + sys_jec[p][(mm+1)*5] * b;

	 sysrochor_mean_arr[p][mm * 5 + k] = sysrochor_mean_arr[p][mm*5]* a + sysrochor_mean_arr[p][(mm+1)*5] * b;
	 sysrochor_reso_arr[p][mm * 5 + k] = sysrochor_reso_arr[p][mm*5]* a + sysrochor_reso_arr[p][(mm+1)*5] * b;

	 sysem_mean_arr_e[p][mm * 5 + k] = sysem_mean_arr_e[p][mm*5]* a + sysem_mean_arr_e[p][(mm+1)*5] * b;
	 sysem_reso_arr_e[p][mm * 5 + k] = sysem_reso_arr_e[p][mm*5]* a + sysem_reso_arr_e[p][(mm+1)*5] * b;
	 

	 


	 cout<<"========starting to interpolate now for mass "<<mass<<"================="<<endl;
	 cout<<"====masses taken are "<<m1<<" and "<<m2<<endl;
	 siginterpolate(sigmodel, p, mm*5, k, a,b, _nuisance, cat);
	 cout<<"BACK after interpolating"<<endl;
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
       

       if(schannel=="mmg"){
	 
	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("delta_muonRochor_mean_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));
	 wspace.factory(Form("delta_muonRochor_sigma_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,delta_muonRochor_mean_chan%d_m%d_cat%d)",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,delta_muonRochor_sigma_chan%d_m%d_cat%d)",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat));	 	 

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 
       }

	///em scale 
       if(schannel=="eeg"){
	 
	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("delta_eleEM_mean_chan%d_m%d_cat%d[0]",p,5*mm+k,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod
	 wspace.factory(Form("delta_eleEM_sigma_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,delta_eleEM_mean_chan%d_m%d_cat%d)",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,delta_eleEM_sigma_chan%d_m%d_cat%d)",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat));
	 
	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 
       }


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
      //datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      
      /*
      //if(cat>=1&&cat<=4){
      datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:newpdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);
	//}
	*/
      
      /*//21 jan, 2017 - use ggF shape for other processes for cats 1-4 since these events yield is quite small
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
      
      ////when no energy/mom scale sys is applied
      /*datacard += TString::Format("shapes  VBFH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_VBF_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ggH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ggH_%d_%s\n", channel, cats, mass,cats);

      datacard += TString::Format("shapes  ZH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ZH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WminusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_WminusH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  WplusH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_WplusH_%d_%s\n", channel, cats, mass,cats);
      datacard += TString::Format("shapes  ttH         * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_sig_ttH_%d_%s\n", channel, cats, mass,cats);
      */
      
      //      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_bgr\n", channel, cats);
      datacard += TString::Format("shapes  bgr       * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:pdf_bgr_%s\n", channel, cats,cats);
      datacard += TString::Format("shapes  data_obs  * for_datacards_hzg_%s_%s_8TeV.root hzg_workspace:data_obs_%s\n", channel, cats,cats);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %s\n", bin);
      datacard += TString::Format("bin            %s\n", cats);
      //datacard += TString::Format("observation    %d\n", observation);
      datacard += TString::Format("observation    %d\n", -1);
      cout<<"wrote about teh observation ... "<<endl;
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      //datacard += TString::Format("bin            %-15s %-15s %-15s\n", bin, bin, bin);
      /*datacard += TString::Format("bin            %-15s %-15s %-15s\n", cats, cats, cats);
      datacard +=                 "process        VBFH            ggH             bgr\n";
      datacard +=                 "process        -1              0               1\n";
      datacard += TString::Format("rate           %-15f %-15f %-15f\n",
				  //			  expected[1][m], expected[0][m], expectation);
				  expected[1][m], expected[0][m], 1.0);
      */

      cout<<"writing the process expectation"<<endl;
      cout<<"ggF "<<expected[0][m]<<endl;
      cout<<"VBF "<<expected[1][m]<<endl;
      cout<<"W+H "<<expected[2][m]<<endl;
      cout<<"W-H "<<expected[3][m]<<endl;
      cout<<"ZH "<<expected[4][m]<<endl;
      cout<<"ttH "<<expected[5][m]<<endl;
      
      /*
      datacard += TString::Format("bin            %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n", cats, cats, cats,cats, cats,cats, cats);
      datacard +=                 "process        ttH       ZH      WminusH        WplusH           VBFH            ggH             bgr\n";
      datacard +=                 "process        -5        -4       -3             -2                -1             0               1\n";

      cout<<"all the expectations ... "<<expected[0][m]<<" "<<expected[1][m]<<" "<<expected[2][m]<<" "<<expected[3][m]<<" "<<expected[4][m]<<" "<<expected[5][m]<<endl;
      
      datacard += TString::Format("rate           %-15f %-15f %-15f %-15f %-15f %-15f %-15f\n",
				  expected[5][m], expected[4][m], expected[3][m], expected[2][m], expected[1][m], expected[0][m], 1.0);

      datacard += "------------------------------------------------------------------------------------------------------------\n";
      */


      datacard += TString::Format("bin            %-15s %-15s\n", cats, cats);
      datacard +=                 "process        ggH             bgr\n";
      datacard +=                 "process        0               1\n";

      cout<<"all the expectations ... "<<expected[0][m]<<" "<<expected[1][m]<<" "<<expected[2][m]<<" "<<expected[3][m]<<" "<<expected[4][m]<<" "<<expected[5][m]<<endl;
      
      datacard += TString::Format("rate        %-15f %-15f\n",expected[0][m], 1.0);

      datacard += "------------------------------------------------------------------------------------------------------------\n";

      // theoretical uncertainties
      /*
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

      //datacard += "lumi_8TeV          lnN     1.027          1.027   1.027    1.027   1.027          1.027            -\n";
      datacard += "lumi_8TeV          lnN     1.027          -\n";
   
      
      /*
      datacard += TString::Format("CMS_IDeff_8TeV     lnN        %-15f          %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigTotSys[5][m], 1+sigTotSys[4][m], 1+sigTotSys[3][m], 1+sigTotSys[2][m], 1+sigTotSys[1][m], 1+sigTotSys[0][m]);


      datacard += TString::Format("CMS_PU_8TeV     lnN        %-15f          %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sigSysPU[5][m], 1+sigSysPU[4][m], 1+sigSysPU[3][m], 1+sigSysPU[2][m], 1+sigSysPU[1][m], 1+sigSysPU[0][m]);


      datacard += TString::Format("CMS_JEC_8TeV     lnN        %-15f          %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_jec[5][m], 1+sys_jec[4][m], 1+sys_jec[3][m], 1+sys_jec[2][m], 1+sys_jec[1][m], 1+sys_jec[0][m]);



      */




      datacard += TString::Format("CMS_IDeff_8TeV     lnN      %-15f         -\n", 1+sigTotSys[0][m]);


      datacard += TString::Format("CMS_PU_8TeV     lnN       %-15f         -\n",1+sigSysPU[0][m]);


      //datacard += TString::Format("CMS_JEC_8TeV     lnN        %-15f          %-15f         %-15f    %-15f          %-15f          %-15f         -\n",1+sys_jec[5][m], 1+sys_jec[4][m], 1+sys_jec[3][m], 1+sys_jec[2][m], 1+sys_jec[1][m], 1+sys_jec[0][m]);

      /*
      int imass = m;
      if(mass>=120&&mass<125) imass = 0;
      if(mass>=125 && mass<130) imass = 5;


      
      if(cat>=1&&cat<=5 && schannel == "mmg"){
	  datacard += TString::Format("delta_muonRochor_mean_chan0_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_mean_arr[0][imass]);
	  datacard += TString::Format("delta_muonRochor_mean_chan1_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_mean_arr[1][imass] );
	  datacard += TString::Format("delta_muonRochor_mean_chan2_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_mean_arr[2][imass]);
	  datacard += TString::Format("delta_muonRochor_mean_chan3_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_mean_arr[3][imass] );
	  datacard += TString::Format("delta_muonRochor_mean_chan4_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_mean_arr[4][imass]);
	  datacard += TString::Format("delta_muonRochor_mean_chan5_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_mean_arr[5][imass] );

	  datacard += TString::Format("delta_muonRochor_reso_chan0_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_reso_arr[0][imass]);
	  datacard += TString::Format("delta_muonRochor_reso_chan1_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_reso_arr[1][imass]);
	  datacard += TString::Format("delta_muonRochor_reso_chan2_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_reso_arr[2][imass]);
	  datacard += TString::Format("delta_muonRochor_reso_chan3_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_reso_arr[3][imass]);
	  datacard += TString::Format("delta_muonRochor_reso_chan4_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_reso_arr[4][imass]);
	  datacard += TString::Format("delta_muonRochor_reso_chan5_m%d_cat%d    param 1 %-15f\n",  m, cat, sysrochor_reso_arr[5][imass]);


	}


      if(cat>=1&&cat<=5 && schannel == "eeg"){
	  datacard += TString::Format("delta_eleEM_mean_chan0_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[0][imass]);
	  datacard += TString::Format("delta_eleEM_mean_chan1_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[1][imass] );
	  datacard += TString::Format("delta_eleEM_mean_chan2_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[2][imass]);
	  datacard += TString::Format("delta_eleEM_mean_chan3_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[3][imass] );
	  datacard += TString::Format("delta_eleEM_mean_chan4_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[4][imass]);
	  datacard += TString::Format("delta_eleEM_mean_chan5_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_mean_arr_e[5][imass] );


	  datacard += TString::Format("delta_eleEM_reso_chan0_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[0][imass]);
	  datacard += TString::Format("delta_eleEM_reso_chan1_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[1][imass]);
	  datacard += TString::Format("delta_eleEM_reso_chan2_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[2][imass]);
	  datacard += TString::Format("delta_eleEM_reso_chan3_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[3][imass]);
	  datacard += TString::Format("delta_eleEM_reso_chan4_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[4][imass]);
	  datacard += TString::Format("delta_eleEM_reso_chan5_m%d_cat%d    param 1 %-15f\n",  m, cat, sysem_reso_arr_e[5][imass]);

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
   
   
   
}



//////////////////////////////////////end of calling lepton tagged mkdatacards////////////

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
  //double lumi = 2689.8; ///pb-1
  double lumi =mylumi; //pb-1

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
   //statAn s("RooGaussStepBernstein", _poldeg); ;
   statAn s("RooBernstein", _poldeg); ;

   // actual work
   //for (int cat = 1; cat <=4; cat++) {
   for (int cat = 1; cat <=5; cat++) {
   //for (int cat = 5; cat <=5; cat++) {
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
   s.simplemkdatacard("eeg", 0);
   cout<<"=========Muon============"<<endl;
   s.simplemkdatacard("mmg", 0);
   

   ////for lepton tagged category
   cout<<"=========Electron LEPTON TAGGED============"<<endl;
   s.simplemkdatacard("eeg_mmg", 6789);
   //cout<<"=========Muon LEPTON TAGGED============"<<endl;
   //s.simplemkdatacard("mmg", 6789);


   ///   
   cout<<"=========DOING for CAT 10============"<<endl;
   s.simplemkdatacard("eeg", 10);
   s.simplemkdatacard("mmg", 10);
   
  

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



///get expectation within some range of the Mlg mass
///////
///scat tells if we want untagged or the tagged one 
///scat is either "untagged" OR "all"

///get expectation within some range of the Mlg mass
///////
///scat tells if we want untagged or the tagged one 
///scat is either "untagged" OR "all"

void getExpectedEvents(int cat, string scat, string channel, double xmin, double xmax, double sigma_away, bool useRange) ///if useRange is set to true, then the events are calculated within that range else within the sig_eff
{
  //cat 0 is the combined category: 1,2,3,4,5

  std::cout<<"Inside getExpectedEvents; cat is "<<cat <<std::endl;

  //const int totNcat = 9;
  const int totNcat = 10;
  CalcSigma t;

  double Estbkg = 0;
  double data = 0;
  double sig125_ggf = 0;
  double sig125_vbf = 0;
  double sig125_lep = 0;

 
  double sigN[totNcat] = {0,0,0,0,0,0,0,0,0,0};
  double sigTotN = 0;
  string infilename;

  double dataN[totNcat] = {0,0,0,0,0,0,0,0,0,0};
  double dataTotN = 0;

  if(channel == "eeg") infilename = "files_ele.list";
  if(channel == "mmg") {
    infilename = "files_mu.list";
  }


  string infile = "";

  double sigma_ggF[totNcat];
  double sigma_VBF[totNcat];
  double sigma_lepTag[totNcat];
  double sigma_BT[totNcat];


  //ggF
  if(channel=="eeg")
    infile = Form("minitree_ele_sigggH_%d_out.root",125);

  else if(channel=="mmg")
    infile = Form("minitree_mu_sigggH_%d_out.root",125);
  
  else cout<<"Channel WRONGLY given, give either eeg or mmg"<<endl;

  cout<<"==========ggF======================="<<endl;
  //for(int ii=1; ii<=4; ii++)
  //for(int ii=1; ii<=5; ii++)
  for(int ii=1; ii<=totNcat; ii++) ///6thfeb
  //for(int ii=cat; ii<=cat; ii++)
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

  cout<<"==========VBF======================="<<endl;
  //for(int ii=1; ii<=4; ii++)
  for(int ii=1; ii<=totNcat; ii++)
  //for(int ii=cat; ii<=cat; ii++)
    {
      double xmin_s = 0;
      double xmax_s = 0;

      //t.Loop(infile,xmin_s, xmax_s,ii);
      t.Loop(infile,xmin, xmax, ii, xmin_s, xmax_s);
      sigma_VBF[ii-1] = (xmax_s - xmin_s)/2.;
    }

  cout<<"For cat "<<cat<<" and mass 125 GeV, sigma_eff for ggF and VBf "<<sigma_ggF[cat-1]<<" "<<sigma_VBF[cat-1]<<endl;
  
  cout<<"using the same sigma for ggF and VBF and lepton tagged ... ggF"<<endl;

  /// minitree_mu_sigZH_WH_ttH_125_out.root
  ///Lepton tagged
  if(channel=="eeg")
    infile = Form("minitree_ele_sigZH_WH_ttH_%d_out.root",125);

  else if(channel=="mmg")
    infile = Form("minitree_mu_sigZH_WH_ttH_%d_out.root",125);
  
  else cout<<"Channel WRONGLY given, give either eeg or mmg"<<endl;

  cout<<"==========LEPTON TAGGED======================="<<endl;
  //for(int ii=1; ii<=4; ii++)
  //for(int ii=1; ii<=5; ii++)
  for(int ii=1; ii<=totNcat; ii++) ///6thfeb
  //for(int ii=cat; ii<=cat; ii++)
    {
      double xmin_s = 0;
      double xmax_s = 0;
      t.Loop(infile,xmin, xmax, ii, xmin_s, xmax_s);
      sigma_lepTag[ii-1] = (xmax_s - xmin_s)/2.;
    }

  cout<<"For cat "<<cat<<" and mass 125 GeV, sigma_eff using all the samples combined "<<sigma_lepTag[cat-1]<<" "<<sigma_lepTag[cat-1]<<endl;

  //////BOOSTED CATEGORY/////
  if(channel=="eeg")
    infile = Form("minitree_ele_sigggH_%d_out.root",125);

  else if(channel=="mmg")
    infile = Form("minitree_mu_sigggH_%d_out.root",125);
  
  else cout<<"Channel WRONGLY given, give either eeg or mmg"<<endl;

  cout<<"==========BOOSTED TAGGED======================="<<endl;
  //for(int ii=1; ii<=4; ii++)
  //for(int ii=1; ii<=5; ii++)
  for(int ii=1; ii<=totNcat; ii++) ///6thfeb
  //for(int ii=cat; ii<=cat; ii++)
    {
      double xmin_s = 0;
      double xmax_s = 0;
      t.Loop(infile,xmin, xmax, ii, xmin_s, xmax_s);
      sigma_BT[ii-1] = (xmax_s - xmin_s)/2.;
    }

  cout<<"For cat "<<cat<<" and mass 125 GeV, sigma_eff using all the samples combined "<<sigma_BT[cat-1]<<" "<<sigma_BT[cat-1]<<endl;



  /////END OF BOOSTED CATEGORY////


  //VBF
  if(channel=="eeg")
    infile = Form("minitree_ele_sigVBF_%d_out.root",125);

  else if(channel=="mmg")
    infile = Form("minitree_mu_sigVBF_%d_out.root",125);
  
  else cout<<"Channel WRONGLY given, give either eeg or mmg"<<endl;


  double low[totNcat];
  double high[totNcat];
  
  //for(int ii=1; ii<=4; ii++){
  for(int ii=1; ii<=totNcat; ii++){
  //for(int ii=cat; ii<=cat; ii++)
    low[ii-1] = (125-sigma_away*sigma_ggF[ii-1]);
    high[ii-1] = (125+sigma_away*sigma_ggF[ii-1]);

    if(ii==4){
      low[ii-1] = (125-sigma_away*sigma_VBF[ii-1]);
      high[ii-1] = (125+sigma_away*sigma_VBF[ii-1]);
    }

    if(ii>4 && ii<=9){
      low[ii-1] = (125-sigma_away*sigma_lepTag[ii-1]);
      high[ii-1] = (125+sigma_away*sigma_lepTag[ii-1]);
    }
    
    if(ii==10){
      low[ii-1] = (125-sigma_away*sigma_BT[ii-1]);
      high[ii-1] = (125+sigma_away*sigma_BT[ii-1]);
      
    }
    
    cout<<"for cat "<<ii<<" low : high "<<low[ii-1]<<" "<<high[ii-1]<<endl;
  }//for(int ii=1; ii<=totNcat; ii++)

  
  ////try to get the bkg pdf and integrate within the 2sigma window
  char ccat[50];
  sprintf(ccat,"%d",cat);
  string stringcat(ccat);
  
  TFile *f = TFile::Open( ("output/datacards/for_datacards_hzg_"+channel+"_cat"+stringcat+"_13TeV.root").c_str() );
  
  RooWorkspace *ws = (RooWorkspace *)f->Get("hzg_workspace");
  RooAbsData *datapdf = ws->data(Form("data_obs_cat%i",cat));
  datapdf->Print("");
  
  //RooRealVar *mass = ws->var("CMS_hzg_mass");
  RooRealVar *mass = ws->var("CMS_hzg_mass_13TeV");
  mass->Print();
  mass->setBins(100);
  
  
  RooAbsPdf *extbkg = ws->pdf(Form("pdf_bgr_cat%i",cat));

  mass->setRange(low[cat-1], high[cat-1]);
  //mass->setRange();
  //extbkg->createIntegral(*mass);
  
  ////one needs to plot the bkg fit to get the proper integral
  RooPlot *plot = mass->frame();

  double tmpbinning = 55;
  datapdf->plotOn(plot,Binning(tmpbinning),RooFit::Name("data"));
  extbkg->plotOn(plot,RooFit::Name("central")); 
  RooCurve *central = plot->getCurve("central"); 

  RooHist *datahist = (RooHist*) plot->getHist("data"); 

  
  



  ///for now - VBF category - take the same as the untagged sigma 
  
  //low[4] = low[3];
  //high[4] = high[3];

    

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
   //double xsec_ggf[] = {47.38*1.1e-3*br,44.14*1.533e-3*br,41.23*1.941e-3*br}; ///NNLO + NNLL - old
   double xsec_ggf[] = {5.222E+01*1.1e-3*br,4.858E+01*1.533e-3*br,4.531E+01*1.941e-3*br};  ///NNNLO - 18th dec
   double xsec_vbf[] = {3.935*1.1e-3*br, 3.782*1.533e-3*br,3.637*1.941e-3*br};

   double xsec_wph[] = {1.565E+00*1.1e-3*br, 1.373E+00*1.533e-3*br, 1.209E+00*1.941e-3*br};

   //double xsec_wph[] = {9.560E-01*1.1e-3*br, 8.399E-01*1.533e-3*br, 7.413E-01*1.941e-3*br};
   //double xsec_wmh[] = {6.092E-01*1.1e-3*br, 5.327E-01*1.533e-3*br, 4.676E-01*1.941e-3*br};
   double xsec_zh[] = {9.939E-01*1.1e-3*br, 8.839E-01*1.533e-3*br, 7.899E-01*1.941e-3*br};
   double xsec_tth[] = {5.697E-01*1.1e-3*br, 5.071E-01*1.533e-3*br, 4.539E-01*1.941e-3*br};   


   //double lumi = 2569.10 *(1+0.023);
   //double lumi = 2689.8; //pb-1
   double lumi =mylumi; //pb-1

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

   

     double trigEff = 1;
     if( channel == "eeg" ) trigEff = trigEff_ele;
     if( channel == "mmg" ) trigEff = trigEff_mu;
    
   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
     if (tree->GetEntry(ev) <= 0)
       FATAL("TTree::GetEntry() failed");
     
     if (cat > 0 && category != cat)
       continue;
     
     //if (cat==0 && !(category>=1 && category<=5) && scat=="all")
     //if (cat==0 && !(category>=1 && category<=9) && scat=="all")
     if (cat==0 && !(category>=1 && category<=10) && scat=="all")
       continue;

     if (cat==0 && !(category>=1 && category<=4) && scat=="untagged")
       continue;

     if (cat==0 && !(category>=6 && category<=9) && scat=="leptagged")
       continue;


     if (cat==5 && !(category==5))
       continue;


     if (cat==10 && !(category==10))
       continue;
     
          
     if(useRange){
       if(hzg_mass<xmin || hzg_mass>xmax) continue;
     }
     
     else{
       if(hzg_mass<low[category-1] || hzg_mass>high[category-1]){
	 //cout<<"cat is : low : high : hzg "<<category<<" "<<low[category-1]<<" "<<high[category-1]<<" "<<hzg_mass<<endl;
	 continue;  
       }
     }
     

     if(hzg_mass<low[category-1] || hzg_mass>high[category-1]) {

       cout<<"CHECK THE CODE ... Range not rejected "<<endl;
     }
     
     
     //cout<<"came here "<<endl;
 
     if(itype==100 || itype==101) ///100 for electron and 101 for double muon
     //if(itype>=100) ///100 for electron and 101 for double muon
       {
	 data += 1;
	 dataN[category-1]++;
	 dataTotN++;
       }


     ////added after CWR  - 26th march
     double fac_gstar = 0.969;
     //double fac_gstar = 1;
     
     

     double tmp_xsec_vbf = 1;
     if(itype==2) tmp_xsec_vbf = xsec_vbf[1]/fac_gstar;

     double tmp_xsec_ggf = 1;
     if(itype==12) tmp_xsec_ggf = xsec_ggf[1]/fac_gstar;

     double tmp_xsec_wp = 1;
     if(itype==15 || itype==18) tmp_xsec_wp = xsec_wph[1]/fac_gstar;

     //double tmp_xsec_wm = 1;
     //if(itype==18) tmp_xsec_wm = xsec_wmh[1];

     double tmp_xsec_zh = 1;
     if(itype==21) tmp_xsec_zh = xsec_zh[1]/fac_gstar;

     double tmp_xsec_tth = 1;
     if(itype==24) tmp_xsec_tth = xsec_tth[1]/fac_gstar;



     
     
     if(itype<0) Estbkg += mcwei*lumi*trigEff;
     if(itype==2) {
       //cout<<"xsec VBF "<<tmp_xsec_vbf<<endl;
       //sig125_vbf += mcwei*tmp_xsec_vbf*lumi/ngen(channel.c_str(),1,125); //vbF
       //sig125_vbf += tmp_xsec_vbf*lumi/ngen(channel.c_str(),1,125); //vbF

       sig125_vbf += tmp_xsec_vbf*lumi*mcwei*trigEff/ngen(channel.c_str(),1,125); //vbF

       sigN[category-1] = sigN[category-1]+tmp_xsec_vbf*lumi*mcwei*trigEff/ngen(channel.c_str(),1,125);
       sigTotN += tmp_xsec_vbf*lumi*mcwei*trigEff/ngen(channel.c_str(),1,125);
     }
     
     if(itype==12){
       //cout<<"xsec ggF "<<tmp_xsec_ggf<<endl;
       //sig125_ggf += mcwei*tmp_xsec_ggf*lumi/ngen(channel.c_str(),0,125); //ggF
       sig125_ggf += tmp_xsec_ggf*lumi*mcwei*trigEff/ngen(channel.c_str(),0,125); //ggF

       sigN[category-1] = sigN[category-1]+tmp_xsec_ggf*lumi*mcwei*trigEff/ngen(channel.c_str(),0,125);
       sigTotN += tmp_xsec_ggf*lumi*mcwei*trigEff/ngen(channel.c_str(),0,125);
     }


     ///////////lepton tagged
     //1. W+H
     if(itype==15 || itype==18) {
       
       double tmpngen = ngen(channel.c_str(),2,125);
       sig125_lep += tmp_xsec_wp*lumi*mcwei*trigEff/tmpngen;

       sigN[category-1] = sigN[category-1]+tmp_xsec_wp*lumi*mcwei*trigEff/tmpngen;
       sigTotN += tmp_xsec_wp*lumi*mcwei*trigEff/tmpngen;
     }

     /*
     //2. W-H
     if(itype==18) {
       
       double tmpngen = ngen(channel.c_str(),3,125);
       sig125_lep += tmp_xsec_wm*lumi*mcwei*trigEff/tmpngen;

       sigN[category-1] = sigN[category-1]+tmp_xsec_wm*lumi*mcwei*trigEff/tmpngen;
       sigTotN += tmp_xsec_wm*lumi*mcwei*trigEff/tmpngen;
     }
     */
     
     //3. ZH
     if(itype==21) {
       
       //double tmpngen = ngen(channel.c_str(),4,125);
       double tmpngen = ngen(channel.c_str(),3,125);
       sig125_lep += tmp_xsec_zh*lumi*mcwei*trigEff/tmpngen;

       sigN[category-1] = sigN[category-1]+tmp_xsec_zh*lumi*mcwei*trigEff/tmpngen;
       sigTotN += tmp_xsec_zh*lumi*mcwei*trigEff/tmpngen;
     }


     //3. ttH
     if(itype==24) {
       
       //double tmpngen = ngen(channel.c_str(),5,125);
       double tmpngen = ngen(channel.c_str(),4,125);
       sig125_lep += tmp_xsec_tth*lumi*mcwei*trigEff/tmpngen;

       sigN[category-1] = sigN[category-1]+tmp_xsec_tth*lumi*mcwei*trigEff/tmpngen;
       sigTotN += tmp_xsec_tth*lumi*mcwei*trigEff/tmpngen;
     }


   }

   //cout<<"Data observed \t SM expectation \t Signal expectation (ggF) \t Signal expectaion(VBF) \t"<<data<<" \t "<<Estbkg<<" \t "<<sig125_ggf<<" \t "<<sig125_vbf<<endl;

   cout<<"Data observed \t SM expectation \t Signal expectation (ggF) \t Signal expectaion(VBF) \t Lepton tagged "<<data<<" \t "<<Estbkg<<" \t "<<sig125_ggf*trigEff<<" \t "<<sig125_vbf*trigEff<<" "<<sig125_lep*trigEff<<endl;

   //cout<<"BKG expected from BKG pdf "<<extbkg->getVal(*mass)<<endl;
   cout<<"BKG expected from BKG pdf "<<central->Integral()<<endl;
   cout<<"BKG expected from DATAHIST  "<<datahist->Integral()<<endl;
   
   if ( (cat==0 && scat=="all") || (cat==0 && scat=="untagged") )
     for(int icat=1; icat<=totNcat; icat++){
       cout<<""<<endl;
       double perc = sigN[icat-1]/sigTotN;
       cout<<"SIGNAL : frac for icat is "<<icat <<" "<<perc<<endl;


       perc = dataN[icat-1]/dataTotN;
       cout<<"DATA : frac for icat is "<<icat <<" "<<perc<<endl;
     }


   datahist->Print();

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
  //double nev_eeg[6][3] = {
  double nev_eeg[5][3] = {
    /*
    { 33463,33267, 66968},
    {33315,33075,32695}
    */

    /*
    {33463+66519, 33267+129713, 33116+66968}, //ggF
    {33315, 33075+66792, 32524}, //VBF
    {3289, 10014, 3173}, //W+H
    {3221, 10039, 3272}, //W-H
    {6675, 18918, 6408}, //ZH
    {4462, 8675, 3925} //ttH
    */


    {66518, 131180, 66966}, //ggF
    {33313, 66790, 33059}, //VBF
    //{3289, 10014, 3236}, //W+H
    //{3072, 9902, 3272}, //W-H
    //combining W+H and W-H
    {3289+3072, 10014+9902, 3236+3272}, //W-H
    {6675, 19303, 6408}, //ZH
    {4461, 8882, 4388} //ttH


  };

  //double nev_mmg[6][3] = {
  double nev_mmg[5][3] = {

    /*
    {33246,33492,66435 },
    {33398,32695,32948}
    */

    /*
    {33245+66353, 33492+129925, 32934+66434}, //ggF
    {33399, 32695+66482, 32948}, //VBF
    {3376, 9962, 3201}, //W+H
    {3237, 9912, 3361}, //W-H
    {6554, 18856, 6338}, //ZH
    {4500, 8621, 3980} //ttH
    */

    {66350, 131294, 66433}, //ggF
    {33399, 66478, 33424}, //VBF
    //{3376, 10014, 3236}, //W+H
    //{3097, 9757, 3361}, //W-H
    ///combining w+h and w-h
    {3376+3097, 10014+9757, 3236+3361}, //W-H
    {6554, 19205, 6338}, //ZH
    {4502, 8823, 4456} //ttH



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
  
  else  FATAL("wrong channel given");   

}




double statAn::guessNew(RooRealVar *mhzg, RooAbsPdf *pdf, RooAbsData *data, double bestPoint, double nllBest, double boundary, double massRangeLow, double massRangeHigh, double crossing, double tolerance){
  
  //bool verbose_ = true;
  
  bool isLowSide;
  double guess, guessNll, lowPoint,highPoint;
  if (boundary>bestPoint) {
    isLowSide=false;
    lowPoint = bestPoint;
    highPoint = boundary;
  }
  else {
    isLowSide=true;
    lowPoint = boundary;
    highPoint = bestPoint;
  }
  //double prevDistanceFromTruth = 1.e6;
  double distanceFromTruth = 1.e6;
  int nIts=0;
  while (TMath::Abs(distanceFromTruth/crossing)>tolerance) {
    
    //prevDistanceFromTruth=distanceFromTruth;
    guess = lowPoint+(highPoint-lowPoint)/2.;
    guessNll = getNLL(mhzg,pdf,data,guess,massRangeLow,massRangeHigh)-nllBest;
    distanceFromTruth = crossing - guessNll;
    
    if (verbose_) {
      cout << "[INFO] "<< Form("\t lP: %7.3f hP: %7.3f xg: %7.3f yg: %7.3f",lowPoint,highPoint,guess,guessNll) << endl;;
      cout << "[INFO] \t ----------- " << distanceFromTruth/crossing << " -------------" << endl;
    }
    
    // for low side. if nll val is lower than target move right point left. if nll val is higher than target move left point right
    // vice versa for high side
    if (isLowSide){
      if (guessNll>crossing) lowPoint = guess;
      else highPoint=guess;
    }
    else {
      if (guessNll>crossing) highPoint = guess;
      else lowPoint=guess;
    }
    nIts++;
    // because these are envelope nll curves this algorithm can get stuck in local minima
    // hacked get out is just to start trying again
    if (nIts>20) {
      return guess;
      lowPoint = TMath::Max(0.,lowPoint-20);
      highPoint += 20;
      nIts=0;
      if (verbose_) cout << "[INFO] RESET:" << endl;
      // if it's a ridicolous value it wont converge so return value of bestGuess
      if (TMath::Abs(guessNll)>2e4) return 0.; 
    }
  }
  return guess;
}






double statAn::getNLL(RooRealVar *mhzg, RooAbsPdf *pdf, RooAbsData *data, double normVal, double massRangeLow, double massRangeHigh){


  double bestFitNll=1.e8;
  double bestFitNorm;
  RooAbsReal *nll;
  RooExtendPdf *extPdf;
  RooRealVar *normVar = new RooRealVar(Form("%snorm",pdf->GetName()),"",0.,1.e6);

  if (massRangeLow>-1. && massRangeHigh>-1.){
  

    mhzg->setRange("errRange",massRangeLow,massRangeHigh);
    extPdf = new RooExtendPdf(Form("%sext",pdf->GetName()),"",*pdf,*normVar,"errRange");
    nll = extPdf->createNLL(*data,Extended()); //,Range(massRangeLow,massRangeHigh));//,Range("errRange"));
  }
  else {
    extPdf = new RooExtendPdf(Form("%sext",pdf->GetName()),"",*pdf,*normVar);
    nll = extPdf->createNLL(*data,Extended());
  }
  
  if (normVal>-1.){
    normVar->setConstant(false);
    normVar->setVal(normVal);
    normVar->setConstant(true);
  }

  RooMinimizer minim(*nll);
  minim.setStrategy(0);
  //minim.minimize("Minuit2","simplex");
  minim.migrad();

  double corrNll = nll->getVal();
  if (normVal>-1.) normVar->setConstant(false);
  
  return 2*corrNll;
}
