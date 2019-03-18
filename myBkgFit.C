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


bool verbose_ = false;

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

RooAbsData* statAn::fill_events(TH1* hMass, const char* filename, double trigEff, int cat, bool usewei){
  TFile* fi = TFile::Open(filename);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");
	
  //TDirectory * dir = (TDirectory*)fi->Get("diPhoAna");
  //TTree *tree;
  //dir->GetObject("DiPhotonTree",tree);
	
  TTree* tree = dynamic_cast<TTree*> (fi->Get("diPhoAna/diPhotonTree"));
  if (!tree) FATAL("TFile::Get() failed");
  
  Int_t category;
  float hzg_mass;
  float puwei;	
	
  rooMass = new RooDataSet("rooMass", "", RooArgSet(*fX));
  for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
    *fX = hzg_mass;	  
    if(hzg_mass<xmin || hzg_mass>xmax) continue;
           
    hMass->Fill(hzg_mass);
    rooMass->add(RooArgSet(*fX));
	  
  }
  delete tree;
  delete fi;
  return rooMass;
}

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

void mkdatacards(){
  gSystem->SetIncludePath("-I$ROOFITSYS/include");
  //gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
  //gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");
  TH1::AddDirectory(0);
  int _poldeg = 3;
  statAn s("RooBernstein", _poldeg);
	
	
	