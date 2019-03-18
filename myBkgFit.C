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
   int W = 800;
   int H = 600;
   
   int H_ref = 600;
   int W_ref = 800;
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;

   fX->setRange("unblindR1",115,blind_min);
   fX->setRange("unblindR2",blind_max,170.0);
   
   //RooPlot *plot = fX->frame(RooFit::Title("M_{Z#gamma} distribution"));
   RooPlot *plot = fX->frame(RooFit::Title("   "));

   cout<<"Printing the original data"<<endl;
   dataObs->Print();

   if(blind) dataObs->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
   else dataObs->plotOn(plot,Binning(nBinsForMass));

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
	
   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   rfit->Print();
   cout<<"==================PRINTED r BEFORE PLOTTING===================== "<<endl;
   
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
       
       oneSigmaBand->SetPoint(p,center,nomBkg);
       twoSigmaBand->SetPoint(p,center,nomBkg);
       oneSigmaBand->SetPointError(p,0.,0.,errLow1,errHigh1);
       twoSigmaBand->SetPointError(p,0.,0.,errLow2,errHigh2);
       
       p++;
       }
     }//if(!plotRooFITErrBands)	     

     TCanvas *c = setTCanvasNicev1("ccat");
     TPad *pad1;
     TPad *pad2;

     if(drawPulldistribution){
       pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0,21);
       pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.350,22);
   
       pad1->SetFillColor(0);
       pad2->SetFillColor(0);  
       
       pad1->Draw();
       pad2->Draw();
	     
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
     
     plot->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
     plot->GetXaxis()->SetTitleSize(0.05);
     plot->GetXaxis()->SetLabelFont(22);
     plot->GetXaxis()->SetTitleFont(22);
   
     plot->GetYaxis()->SetTitle("Events / GeV");
     plot->GetYaxis()->SetTitleSize(0.05);

     plot->GetYaxis()->SetTitleOffset(1);
     plot->GetYaxis()->SetLabelFont(22);
     plot->GetYaxis()->SetTitleFont(22);
   
     
   plot->GetXaxis()->SetLabelSize(0.05);
   plot->GetYaxis()->SetLabelSize(0.05);

   plot->GetXaxis()->SetTitleSize(0.05);
   plot->GetYaxis()->SetTitleSize(0.05);

   TLatex *lat = new TLatex();
   lat->SetNDC();
   lat->SetTextFont(42);

   TMathText *latm = new TMathText();
   latm->SetNDC();
   latm->SetTextFont(42);

   string text = "";
   if (1) text = "#gamma#gamma";

   plot->Draw();

   if(!plotRooFITErrBands){  
     oneSigmaBand->SetFillColor(kGreen+1);
     oneSigmaBand->SetLineColor(kGreen+1);
          
     twoSigmaBand->SetFillColor(kOrange);
     twoSigmaBand->SetLineColor(kOrange);	   

      twoSigmaBand->Draw("L3 same");
     oneSigmaBand->Draw("L3 same");
   }

   c->Modified();
   c->Update();

   plot->Draw("same");

   c->Modified();
   c->Update();

   lat->DrawLatex(0.35,0.85,Form("#font[22]{#scale[0.85]{H#rightarrow #gamma#gamma#rightarrow %s }}",text.c_str()));
   
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

   int iPeriod = 4;
   int iPos = 11;

   //CMS_lumi( c, iPeriod, iPos );

   c->Update();
   c->Modified();
   c->Update();
   
   c->RedrawAxis();
   c->Update();

   TH1F *h1 = new TH1F("h1","",1,1,2);
   h1->SetFillColor(kGreen-4);
   TH1F *h2 = new TH1F("h2","",1,1,2);
   h2->SetFillColor(kYellow-4);

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

   leg->AddEntry(oneSigmaBand,"#pm1 st. dev.","f");
   leg->AddEntry(twoSigmaBand,"#pm2 st. dev.","f");

   leg->AddEntry(hMasssig[iproc][5], "Expected signal #times 10","l");

   leg->Draw();
   c->Update();
   c->Modified();
   c->Update();

   if(drawPulldistribution){
     RooPlot* frame2 = fX->frame(Title("Pull distribution")) ;
     frame2->addPlotable(hpull,"P") ;
   
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

}
   

//*******************************************************************************************************************************
void myBkgFit(){
  gSystem->SetIncludePath("-I$ROOFITSYS/include");
  //gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
  //gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");
  TH1::AddDirectory(0);
  int _poldeg = 3;
  statAn s("RooBernstein", _poldeg);

}	
	
	
