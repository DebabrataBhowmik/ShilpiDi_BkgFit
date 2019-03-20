#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAxis.h"

using namespace RooFit ;

int nBinsForMass = 55;

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooAbsData *data, std::string name);

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

double getMax(TGraph *g, int n, double max){
  double tmpmax = TMath::MaxElement(n,g->GetY());
  if(tmpmax > max) max = tmpmax;
  return max;
}

double getMin(TGraph *g, int n, double min){
  double tmpmin = TMath::MinElement(n,g->GetY()); 
  if(tmpmin > min) min = tmpmin;
  return min;
}

void plot(string channel, bool blind){

 gSystem->SetIncludePath("-I$ROOFITSYS/include");
 
 
 
 //TFile *f = TFile::Open( ("output/datacards/for_datacards_hzg_"+channel+"_cat"+scat+"_8TeV.root").c_str() );
 TFile *f = TFile::Open("output_WorkSpace_Exp_final.root");
 
 RooWorkspace *ws = (RooWorkspace *)f->Get("mgg_workspace");
 RooAbsData *data = ws->data("data_obs_cat1");
 data->Print("");
 
 RooRealVar *mass = ws->var("CMS_mgg_mass");
 mass->Print();
 
double blind_min = 120;
double blind_max = 130;

mass->setRange("unblindR1",100,blind_min);
mass->setRange("unblindR2",blind_max,170.0);
      
RooPlot *plot = mass->frame(RooFit::Title("M_{#gamma#gamma} distribution"));
if(blind) data->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
else data->plotOn(plot,Binning(nBinsForMass));

TH1F *h = (TH1F*)data->createHistogram("CMS_mgg_mass");

h->SetMarkerStyle(20);
h->SetMarkerColor(1);
h->SetLineColor(1);


cout<<"nbins "<<h->GetNbinsX()<<endl;
for(int ii=1; ii<h->GetNbinsX(); ii++){
	
	int ibin = h->FindBin(blind_min);
	int jbin = h->FindBin(blind_max);
	
	if(ii>=ibin && ii<=jbin) h->SetBinContent(ii,0);
	if(ii>=ibin && ii<=jbin) h->SetBinError(ii,0);
}

 cout<<"Looped over bins"<<endl;
 
 RooDataHist datahist("data", "data", *mass, h);
 
TLegend *legbkg = new TLegend(0.8190955,0.7281022,0.9535176,0.9160584,NULL,"brNDC"); 
  
RooAbsPdf *extbkg[50]; 
int icolbkg[] = {1,2,3,4,5,6,8,9,41,42,43,40,17,30,46,47,40,32,30,30};

RooFitResult *rbkg;

cout<<"plotting Exp "<<endl;

extbkg[1] = ws->pdf("pdf_Exp_bgr_cat1");
extbkg[1]->plotOn(plot,LineColor(icolbkg[1]), RooFit::Name(Form("pdf_Exp")));
legbkg->AddEntry(plot->findObject("pdf_Exp"),"Exp","l");

f->Close();
 
//plot the Pow
 TFile *f2 = TFile::Open("output_WorkSpace_pow_test2.root");
 RooWorkspace *ws2 = (RooWorkspace *)f->Get("mgg_workspace");
 
extbkg[2] = ws->pdf("pdf_Pow_bgr_cat1");
extbkg[2]->plotOn(plot,LineColor(icolbkg[2]), RooFit::Name(Form("pdf_Pow")));
legbkg->AddEntry(plot->findObject("pdf_Pow"),"Pow","l");
 
f2->Close();
 
//plot the Bernstein
cout<<"Doing Bernstein/Roopolynomial now "<<endl;

TFile *f3 = TFile::Open("output_WorkSpace_RooPolynomial_final.root");
RooWorkspace *ws3 = (RooWorkspace *)f->Get("mgg_workspace");
extbkg[3] = ws->pdf("pdf_RooPolynomial_bgr_cat1");
extbkg[3]->plotOn(plot,LineColor(icolbkg[3]), RooFit::Name(Form("pdf_Bernstein"))); 
legbkg->AddEntry(plot->findObject("pdf_Laurent"),"Bernstein","l"); 
 
TGraphAsymmErrors onesigma;
TGraphAsymmErrors twosigma;

char *filename = new char[100];
char dirName[100] = "plots";

TCanvas *c = setTCanvasNicev1("ccat");

gPad->SetLeftMargin(0.15);
plot->GetYaxis()->SetTitleOffset(1.6);

if(blind) plot->SetMinimum(0.0001);

double max = plot->GetMaximum();
plot->SetMaximum(max*1.5);
if(blind) plot->SetMinimum(0.0001);

plot->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");

plot->Draw();

c->Update();
c->Modified();
c->Update();

legbkg->Draw();

c->SaveAs("Final_mgg_Fit.png");

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooAbsData *data, std::string name){


  TRandom3 *RandomGen = new TRandom3();
  
  double prob;
  int ntoys = 500;

  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);
  //norm.removeRange();

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  // get The Chi2 value from the data
  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass),Name("data"));

  pdf->plotOn(plot_chi2,Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  // The first thing is to check if the number of entries in any bin is < 5 
  // if so, we don't rely on asymptotic approximations
 
  if ((double)data->sumEntries()/nBinsForMass < 5 ){

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    // store pre-fit params 
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();
 
    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
      //  std::cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << std::flush;
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);
      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0),RooFit::SumW2Error(kTRUE)); //FIXME

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t);
      pdf->plotOn(plot_t);//,RooFit::NormRange("fitdata_1,fitdata_2"));

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
      delete plot_t;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    sort(toy_chi2.begin(), toy_chi2.end()); 
    
    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForMass-np),toyhist.GetMaximum(),chi2*(nBinsForMass-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    // back to best fit 
    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
  }
  std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}





 
 
