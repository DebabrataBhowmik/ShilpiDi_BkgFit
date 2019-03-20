////usage: root -l -b plot.C\(\"eeg\"\)
///
///   gSystem->SetIncludePath("-I$ROOFITSYS/include");

///   gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
///   gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");

/*
// C a l c u l a t e   c h i ^ 2 
  // ------------------------------

  // Show the chi^2 of the curve w.r.t. the histogram
  // If multiple curves or datasets live in the frame you can specify
  // the name of the relevant curve and/or dataset in chiSquare()
  cout << "chi^2 = " << frame1->chiSquare() << endl ;





// S h o w   r e s i d u a l   a n d   p u l l   d i s t s
  // -------------------------------------------------------

  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = frame1->residHist() ;

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = frame1->pullHist() ;

  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* frame2 = x.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame3 = x.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P") ;

*/


////https://github.com/cms-analysis/flashggFinalFit/blob/master/Background/test/fTest.cpp

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

#include "CMS_lumi.C"

using namespace RooFit ;


//int nBinsForMass = 100;
//int nBinsForMass = 35;
int nBinsForMass = 55;
//int nBinsForMass = 60;

//double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name);
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


  //int ncats = 5;
  //const int ncats = 6;
  //int cat_arr[ncats] = {1,2,3,4,5,6789};


  const int ncats = 7;
  //int cat_arr[ncats] = {1,2,3,4,5,6789,10};
  int cat_arr[ncats] = {1,2,3,4,5,10,6789};

  //int ncats = 2;
  //int ncats = 1;
  
   // FIXME: this line is required on lxplus and pcncu1X machines
   gSystem->SetIncludePath("-I$ROOFITSYS/include");

   gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
   gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");

  ///channel is eeg or mmg


   //blind = true;
   

   int gcol[5] = {12,46,38,8,9};
   int gstyle[5] = {20,21,22,23,34};

   int inicat = 0;
   //int totncat = 5;
   int totncat = 6;
   if(channel=="eeg_mmg"){
     //inicat = 5;
     //totncat = 6;

     inicat = 6;
     totncat = 7;

   }

   //for(int icat=1; icat<ncats; icat++){
   for(int iicat=inicat; iicat<totncat; iicat++){
     
     int icat = cat_arr[iicat];
     
      char ccat[50];
      sprintf(ccat,"%d",icat);
      string scat(ccat);
      
      TFile *f = TFile::Open( ("output/datacards/for_datacards_hzg_"+channel+"_cat"+scat+"_8TeV.root").c_str() );
      
      RooWorkspace *ws = (RooWorkspace *)f->Get("hzg_workspace");
      RooAbsData *data = ws->data(Form("data_obs_cat%i",icat));
      data->Print("");
      
      RooRealVar *mass = ws->var("CMS_hzg_mass");
      mass->Print();
      
      ///define range
      double blind_min = 120;
      double blind_max = 135;
      //mass->setRange("signal",124.0,126.0);
      mass->setRange("unblindR1",100,blind_min);
      mass->setRange("unblindR2",blind_max,170.0);
      
      RooPlot *plot = mass->frame(RooFit::Title("M_{Z#gamma} distribution"));
      if(blind) data->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
      else data->plotOn(plot,Binning(nBinsForMass));
      //data->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
      //data->plotOn(plot,CutRange("unblindR1"),Binning(nBinsForMass));
      
      TH1F *h = (TH1F*)data->createHistogram("CMS_hzg_mass");

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
      //h->Rebin(2);
      RooDataHist datahist("data", "data", *mass, h);
      //datahist.plotOn(plot);

      cout<<"declaring legend "<<endl;
      
      TLegend *legbkg = new TLegend(0.8190955,0.7281022,0.9535176,0.9160584,NULL,"brNDC");      
      

      RooAbsPdf *extbkg[50];
      //int icolbkg[] = {1,2,3,4,5,6,8,9,41,42,43,44,45,30,46,47,40,32,30,30};
      int icolbkg[] = {1,2,3,4,5,6,8,9,41,42,43,40,17,30,46,47,40,32,30,30};
      //int totNdeg = 5;
      int totNdeg = 4;
      //int totNdeg = 2;
      
      int ibkg = 2;

      RooFitResult *rbkg;
      
      /*
      cout<<"plotting ibkg "<<ibkg<<endl;
      //pdf_2_bgr_cat1
      extbkg[ibkg] = ws->pdf(Form("pdf_bgr_cat%i",icat));
      //pdf_bgr_cat1_fitresult_2
      RooFitResult *rbkg = (RooFitResult *)ws->genobj(Form("pdf_bgr_cat%i_fitresult",icat));
      if(rbkg==NULL) cout<<"Null pointer"<<endl;
      //rbkg->Dump();
      
      
      rbkg->Print();
      rbkg->correlationMatrix().Print();
      
      extbkg[ibkg]->plotOn(plot,LineColor(icolbkg[0]), RooFit::Name(Form("pdf")));
      
      legbkg->AddEntry(plot->findObject(Form("pdf")),Form("pol(%d)",ibkg),"l");
      
      cout<<"done plotting bkg "<<endl;
      */

      
      cout<<"plotting Exp "<<endl;
      //plot the Exp
      extbkg[totNdeg+1] = ws->pdf(Form("pdf_Exp_bgr_cat%i",icat));
      cout<<"Got Exp "<<endl;
      extbkg[totNdeg+1]->plotOn(plot,LineColor(icolbkg[1]), RooFit::Name(Form("pdf_Exp")));
      if(channel=="eeg" || channel=="eeg_mmg") legbkg->AddEntry(plot->findObject("pdf_Exp"),"Exp","l");
      if(channel=="mmg") legbkg->AddEntry(plot->findObject("pdf_Exp"),"Exp(2)","l");

      
      //plot the Pow
      cout<<"Doing Pow now "<<endl;
      extbkg[totNdeg+3] = ws->pdf(Form("pdf_Pow_bgr_cat%i",icat));
      extbkg[totNdeg+3]->plotOn(plot,LineColor(icolbkg[2]), RooFit::Name(Form("pdf_Pow")));
      legbkg->AddEntry(plot->findObject("pdf_Pow"),"Pow","l");
      

      //plot the Laurent
      cout<<"Doing Laurent now "<<endl;
      extbkg[totNdeg+4] = ws->pdf(Form("pdf_Laurent_bgr_cat%i",icat));
      extbkg[totNdeg+4]->plotOn(plot,LineColor(icolbkg[3]), RooFit::Name(Form("pdf_Laurent"))); ///commented on 22jan, 2017
      legbkg->AddEntry(plot->findObject("pdf_Laurent"),"Laurent","l");


      ////Now plot the extra Bernsteins
      for(int ideg=1; ideg<=5; ideg++){
	cout<<"plotting Other Berns now "<<endl;
	//plot the Exp
	//extbkg[totNdeg+ideg+4] = ws->pdf(Form("pdf_Bern%d_bgr_cat%i",ideg,icat));
	//pdf_Bern5_bgr_cat1 
	cout<<"looking for "<<Form("pdf_Bern%d_bgr_cat%i",ideg,icat)<<endl;
	
	//if(!ws->FindObject(Form("pdf_Bern%d_bgr_cat%i",ideg,icat)) ) continue;
	
	extbkg[totNdeg+ideg+4] = ws->pdf(Form("pdf_Bern%d_bgr_cat%i",ideg,icat));
	if(!extbkg[totNdeg+ideg+4]) continue;
	cout<<"Got extra Bern "<<ideg<<endl;
	
	extbkg[totNdeg+ideg+4]->plotOn(plot,LineColor(icolbkg[ideg+3]), RooFit::Name(Form("pdf_Bern%d",ideg)));
	legbkg->AddEntry(plot->findObject(Form("pdf_Bern%d",ideg)),Form("Bern(%d)",ideg),"l");
      }

      
      ///plot the Sech * Bern
      //pdf_SechBern2_bgr_cat
      /*
      int ibkg;
      for(int jbkg=2;jbkg<=totNdeg; jbkg++){
      //for(int jbkg=2;jbkg<=5; jbkg++){

	ibkg = totNdeg + 4 + jbkg;
	
	cout<<"plotting ibkg "<<ibkg<<endl;
	//pdf_2_bgr_cat1
	extbkg[ibkg] = ws->pdf(Form("pdf_SechBern%d_bgr_cat%i",jbkg,icat));
	cout<<"got the bkg "<<endl;
	//pdf_bgr_cat1_fitresult_2
	RooFitResult *rbkg = (RooFitResult *)ws->genobj(Form("pdf_bgr_cat%i_fitresult_SechBern%i",icat,jbkg));
	if(rbkg==NULL) cout<<"Null pointer"<<endl;
	//rbkg->Dump();
	
	
	rbkg->Print();
	rbkg->correlationMatrix().Print();
	
	extbkg[ibkg]->plotOn(plot,LineColor(icolbkg[ibkg-2]), RooFit::Name(Form("pdf_SechBern%d",jbkg)));

	legbkg->AddEntry(plot->findObject(Form("pdf_SechBern%d",jbkg)),Form("Sech * Bern(%d)",jbkg),"l");
	
	cout<<"done plotting bkg "<<endl;
      }//for(int ibkg=2;ibkg<=totNdeg; ibkg++)
      ///end of plotting sech * bern
      */
      /*
      cout<<"ibkg "<<ibkg<<endl;
      cout<<"2*totNdeg+5 "<<2*totNdeg+5<<endl;
      //plot the Sech * Pow
      cout<<"Doing Sech * Pow now "<<endl;
      extbkg[2*totNdeg+5] = ws->pdf(Form("pdf_SechPow_bgr_cat%i",icat));
      extbkg[2*totNdeg+5]->plotOn(plot,LineColor(icolbkg[2*totNdeg+5]), RooFit::Name(Form("pdf_SechPow")));
      legbkg->AddEntry(plot->findObject("pdf_SechPow"),"Sech(x)Pow","l");
      */

      ///try using TGraphAsymtricalErro
      TGraphAsymmErrors onesigma;
      TGraphAsymmErrors twosigma;

      //getErroBands(onesigma, twosigma, mass, fit, tmpCurve, data, testFrame, year, lep)
	
      
      char *filename = new char[100];
      char dirName[100] = "plots";
      
      TCanvas *c = setTCanvasNicev1("ccat");
      //c->Divide(1,2);

      gPad->SetLeftMargin(0.15);
      plot->GetYaxis()->SetTitleOffset(1.6);
      //plot->SetMinimum(0);
      if(blind) plot->SetMinimum(0.0001);

      double max = plot->GetMaximum();
      plot->SetMaximum(max*1.5);
      //plot->SetMinimum(0);
      if(blind) plot->SetMinimum(0.0001);

      plot->GetXaxis()->SetTitle("M_{Z#gamma} [GeV]");

      //h->Draw();
      plot->Draw();
      //h->Draw("same");
      
      int iPeriod = 4;
      int iPos = 11;
      //int iPos = 10;
      //int iPos = 40;
      CMS_lumi( c, iPeriod, iPos );
      c->Update();
      c->Modified();
      c->Update();


      legbkg->Draw();
      //cout<<"saving .gif"<<endl;
      //sprintf(filename,"%s/%s.gif",dirName, (channel+"_cat"+scat).c_str());
      cout<<"saving .pdf"<<endl;
      sprintf(filename,"%s/%s.pdf",dirName, (channel+"_cat"+scat).c_str());
      c->Print(filename);
      cout<<"saving .png"<<endl;
      sprintf(filename,"%s/%s.png",dirName, (channel+"_cat"+scat).c_str());
      c->Print(filename);
      
   
      //ggF
      RooPlot *plotsig[3];
      isig = 0;
      for(int imass=120; imass<=130; imass=imass+5){
	plotsig[isig] = mass->frame();
	cout<<"signal is "<<Form("signaldata_ggH_%d_cat%d",imass,icat) <<endl;
	RooAbsData *sig = ws->data( Form("signaldata_ggH_%d_cat%d",imass,icat) );
	RooAbsPdf *sigpdf = ws->pdf( Form("pdf_sig_ggH_%d_cat%d",imass,icat) );
	
	cout<<"Mass is "<<imass<<endl;
	cout<<"========printing info about the signal PDF========="<<endl;
	sigpdf->Print();
	cout<<"========printed info about the signal PDF========="<<endl;
	
	sig->plotOn(plotsig[isig]);
	sigpdf->plotOn(plotsig[isig]);

	RooAbsPdf *subpdf1 = ws->pdf(Form("subpdf1_sig_ggH_%d_cat%d",imass,icat));
	RooAbsPdf *subpdf2 = ws->pdf(Form("subpdf2_sig_ggH_%d_cat%d",imass,icat));
	//sigpdf->plotOn(plotsig[isig],Components(Form("subpdf1_sig_ggH_%d_cat%d",imass,icat)),LineColor(kRed),LineStyle(kDotted)) ;
	//sigpdf->plotOn(plotsig[isig],Components(Form("subpdf2_sig_ggH_%d_cat%d",imass,icat)),LineColor(kGreen),LineStyle(kDotted)) ;

	sigpdf->plotOn(plotsig[isig],Components(Form("subpdf1_sig_ggH_%d_cat%d",imass,icat)),LineColor(kRed),LineStyle(9)) ;
	sigpdf->plotOn(plotsig[isig],Components(Form("subpdf2_sig_ggH_%d_cat%d",imass,icat)),LineColor(kGreen),LineStyle(9)) ;
	
	///save now
	
	TCanvas *c = setTCanvasNicev1("csigggF");
	c->SetLogy();
	//plot->SetMinimum(0);
	plotsig[isig]->Draw();
	c->Print(Form("plots/sig%s_ggF_mass%d_cat%d.gif",channel.c_str(),imass,icat));
	c->Print(Form("plots/sig%s_ggF_mass%d_cat%d.pdf",channel.c_str(),imass,icat));
	
	c->Print(Form("plots/sig%s_ggF_mass%d_cat%d.C",channel.c_str(),imass,icat));
	
	isig++;
      }

    
      
    }//for(int icat=0; icat<4; icat++)

   ///plot the parameter plots here

   cout<<"Now plotting the parameters from different "<<endl;
   

   
   
}




//double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name){
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

