#include <iostream>

//const int nprocess = 2;
const int nprocess = 6;
const int nsig = 3;
  const int nmass = 30;

class statAn {



 public:


  statAn(string _model, int _poldeg) ;
  ~statAn();
  
  void simplemkdatacard(const char* channel = "eeg", int cat = 1);
  
  void mybkgfit(double xmin, double xmax, int icat);

  RooAbsData* fill_events(TH1* hMass, const char* filename, double trigEff, int cat = -1, bool usewei = true);

  RooAbsData* fill_events(TH1* hMass, const char* filename1, const char* filename2, double trigEff1, double trigEff2, int cat=-1, bool usewei=true);

  double normfactor_SMHiggs_8TeV(const char* channel, int p, int mass) ;

  void mkdatacard(const char* channel = "eeg", int cat = 1);
  void smallcopymkdatacard(const char* channel = "eeg", int cat = 1);
  void mysigfit(string model,double xmin, double xmax, bool nuisance, int np, int nm, int icat, string schannel);

  void mysigfunc(string model);
  
  void siginterpolate(string model, int np, int nm, int nk, double a, double b, bool nuisance, int icat);

  void setSigParConst(string model, bool nuisance);
  
  void deleteSigPar(bool nuisance);

  std::string getbkgfitParName(int ipar);

  void getIDsys(const char* filename, double trigEff, double &sys_lep, double &sys_pho, int cat=-1, bool usewei=true);

  void biasStudymkdatacards(const char* channel, int cat);
  int getNbkgPar();
  
  RooAbsData* dataObs;
  RooAbsPdf* bgrfit;
  RooFitResult *rfit;
  RooRealVar  *fX;
  
  RooRealVar  *pol[64];

  RooFormulaVar  *sqpol[64];

  RooArgList *coeflist;
  RooRealVar*   fPar[64];
  RooRealVar*   fPar_sig[64];
  
  RooRealVar *mean;
  RooRealVar *sigma;
  RooRealVar *stepval;

  RooRealVar *mu;
  RooRealVar *lambda;
  RooRealVar *alphabkg;
  RooRealVar *betabkg;


  

  ////signal paramters

  RooRealVar *mean1;
  RooRealVar *sigma1;
  RooRealVar *alpha;
  RooRealVar *power;
  //RooRealVar *delta21;
  RooRealVar *s21;
  RooRealVar *frac;
  RooAbsReal *sigmean;
  RooAbsReal *sigsigma;

  //RooAbsReal *mean2;
  RooAbsReal *sigma2;
  
  RooRealVar* scale;
  RooRealVar* resol;

  RooRealVar* stepbkg;
  
  RooCBShape*  pdf1;
  RooGaussian* pdf2;
  RooAddPdf *sigfit1;

  RooAbsData* rooMassSig[nprocess][nmass];
  RooAbsPdf* sigfit[nprocess][nmass]; 
  RooFitResult *rsig[nprocess][nmass];
  
  

  string model;
  int poldeg;
  
  

};


statAn::statAn(string _model,int _poldeg):model(_model), poldeg(_poldeg)
{
  
}



  
statAn::~statAn()
{}


