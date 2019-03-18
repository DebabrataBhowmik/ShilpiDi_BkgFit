#include <iostream>

//const int nprocess = 6;
const int nprocess = 5;
const int nsig = 3;
  const int nmass = 30;

class statAn {



 public:


  statAn(string _model, int _poldeg) ;
  ~statAn();
  
  void simplemkdatacard(const char* channel = "eeg", int cat = 1);
  void simplemkdatacard_LT(const char* channel = "eeg", int cat = 1);
  
  void mybkgfit(double xmin, double xmax, int icat, string schannel);

  RooAbsData* fill_events(TH1* hMass, const char* filename, double trigEff, int cat = -1, bool usewei = true);


  RooAbsData* fill_events(TH1* hMass, const char* filename1, const char* filename2, double trigEff1, double trigEff2, int cat=-1, bool usewei=true);

  RooAbsData* fill_events(TH1* hMass, const char* filename1, const char* filename2, const char* filename3, const char* filename4, double trigEff1, double trigEff2, int cat=-1, bool usewei=true);

  double normfactor_SMHiggs_8TeV(const char* channel, int p, int mass) ;
  //  double ngen_Tot(const char* channel, int p, int mass) ;

  void mkdatacard(const char* channel = "eeg", int cat = 1);
  void smallcopymkdatacard(const char* channel = "eeg", int cat = 1);
  void mysigfit(string model,double xmin, double xmax, bool nuisance, int np, int nm, int icat, string schannel);

  void mysigfunc(string model);
  
  void siginterpolate(string model, int np, int nm, int nk, double a, double b, bool nuisance, int icat);

  void setSigParConst(string model, bool nuisance);
  
  void deleteSigPar(bool nuisance);

  std::string getbkgfitParName(int ipar);

  void getIDsys(const char* filename, double trigEff, double &sys_lep, double &sys_pho, int cat=-1, bool usewei=true);

  void getIDsys(const char* filename1, const char* filename2, double trigEff1 , double trigEff2, double &sys_lep, double &sys_pho, int cat=-1, bool usewei=true);
  
  void getPUsys(const char* filename, double trigEff, double &sys_pu, int cat=-1, bool usewei=true);
  
  void getPUsys(const char* filename1, const char* filename2, double trigEff1 , double trigEff2, double &sys_pu, int cat=-1, bool usewei=true);


  void getHLTsys(const char* filename, double trigEff, double &sys_lep, int cat=-1, bool usewei=true);

  void getHLTsys(const char* filename1, const char* filename2, double trigEff1 , double trigEff2, double &sys_lep, int cat=-1, bool usewei=true);
  
  double getNLL(RooRealVar *mhzg, RooAbsPdf *pdf, RooAbsData *data, double normVal=-1., double massRangeLow=-1., double massRangeHigh=-1.);
  
  double guessNew(RooRealVar *mhzg, RooAbsPdf *pdf, RooAbsData *data, double bestPoint, double nllBest, double boundary, double massRangeLow, double massRangeHigh, double crossing, double tolerance);
  
  int getNbkgPar();
  
  RooAbsData* dataObs;
  RooAbsPdf* bgrfit;
  RooFitResult *rfit;
  RooRealVar  *fX;
  


  RooFormulaVar  *sqpol[64];
  
  RooRealVar  *pol[64];
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
  
  TH1D *hMasssig[nprocess][nmass];

  string model;
  int poldeg;
  
  

};


statAn::statAn(string _model,int _poldeg):model(_model), poldeg(_poldeg)
{
  
}



  
statAn::~statAn()
{}


