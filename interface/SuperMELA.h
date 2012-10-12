#ifndef MELA_SuperMela_h
#define MELA_SuperMela_h


#include <Riostream.h>
#include <string>
#include <fstream>

#include "TLorentzVector.h"
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"



class SuperMELA {

 public:
  SuperMELA(double mH=120,string channel="4mu",int LHCsqrts=8);
  ~SuperMELA();
  void init();
  void computeKD(double m4l,double PSigMelaIn,double PBkgMelaIn,   //in
		 double &superMELA,double &MELA,double &Psig,double &Pbkg);//out

//use the following only if you recalc MELA on the fly
  void computeKD(double m4l,bool use4vectors,
		 double &superMELA,double &MELA,double &Psig,double &Pbkg);//out


  //setters
  void  RecalculateMELA(bool recoMELA=true){recalculateMELA_=recoMELA;init();}//default is to take MELA as external input
  void SetVerbosity(bool verb=true){verbose_=verb;}
  void SetDecayChannel(string myChan);
  void SetMH(double myMH){
    mHVal_=myMH;
    mH_rrv_->setVal(mHVal_);
    if(verbose_)std::cout<<"Setting MH to "<<mHVal_<<std::endl;
    init();
  }

  void SetMELAProbabilities(double myPsig,double myPbkg){

    mela_psig_=myPsig;
    mela_pbkg_=myPbkg;
    if(verbose_)std::cout<<"Setting MELA Psig -> "<<mela_psig_<<"   Pbkg -> "<<mela_pbkg_<<std::endl;

  }

  void SetDecayKinematics(TLorentzVector Z1_lept1, int Z1_lept1Id,
			  TLorentzVector Z1_lept2, int Z1_lept2Id,
			  TLorentzVector Z2_lept1, int Z2_lept1Id,
			  TLorentzVector Z2_lept2, int Z2_lept2Id){
    Z1_lept1_=Z1_lept1;
    Z1_lept2_=Z1_lept2;
    Z2_lept1_=Z2_lept1;
    Z2_lept2_=Z2_lept2;

    Z1_lept1Id_=Z1_lept1Id;
    Z1_lept2Id_=Z1_lept2Id;
    Z2_lept1Id_=Z2_lept1Id;
    Z2_lept2Id_=Z2_lept2Id;
  }
  void SetDecayKinematics(float m1,float m2,float hs,float h1,float h2,float phi,float phistar1){
    m1_=m1;     m2_=m2;
    hs_=hs;
    h1_=h1;     h2_=h2;
    phi_=phi;     phistar1_=phistar1;
  }

  void SetPathToCards(string dirToCards){ pathToCards_=dirToCards;
    if(verbose_)std::cout<<"New path to cards is "<<pathToCards_.c_str()<<std::endl;}


 private:

  void readSigParsFromFile(string &str_mean_CB,string &str_sigma_CB ,string &str_n_CB ,string &str_alpha_CB);
  void readBkgParsFromFile(std::vector<double> &apars );
  void calc_mZZ_range(const double mHVal,double &low_M,double &high_M);
  std::pair<double,double>  superMelaLikelihoodDiscriminant (double m4l,double melaPsig,double melaPbkg);
  bool checkChannel();
  ///data members
  float mela_psig_,mela_pbkg_;
  double mHVal_;
  double sqrts_;
  double lowMH_,highMH_;
  string strChan_; int ch_;
  bool verbose_;
  string pathToCards_;
  //signal m4l shape
  RooRealVar *m4l_rrv_;//this one is the variable!
  RooRealVar *mH_rrv_;//this one is a fixed param !
  RooRealVar *mean_dummy_,*sigma_dummy_,*alpha_dummy_,*n_dummy_;
  RooFormulaVar *n_CB_, *alpha_CB_,*mean_CB_,*sigma_CB_,*meanTOT_CB_;//,*gamma_BW_;
  RooCBShape *sig_CB_;
  double norm_sig_CB_;

  //qqZZ background m4l shape
  RooRealVar *a0_qqZZ_, *a1_qqZZ_, *a2_qqZZ_, *a3_qqZZ_,
             *a4_qqZZ_, *a5_qqZZ_, *a6_qqZZ_, *a7_qqZZ_,
             *a8_qqZZ_, *a9_qqZZ_, *a10_qqZZ_,*a11_qqZZ_,
             *a12_qqZZ_,   *a13_qqZZ_  ;
  RooqqZZPdf_v2 *qqZZ_pdf_;       
  double norm_bkg_qqZZ_;

  //needed for recalculating MELA
  bool recalculateMELA_;
  //  bool use4vectors_;
  Mela *melaProd_;
  TLorentzVector Z1_lept1_,Z1_lept2_,Z2_lept1_,Z2_lept2_;
  int Z1_lept1Id_, Z1_lept2Id_, Z2_lept1Id_, Z2_lept2Id_;
  float m1_,m2_,hs_,h1_,h2_,phi_,phistar1_;
};


#endif
