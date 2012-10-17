#ifndef MELA_SuperMela_proto_h
#define MELA_SuperMela_proto_h

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "TString.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "TFile.h"

class SuperMela_proto{

 public:

  // methods
  SuperMela_proto(float mH=125.,TString workspaceInputDir="../test/",
		  TString channel="4mu",int LHCsqrts=8,bool verbose_=false);
  ~SuperMela_proto(){
    delete wsp;
    delete sigNorm;
    delete bkgNorm;
  };

  void init();
  
  void computeKD(float m4l,float PSigMelaIn,float PBkgMelaIn, //in
		 float &superMELA,float &Psig, float &Pbkg);  //out
  
  void setDecayKinematics(float m1,float m2,
			  float hs,float h1,float h2, 
			  float phi,float phi1);
  
  void setMH(float mH);

  void setDecayChannel(TString channel);

  void setSigCBmeanSyst(float syst);
  void setSigCBsigmaSyst(float syst);
  void setVerbose(bool yesNo){verbose=yesNo;};

 private:

  // data members
  RooWorkspace* wsp;
  TString chan;
  TString wkspDir;
  int sqrts;
  float m1_,m2_,h1_,h2_,hs_,phi_,phi1_;
  RooAbsReal *sigNorm,*bkgNorm;
  bool verbose;
  
};

#endif
