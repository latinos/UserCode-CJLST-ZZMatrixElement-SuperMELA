#include "SuperMela_proto.h"


SuperMela_proto::SuperMela_proto(float mH,TString workspaceInputDir,
				 TString channel,int LHCsqrts, bool verbose_){

  verbose=verbose_;
  if(verbose)
    cout << "SuperMela_proto::SuperMela_proto" << endl;

  wkspDir=workspaceInputDir;
  chan=channel;
  sqrts=LHCsqrts;

  init();

  if(!wsp){
    cout << "SuperMela_proto::SuperMela_proto - error in init(), workspace currently NULL" << endl;
    return;
  }

  if(verbose) cout << "setting mH to " << mH << endl;

  wsp->var("MH")->setVal(mH);

}

void SuperMela_proto::init(){

  if(verbose)
    cout << "SuperMela_proto::init" << endl;

  TString workspacePath=wkspDir+"hzz4l_"+chan+"S_";
  workspacePath+=sqrts;
  workspacePath+="TeV.input.root";

  if(verbose) cout << "reading workspace from " << workspacePath << endl;

  TFile f(workspacePath);

  wsp = (RooWorkspace*) f.Get("w");

  if(!wsp){
    cout << "SuperMela_proto::init - error, wsp is NULL" << endl;
    return;
  }

  sigNorm = wsp->pdf("signalCB_ggH")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));
  bkgNorm = wsp->pdf("bkg_qqzzTmp")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));
  
  if(verbose) cout << "sigNorm " << sigNorm << " bkgNorm " << bkgNorm << endl;

}

void SuperMela_proto::setMH(float mH){

  if(verbose)
    cout << "SuperMela_proto::setMH" << endl;

  if(!wsp){
    cout << "SuperMela_proto::setMH - error, wsp is NULL" << endl;
    return;
  }
  
  wsp->var("MH")->setVal(mH);
  sigNorm = wsp->pdf("signalCB_ggH")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));
  bkgNorm = wsp->pdf("bkg_qqzzTmp")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));
 
}

void SuperMela_proto::setDecayChannel(TString channel){

  if(verbose)
    cout << "SuperMela_proto::setDecayChannel" << endl;

  chan=channel;
  init();

}

void SuperMela_proto::setDecayKinematics(float m1,float m2,
				    float hs,float h1,float h2,
				    float phi,float phi1){

  if(verbose)
    cout << "SuperMela_proto::setDecayKinematics" << endl;

  m1_=m1;
  m2_=m2;
  hs_=hs;
  h1_=h1;
  h2_=h2;
  phi_=phi;
  phi1_=phi1;

}

void SuperMela_proto::computeKD(float m4l,float PSigMelaIn,float PBkgMelaIn, //in
			   float &superMELA,float &Psig, float &Pbkg){  //out

  if(verbose)
    cout << "SuperMela_proto::computeKD" << endl;

  if(!wsp){
    cout << "SuperMela_proto::computeKD - error, wsp is NULL" << endl;
    superMELA=-99; 
    Psig=-99;
    Pbkg=-99;
    return;
  }
  
  wsp->var("CMS_zz4l_mass")->setVal(m4l);
  Psig = PSigMelaIn*wsp->pdf("signalCB_ggH")->getVal()/sigNorm->getVal();
  Pbkg = PBkgMelaIn*wsp->pdf("bkg_qqzzTmp")->getVal()/bkgNorm->getVal();
  superMELA=Psig/(Psig+Pbkg);

}

void SuperMela_proto::setSigCBmeanSyst(float syst){

  if(!wsp){
    cout << "SuperMela_proto::setSigCBmeanSyst - error, wsp is NULL" << endl;
    return;
  }

  if(verbose){
    cout << "SuperMela_proto::setSigCBmeanSyst" << endl;
    cout << "CMS_zz4l_mean_m_sig: " << wsp->var("CMS_zz4l_mean_m_sig") << endl;
    cout << "CMS_zz4l_mean_e_sig: " << wsp->var("CMS_zz4l_mean_e_sig") << endl;
  }

  if(wsp->var("CMS_zz4l_mean_m_sig"))
     wsp->var("CMS_zz4l_mean_m_sig")->setVal(syst);
  if(wsp->var("CMS_zz4l_mean_e_sig"))
     wsp->var("CMS_zz4l_mean_e_sig")->setVal(syst);

  sigNorm = wsp->pdf("signalCB_ggH")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));
  bkgNorm = wsp->pdf("bkg_qqzzTmp")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));

}

void SuperMela_proto::setSigCBsigmaSyst(float syst){

  if(!wsp){
    cout << "SuperMela_proto::setSigCBsigmaSyst - error, wsp is NULL" << endl;
    return;
  }
  if(verbose){
    cout << "SuperMela_proto::setSigCBsigmaSyst" << endl;
    cout << "CMS_zz4l_sigma_m_sig: " << wsp->var("CMS_zz4l_sigma_m_sig") << endl;
    cout << "CMS_zz4l_sigma_e_sig: " << wsp->var("CMS_zz4l_sigma_e_sig") << endl;
  }

  if(wsp->var("CMS_zz4l_sigma_m_sig"))
     wsp->var("CMS_zz4l_sigma_m_sig")->setVal(syst);
  if(wsp->var("CMS_zz4l_sigma_e_sig"))
     wsp->var("CMS_zz4l_sigma_e_sig")->setVal(syst);

  sigNorm = wsp->pdf("signalCB_ggH")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));
  bkgNorm = wsp->pdf("bkg_qqzzTmp")->createIntegral(RooArgSet(*wsp->var("CMS_zz4l_mass")));

}
 
