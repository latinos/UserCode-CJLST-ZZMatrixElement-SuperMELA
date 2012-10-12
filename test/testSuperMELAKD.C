#include <Riostream.h>
#include <string>

#include "ZZMatrixElement/SuperMELA/interface/SuperMELA.h"
#include "ZZMatrixElement/MELA/interface/PseudoMELA.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void testSuperMELAKD(){
  //load MELA and SuperMELA libraries before compililng this
  //.x loadMELA.C

  string str_sqrts="8TeV";//"7TeV";//"8TeV"
  const int nSamples=1;//10 for 8 TeV, 9 for 7TeV

  string chan[3]={"4mu","4e","2e2mu"};
  //  string files[10]={"HZZ4lTree_ZZTo4mu","HZZ4lTree_ZZTo4tau","HZZ4lTree_ZZTo4e","HZZ4lTree_ZZTo2e2mu","HZZ4lTree_ZZTo2e2tau","HZZ4lTree_ZZTo2mu2tau","HZZ4lTree_ggZZ2l2l","HZZ4lTree_ggZZ4l","HZZ4lTree_H125","HZZ4lTree_H126"};
  //  string files[1]={"HZZ4lTree_jhuPseH125"};
  string files[2]={"HZZ4lTree_H125","HZZ4lTree_H126"};

  int ich=0;//do 4mu
  int ifile=0; //do SM Higgs 125 GeV
 

  string dirName="/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_31082012/JHU_"+str_sqrts+"/"+chan[ich]+"/";
  cout<<"\n----------\nProcessing "<<files[ifile].c_str()<<"  "<<chan[ich].c_str()<<endl;
  string fileName=dirName+files[0]+".root";
  TFile *fIn=new TFile(fileName.c_str(),"READ");
  TTree *sigTree=(TTree*)fIn->Get("SelectedTree");
  
  //string fileName="./HZZ4lTree_H125_withDiscriminants.root";
  //   string fileName="./HZZ4lTree_ZZTo4mu_withDiscriminants.root";
  // TFile *fIn=new TFile(fileName.c_str(),"READ");
  // TTree *sigTree=(TTree*)fIn->Get("angles");
  


 float mzz,melapsig,melapbkg;
 float oldSMD,oldSMDPsig,oldSMDPbkg;
 float m1,m2,hs,h1,h2,phi,phi1,oldD,w,w_noxsec;//,pt4l,Y4l
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1);
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("helphi",&phi);
  sigTree->SetBranchAddress("phistarZ1",&phi1);
  sigTree->SetBranchAddress("ZZLD",&oldD);
  sigTree->SetBranchAddress("MC_weight",&w);
  sigTree->SetBranchAddress("MC_weight_noxsec",&w_noxsec);
  //  sigTree->SetBranchAddress("mela_psig",&melapsig);
  // sigTree->SetBranchAddress("mela_pbkg",&melapbkg);
  // sigTree->SetBranchAddress("supermelaLD",&oldSMD);
  // sigTree->SetBranchAddress("supermela_psig",&oldSMDPsig);
  //sigTree->SetBranchAddress("supermela_pbkg",&oldSMDPbkg);

  double smd,mela,psig,pbkg;
  float psmela,psigps,pbkgps;

 string outFileName=dirName+files[ifile]+"_withSMD.root";
 TFile *fout=new TFile(outFileName.c_str(),"RECREATE");
 TTree *outTree=new TTree("SelectedTree","SelectedTree");
 outTree->Branch("Z2Mass",&m2,"Z2Mass/F");
 outTree->Branch("Z1Mass",&m1,"Z1Mass/F");
 outTree->Branch("ZZMass",&mzz,"ZZMass/F");
 outTree->Branch("costhetastar",&hs,"costhetastar/F");
 outTree->Branch("helcosthetaZ1",&h1,"helcosthetaZ1/F");
 outTree->Branch("helcosthetaZ2",&h2,"helcosthetaZ2/F");
 outTree->Branch("helphi",&phi,"helphi/F");
 outTree->Branch("phistarZ1",&phi1,"phistarZ1/F");
 outTree->Branch("ZZLD",&mela,"ZZLD/D");
 outTree->Branch("superLD",&smd,"superLD/D");
 outTree->Branch("pseudoLD",&psmela,"pseudoLD/F");
 outTree->Branch("MC_weight",&w,"MC_weight/F");
 outTree->Branch("MC_weight_noxsec",&w_noxsec,"MC_weight_noxsec/F");

 //PseudoMELA *mypsLD=new PseudoMELA();

 SuperMELA *mySMD=new SuperMELA(125.0,chan[ich],8);
 // // //you can also change channel and mH after having created the object SuperMELA
 // mySMD->SetDecayChannel(chan[1]);
 // mySMD->SetMH(126.0);
 mySMD->SetPathToCards("/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/spin/CMSSW_5_2_5/src/HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/SM_inputs_"+str_sqrts+"/");
 // // // by default, expect to receive pre-calculated Psig and Pbkg from standard MELA
 // // // set RecalculateMELA to true for recalculating it (you must pass the 4-vectors in that case)
 // mySMD->RecalculateMELA(true);
 //mySMD->SetVerbosity(true);

 mySMD->init();

 // // // example passing pre-calculated values of Psig and Pbkg 
 mySMD->computeKD(124.43, 0.8, 0.33,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (1):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 mySMD->computeKD(126.91, 0.52, 0.68,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (2):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 mySMD->computeKD(136.19, 0.8, 0.33,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (3):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 mySMD->computeKD(136.19, 0.21, 0.87,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (3):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 

 // // // example recalculating Mela
 mySMD->RecalculateMELA(true);
 mySMD->SetVerbosity(false);
 int i=1;//entry of input tree
 sigTree->GetEntry(i);
 if(mzz>180.0||mzz<100.0){
   cout<<"Calculation not available outside outside m4l range [100, 180]. m4l in input is "<<mzz<<endl;
 }
 else{

   // // // tell to SuperMELA the angles needed for recalculating MELA
   mySMD->SetDecayKinematics(m1,m2,hs,h1,h2,phi,phi1);
   mySMD->computeKD(mzz,false,  smd,mela,psig,pbkg);
   cout<<"Entry #"<<i<<" OLD-MELA="<<oldD<<"  New-MELA="<<mela<<"  SuperMELA="<<smd<<endl;
   
   // mySMD->computeKD(mzz, melapsig, melapbkg,   smd,mela,psig,pbkg);
   // mypsLD->computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,psmela,psigps,pbkgps);

   outTree->Fill();
   
 }


 fout->cd();
 outTree->Write();
 delete fout;



 delete  mySMD;
 //delete  mypsLD;

}//end main
