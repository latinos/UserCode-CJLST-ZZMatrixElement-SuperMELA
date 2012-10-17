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
#include "TRandom3.h"

void prodSuperMELAKD(){
  //load MELA and SuperMELA libraries
  //.x loadMela.C

  string str_sqrts="8TeV";//"7TeV";//"8TeV"
  const int nSamples=10;//10 for 8 TeV, 9 for 7TeV

  string chan[3]={"4mu","4e","2e2mu"};
    string files[10]={"HZZ4lTree_ZZTo4mu","HZZ4lTree_H125","HZZ4lTree_ZZTo4tau","HZZ4lTree_ZZTo4e","HZZ4lTree_ZZTo2e2mu","HZZ4lTree_ZZTo2e2tau","HZZ4lTree_ZZTo2mu2tau","HZZ4lTree_ggZZ2l2l","HZZ4lTree_ggZZ4l","HZZ4lTree_H126"};
    //string files[1]={"HZZ4lTree_jhuPseH125"};
  //string files[2]={"HZZ4lTree_H125","HZZ4lTree_H126"};

  TRandom3 *myR=new TRandom3(4887);

  for(int ich=0;ich<3;ich++){

      if(ich!=0)continue;

    string dirName="/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_31082012/PRODFSR_"+str_sqrts+"/"+chan[ich]+"/";
    //string dirName="/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_31082012/JHU_"+str_sqrts+"/"+chan[ich]+"/";



    // for(int ifile=0;ifile<nSamples;ifile++){
 for(int ifile=0;ifile<2;ifile++){

    cout<<"\n----------\nProcessing "<<files[ifile].c_str()<<"  "<<chan[ich].c_str()<<endl;

    bool isSignal=(files[ifile].find("H125")!=string::npos);




  string fileName=dirName+files[ifile]+".root";
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
  // sigTree->SetBranchAddress("mela_psig",&melapsig);
  // sigTree->SetBranchAddress("mela_pbkg",&melapbkg);
  // sigTree->SetBranchAddress("supermelaLD",&oldSMD);
  // sigTree->SetBranchAddress("supermela_psig",&oldSMDPsig);
  //sigTree->SetBranchAddress("supermela_pbkg",&oldSMDPbkg);

  double smd, mela,psig,pbkg, melapsigOut,melapbkgOut;
  double smdSyst1Up, smdSyst1Down, smdSyst2Up, smdSyst2Down, melaTmp,psigTmp,pbkgTmp;
  float psmela,psigps,pbkgps;

 string outFileName=dirName+files[ifile]+"_withSMDv2.root";
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
 outTree->Branch("ZZLD_PSig",&melapsigOut,"ZZLD_PSig/D");
 outTree->Branch("ZZLD_PBkg",&melapbkgOut,"ZZLD_PBkg/D");
 outTree->Branch("superLD",&smd,"superLD/D");
 outTree->Branch("pseudoLD",&psmela,"pseudoLD/F");
 outTree->Branch("MC_weight",&w,"MC_weight/F");
 outTree->Branch("MC_weight_noxsec",&w_noxsec,"MC_weight_noxsec/F");
 outTree->Branch("superLD_syst1Up",&smdSyst1Up,"superLD_syst1Up/D");
 outTree->Branch("superLD_syst1Down",&smdSyst1Down,"superLD_syt1Down/D");
 outTree->Branch("superLD_syst2Up",&smdSyst2Up,"superLD_syst2Up/D");
 outTree->Branch("superLD_syst2Down",&smdSyst2Down,"superLD_syt2Down/D");
 // outTree->Branch("ZZLD",&oldD,"");

 PseudoMELA *mypsLD=new PseudoMELA();

 SuperMELA *mySMD=new SuperMELA(125.0,chan[ich],8);
 // mySMD->SetDecayChannel(chan[1]);
 // mySMD->SetMH(126.0);
 mySMD->SetPathToCards("/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/spin/CMSSW_5_2_5/src/HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/SM_inputs_"+str_sqrts+"/");
 // mySMD->RecalculateMELA(true);
 //mySMD->SetVerbosity(true);
 mySMD->init();

 /*
 mySMD->computeKD(124.43, 0.8, 0.33,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (1):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 mySMD->computeKD(126.91, 0.52, 0.68,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (2):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 mySMD->computeKD(136.19, 0.8, 0.33,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (3):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 mySMD->computeKD(136.19, 0.21, 0.87,   smd,mela,psig,pbkg);
 cout<<"testKDSuperMELA (3):  SuperMELA="<<smd<<"  Psig="<<psig<<"   Pbkg="<<pbkg<<"  MELA="<<mela<<endl;
 */

 mySMD->RecalculateMELA(true);
 mySMD->SetVerbosity(false);

 double meanCB_err=mySMD->GetSigShapeSystematic("meanCB");
 double sigmaCB_err=mySMD->GetSigShapeSystematic("sigmaCB");
 cout<<"Signal shape syst factors are MeanSystErr="<<meanCB_err<<"  SigmaSystErr="<<sigmaCB_err<<endl;
 cout<<"Looping on test TTree"<<endl;
 int nDiff=0;
 int maxEntries=1000;//sigTree->GetEntries();//1000;//sigTree->GetEntries()
 TH1F *hSMD=new TH1F("hsmd"," SuperMELA",200,0.0,1.0);

 for(int i=0;i<maxEntries;i++){
   sigTree->GetEntry(i);
   if(i%1000==0)cout<<"Entry #"<<i<<endl;
   if(mzz>180.0||mzz<100.0)continue;




   /*****
    //for regular way, no recalculation of MELA
    melapsigOut=melapsig;
    melapbkgOut=melapbkg;
   mySMD->computeKD(mzz, melapsigOut, melapbkgOut, smd,mela,psig,pbkg);

   double mzzTmpSig=0.0, mzzTmpBkg=double(mzz);
   double melaTmp2;
   mzzTmpSig=double( mzz*(1.0+meanCB_err) );
   std::pair<double,double> mzzTmpPair = make_pair(mzzTmpSig, mzzTmpBkg);//Psig and Pbkg
   mySMD->computeKD(mzzTmpPair, melapsig, melapbkg, smdSyst1Up, melaTmp2,psigTmp,pbkgTmp);
 
   mzzTmpSig=double( mzz*(1.0-meanCB_err) );
   mzzTmpPair = make_pair(mzzTmpSig, mzzTmpBkg);//Psig and Pbkg
   mySMD->computeKD(mzzTmpPair, melapsig, melapbkg, smdSyst1Down, melaTmp2,psigTmp,pbkgTmp);
 
   mzzTmpSig=mzz* myR->Gaus(1.0,sigmaCB_err);
   mzzTmpPair = make_pair(mzzTmpSig,mzzTmpBkg);//Psig and Pbkg
   mySMD->computeKD(mzzTmpPair, melapsig, melapbkg, smdSyst2Up, melaTmp2,psigTmp,pbkgTmp);
   smdSyst2Down=smdSyst2Up;
   *****/


   //////for recalculating smd
   mySMD->RecalculateMELA(true);
   mySMD->SetDecayKinematics(m1,m2,hs,h1,h2,phi,phi1);
   mySMD->computeKD(mzz,false, smd,psig,pbkg,mela, melapsigOut, melapbkgOut);
   mypsLD->computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,psmela,psigps,pbkgps);   


   double mzzTmpSig=0.0, mzzTmpBkg=double(mzz);
   double melaTmp2;
   if(isSignal){//signal sample at MH=125 GeV
   mySMD->RecalculateMELA(false);
   mzzTmpSig=double( mzz*(1.0+meanCB_err) );
   if(mzzTmpSig>180.0 || mzzTmpSig<100){
     cout<<"Entry #"<< i<<"   ATTENTION: changing scale up has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
     mzzTmpSig=mzz;
   }
   mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst1Up, melaTmp2,psigTmp,pbkgTmp); 
   if( smdSyst1Up<0.0)cout<<"Entry #"<< i<<"   New SMD with scale Up is "<<smdSyst1Up<<endl;
   
   mzzTmpSig=double( mzz*(1.0-meanCB_err) );
   if(mzzTmpSig>180.0 || mzzTmpSig<100){
     cout<<"Entry #"<< i<<"   ATTENTION: changing scale down has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
     mzzTmpSig=mzz;
   }
   mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst1Down, melaTmp2,psigTmp,pbkgTmp);
   mzzTmpSig=mzz* myR->Gaus(1.0,sigmaCB_err);
   if(mzzTmpSig>180.0 || mzzTmpSig<100){
     cout<<"Entry #"<< i<<"   ATTENTION: smearing scale has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
     mzzTmpSig=mzz;
   }
   mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst2Up, melaTmp2,psigTmp,pbkgTmp);
   smdSyst2Down=smdSyst2Up;
   }//end if signal smaple at 125
   else{
     smdSyst1Up=smd;
     smdSyst1Down=smd;
     smdSyst2Up=smd;
     smdSyst2Down=smd;
   }


   // cout<<"Entry #"<<i<<" OLD-MELA="<<oldD<<"  New-MELA="<<mela<<"  SuperMELA="<<smd<<endl;

   // cout<<"Entry #"<<i<<"     SuperMELA="<<smd<<"  SuperMELA-PSig="<<psig<<"   SuperMELA-Pbkg="<<pbkg<<endl;
   // cout<<"Entry #"<<i<<"(OLD)SuperMELA="<<oldSMD<<"  SuperMELA-PSig="<<oldSMDPsig<<"   SuperMELA-Pbkg="<<oldSMDPbkg<<endl;
   //   if(fabs(mela-oldD)/oldD>0.10){//detailed dump for understanding differences
   //  cout<<"Entry #"<<i<<" mZZ="<<mzz<<"  OLD-MELA="<<melapsig/(melapsig+melapbkg)<<" OLD-MELA-PSig="<<melapsig<<"   OLD-MELA-PBkg="<< melapbkg<<" M1="<<m1<<"   M2="<<m2<<"  CosTStar="<<hs<<"  CosTheta1="<<h1<<"CosTheya2="<<h2<<"  Phi="<<phi<<"  PhiStar="<<phi1<<std::endl;
   //  nDiff++; 
   //   }

 
   hSMD->Fill(smd);
   outTree->Fill();
 }//end loop on events

 cout<<"Finished to loop on "<<maxEntries<<" events. "<<nDiff<<" have the new SMD different by >10% respect to what is stored i nthe original TTree"<<std::endl;


 fout->cd();
 outTree->Write();
 delete fout;
 TCanvas *c2=new TCanvas("cc2","CC2",1000,1000);
 c2->cd();
 hSMD->SetXTitle("SuperMELA");
 hSMD->SetYTitle("# events");
 hSMD->Draw();
 c2->SaveAs(("can_"+files[ifile]+"_SuperMELA.root").c_str());


 delete  mySMD;
 delete  mypsLD;
  }//end loop on ifile (samples to process for this channel)
  }//end loop on ich (channels)

}//end main
