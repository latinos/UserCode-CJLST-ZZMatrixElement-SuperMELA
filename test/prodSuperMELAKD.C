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
//#include <boost/filesystem.hpp> //it would be nice to have it properly installed on afs...

bool checkFileExistence(string file){
 
  //  if ( !boost::filesystem::exists( file.c_str() ) ) return false;
  // return true;
  
  ifstream myfile(file.c_str());
  if(!myfile.good()) return false;
  return true;
}

void prodSuperMELAKD(){
  //load MELA and SuperMELA libraries
  //.x loadMela.C

  string str_sqrts="8TeV";//"7TeV";//"8TeV"
  string chan[3]={"4mu","4e","2e2mu"};

  const string genType= "JHU";//"PRODFSR" "JHU"
  const int nSamples=1;//10 for 8 TeV, 9 for 7TeV
  //string files[10]={"HZZ4lTree_ZZTo4mu","HZZ4lTree_H125","HZZ4lTree_ZZTo4tau","HZZ4lTree_ZZTo4e","HZZ4lTree_ZZTo2e2mu","HZZ4lTree_ZZTo2e2tau","HZZ4lTree_ZZTo2mu2tau","HZZ4lTree_ggZZ2l2l","HZZ4lTree_ggZZ4l","HZZ4lTree_H126"};
  string files[1]={"HZZ4lTree_jhuPseH125"};


  TRandom3 *myR=new TRandom3(4887);

  for(int ich=0;ich<3;ich++){

    // if(ich!=2)continue;
 

    string chanDir=chan[ich];
    if(chanDir=="2e2mu")chanDir="2mu2e";

    string dirSqrtS=(str_sqrts=="7TeV"? genType : genType+"_8TeV");
    string dirName="root://lxcms02//data/Higgs/rootuplesOut/191012/"+dirSqrtS+"/"+chanDir+"/";
    string outDirName="/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/Trees_191012_sf0.5/"+genType+"_"+str_sqrts+"/"+chan[ich]+"/";
    for(int ifile=0;ifile<nSamples;ifile++){     
      // if(ifile!=1)continue;

    cout<<"\n----------\nProcessing "<<files[ifile].c_str()<<"  "<<chan[ich].c_str()<<endl;

    bool isSignal=(files[ifile].find("H125")!=string::npos)||(files[ifile].find("H126")!=string::npos);

  string fileName=dirName+files[ifile]+".root";

  //  if(!checkFileExistence(fileName)){
  //  cout<<"\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!\n\nWARNING : File "<<fileName.c_str()<<" NOT FOUND ! Skipping it"<<"\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n"<<endl;
  //  continue;
  // }
  TFile *fIn=TFile::Open(fileName.c_str(),"READ");
  TTree *sigTree=(TTree*)fIn->Get("SelectedTree");
  
  //string fileName="./HZZ4lTree_H125_withDiscriminants.root";
  //   string fileName="./HZZ4lTree_ZZTo4mu_withDiscriminants.root";
  // TFile *fIn=new TFile(fileName.c_str(),"READ");
  // TTree *sigTree=(TTree*)fIn->Get("angles");
  


 float mzz,melapsig,melapbkg;
 float oldSMD,oldSMDPsig,oldSMDPbkg;
 float m1,m2,hs,h1,h2,phi,phi1,oldD,w,w_noxsec;//,pt4l,Y4l
 float oldPSD,oldGravD;
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("Z1Mass",&m1);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1);
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("helphi",&phi);
  sigTree->SetBranchAddress("phistarZ1",&phi1);
  sigTree->SetBranchAddress("ZZLD",&oldD);
  sigTree->SetBranchAddress("ZZpseudoLD",&oldPSD);
  sigTree->SetBranchAddress("ZZgravLD",&oldGravD);
  sigTree->SetBranchAddress("MC_weight",&w);
  sigTree->SetBranchAddress("MC_weight_noxsec",&w_noxsec);
  sigTree->SetBranchAddress("ZZLD_PSig",&melapsig);
  sigTree->SetBranchAddress("ZZLD_PBkg",&melapbkg);
  // sigTree->SetBranchAddress("supermelaLD",&oldSMD);
  // sigTree->SetBranchAddress("supermela_psig",&oldSMDPsig);
  //sigTree->SetBranchAddress("supermela_pbkg",&oldSMDPbkg);

  double smd, mela,psig,pbkg, melapsigOut,melapbkgOut;
  double smdSyst1Up, smdSyst1Down, smdSyst2Up, smdSyst2Down, melaTmp,psigTmp,pbkgTmp;
  float psmela,psigps,pbkgps,gravimela;

 string outFileName=outDirName+files[ifile]+"_withSMD_doubleCBonly.root";
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
 outTree->Branch("graviLD",&gravimela,"graviLD/F");
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
 mySMD->SetPathToCards("/afs/cern.ch/user/b/bonato/work/PhysAnalysis/HZZ4L/spin/SuperMELA/CMSSW_5_3_3_patch3/src/HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/SM_inputs_"+str_sqrts+"/");
 // mySMD->RecalculateMELA(true);
 mySMD->SetVerbosity(false);
 mySMD->init();



 double meanCB_err=mySMD->GetSigShapeSystematic("meanCB");
 double sigmaCB_err=mySMD->GetSigShapeSystematic("sigmaCB");
 cout<<"Signal shape syst factors are MeanSystErr="<<meanCB_err<<"  SigmaSystErr="<<sigmaCB_err<<endl;
 cout<<"Looping on test TTree"<<endl;
 int nDiff=0;
 int maxEntries=sigTree->GetEntries();//1000;//sigTree->GetEntries()
 TH1F *hSMD=new TH1F("hsmd"," SuperMELA",200,0.0,1.0);
 mySMD->SetVerbosity(false);

 for(int i=0;i<maxEntries;i++){
   sigTree->GetEntry(i);
   if(i%4000==0)cout<<"Entry #"<<i<<endl;
   if(mzz>180.0||mzz<100.0)continue;

   psmela=oldPSD;
   gravimela=oldGravD;
   
   ///////////////////
    //for regular way, no recalculation of MELA
   //  melapsigOut=melapsig;///melapsig is wrong in stat trees of 17/10/2012
    melapbkgOut=melapbkg;
    melapsigOut=(oldD*melapbkgOut)/(1-oldD);
    float melaNew=melapsigOut/(melapsigOut+melapbkgOut);
    if(melaNew!=oldD)cout<<"ERROR calc PSig: "<<melaNew<<"  vs "<<oldD<<endl;

    mySMD->computeKD(mzz, melapsigOut, melapbkgOut, smd,psig,pbkg,mela);
    if(fabs(mela-double(melaNew))>1.0e-6) cout<<"ERROR calc MELA in SMD: "<<melaNew<<"  vs "<<mela<<"  Diff="<<mela-double(melaNew)<<endl;
    mela=oldD;
   double mzzTmpSig=0.0, mzzTmpBkg=double(mzz);
   double melaTmp2;

   if(isSignal){//signal sample at MH=125 GeV
     mzzTmpSig=double( mzz*(1.0+meanCB_err) );
     if(mzzTmpSig>180.0 || mzzTmpSig<100){
       //cout<<"Entry #"<< i<<"   ATTENTION: changing scale up has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
       mzzTmpSig=mzz;
     }
     mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst1Up,psigTmp,pbkgTmp, melaTmp2);
 
     mzzTmpSig=double( mzz*(1.0-meanCB_err) );
     if(mzzTmpSig>180.0 || mzzTmpSig<100) mzzTmpSig=mzz;
     mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst1Down,psigTmp,pbkgTmp, melaTmp2);
     
     mzzTmpSig= myR->Gaus(1.0,sigmaCB_err);
     if(mzzTmpSig>180.0 || mzzTmpSig<100) mzzTmpSig=mzz;
     mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst2Up,psigTmp,pbkgTmp, melaTmp2);
     smdSyst2Down=smdSyst2Up;
   }//end if isSignal
   else{
     smdSyst1Up=smd;
     smdSyst1Down=smd;
     smdSyst2Up=smd;
     smdSyst2Down=smd;
   }

   /*************
    ////////////////////
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
     //cout<<"Entry #"<< i<<"   ATTENTION: changing scale up has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
     mzzTmpSig=mzz;
   }
   mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst1Up, melaTmp2,psigTmp,pbkgTmp); 
   if( smdSyst1Up<0.0)cout<<"Entry #"<< i<<"   New SMD with scale Up is "<<smdSyst1Up<<endl;
   
   mzzTmpSig=double( mzz*(1.0-meanCB_err) );
   if(mzzTmpSig>180.0 || mzzTmpSig<100){
     // cout<<"Entry #"<< i<<"   ATTENTION: changing scale down has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
     mzzTmpSig=mzz;
   }
   mySMD->computeKD(mzzTmpSig, melapsigOut, melapbkgOut, smdSyst1Down, melaTmp2,psigTmp,pbkgTmp);
   mzzTmpSig= myR->Gaus(1.0,sigmaCB_err);
   if(mzzTmpSig>180.0 || mzzTmpSig<100){
     //  cout<<"Entry #"<< i<<"   ATTENTION: smearing scale has moved mzz outside [100,180]: "<<mzzTmpSig<<endl;
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
   *******/

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

 cout<<"Writing to "<<fout->GetName()<<"   "<<endl;
 fout->cd();
 outTree->Write();
 // cout<<"and deleting file"<<endl;
 delete fout;

 delete  mySMD;
 delete  mypsLD;


 /*
 cout<<"canvas"<<endl;
 TCanvas *c2=new TCanvas("cc2","CC2",1000,1000);
 c2->cd();
 cout<<"canvas 2"<<endl;
 hSMD->SetXTitle("SuperMELA");
 hSMD->SetYTitle("# events");
 hSMD->Draw();
 cout<<"histo drawn"<<endl;
 c2->SaveAs(("can_"+files[ifile]+"_SuperMELA.root").c_str());
 cout<<"saved canbvas"<<endl;
 delete c2;
 cout<<"deleting last pointers"<<endl;

 delete hSMD;
 */

   }//end loop on ifile (samples to process for this channel)
  }//end loop on ich (channels)

}//end main
