//load MELA and SuperMELA libraries before compililng this
//.x loadMELA.C

#include <string>

TString inputDir="root://lxcms02//data/Higgs/rootuplesOut/131012/PRODFSR_8TeV/";
string chan="4mu";
TString discriminantName="ZZpseudoLD";
bool debug=false;

void setInputDir(TString temp){inputDir=temp;}
void setChan(string temp){chan=temp;}

void generateTemplates(TString inputFile="HZZ4lTree_H125.root",TString tag="signal_4mu"){

  Mela myMELA;

  if(chan=="2mu2e") chan="2e2mu";

  SuperMela_proto mySMD(125.0,"./",chan,8);

  SuperMela_proto mySMD_sigmaUp(125.0,"./",chan,8);
  mySMD_sigmaUp.setSigCBsigmaSyst(.2);

  SuperMela_proto mySMD_meanUp(125.0,"./",chan,8);
  mySMD_meanUp.setSigCBmeanSyst(.004);

  SuperMela_proto mySMD_meanDown(125.0,"./",chan,8);
  mySMD_meanDown.setSigCBmeanSyst(-.004);

  TChain* tree = new TChain("SelectedTree");
  tree->Add(inputDir+inputFile);
  
  float mzz,m1,m2,hs,h1,h2,phi,phi1;
  float w, mela;
  float psig_smd, pbkg_smd, smd; 
  float psig_smd_sigmaUp, pbkg_smd_sigmaUp, smd_sigmaUp; 
  float psig_smd_meanUp, pbkg_smd_meanUp, smd_meanUp; 
  float psig_smd_meanDown, pbkg_smd_meanDown, smd_meanDown; 
  float psig, pbkg,kd;

  tree->SetBranchAddress("ZZMass",&mzz);
  tree->SetBranchAddress("Z1Mass",&m1);
  tree->SetBranchAddress("Z2Mass",&m2);
 
  tree->SetBranchAddress("helcosthetaZ2",&h1);
  tree->SetBranchAddress("helcosthetaZ1",&h2);
  tree->SetBranchAddress("costhetastar",&hs);
  tree->SetBranchAddress("helphi",&phi);
  tree->SetBranchAddress("phistarZ1",&phi1);

  tree->SetBranchAddress("MC_weight_noxsec",&w);
  tree->SetBranchAddress(discriminantName,&mela);

  TH2F* SMDtemplate = new TH2F("SMDtemplate","SMDtemplate",30,0,1,30,0,1);
  TH2F* SMDtemplate_sigmaUp = new TH2F("SMDtemplate_sigmaUp","SMDtemplate_sigmaUp",30,0,1,30,0,1);
  TH2F* SMDtemplate_meanUp = new TH2F("SMDtemplate_meanUp","SMDtemplate_meanUp",30,0,1,30,0,1);
  TH2F* SMDtemplate_meanDown = new TH2F("SMDtemplate_meanDown","SMDtemplate_meanDown",30,0,1,30,0,1);

  for( int iEvt=0; iEvt<tree->GetEntries(); iEvt++){

    tree->GetEntry(iEvt);

    if(iEvt%1000==0) cout << iEvt << "/" << tree->GetEntries() << endl;

    if(mzz>105&&mzz<140){
      myMELA.computeKD(mzz,m1,m2,hs,h1,h2,phi,phi1,kd,psig,pbkg);
      mySMD.computeKD(mzz,psig,pbkg,smd,psig_smd,pbkg_smd);
      mySMD_sigmaUp.computeKD(mzz,psig,pbkg,smd_sigmaUp,psig_smd_sigmaUp,pbkg_smd_sigmaUp);
      mySMD_meanUp.computeKD(mzz,psig,pbkg,smd_meanUp,psig_smd_meanUp,pbkg_smd_meanUp);
      mySMD_meanDown.computeKD(mzz,psig,pbkg,smd_meanDown,psig_smd_meanDown,pbkg_smd_meanDown);

      /*
      cout << "smd: " << smd << endl;
      cout << "smd_sigmaUp: " << smd_sigmaUp << endl;
      cout << "smd_meanUp: " << smd_meanUp << endl;
      cout << "smd_meanDown: " << smd_meanDown << endl;
      */

      SMDtemplate->Fill(smd,mela,w);
      SMDtemplate_sigmaUp->Fill(smd_sigmaUp,mela,w);
      SMDtemplate_meanUp->Fill(smd_meanUp,mela,w);
      SMDtemplate_meanDown->Fill(smd_meanDown,mela,w);

    }
  }

  TCanvas* can = new TCanvas("can","can",500,500);
  gStyle->SetOptStat(0);

  SMDtemplate->Draw("COLZ");
  
  TFile *outFile = new TFile("D"+tag+".root","RECREATE");
  SMDtemplate->Write("h_mzzD");
  SMDtemplate_sigmaUp->Write("h_mzzD_sigmaUp");
  SMDtemplate_meanUp->Write("h_mzzD_meanUp");
  SMDtemplate_meanDown->Write("h_mzzD_meanDown");
  outFile->Close();

  delete SMDtemplate;
  delete SMDtemplate_sigmaUp;
  delete SMDtemplate_meanUp;
  delete SMDtemplate_meanDown;

}

void makeAllTemplates(){

  setInputDir("root://lxcms02//data/Higgs/rootuplesOut/131012/PRODFSR_8TeV/4mu/");
  setChan("4mu");
  generateTemplates("HZZ4lTree_H125.root","signal_4mu");
  generateTemplates("HZZ4lTree_ZZTo*.root","background_qqZZ_4mu");

  setInputDir("root://lxcms02//data/Higgs/rootuplesOut/131012/PRODFSR_8TeV/4e/");
  setChan("4e");
  generateTemplates("HZZ4lTree_H125.root","signal_4mu");
  generateTemplates("HZZ4lTree_ZZTo*.root","background_qqZZ_4e");
 
  setInputDir("root://lxcms02//data/Higgs/rootuplesOut/131012/PRODFSR_8TeV/2mu2e/");
  setChan("2e2mu");
  generateTemplates("HZZ4lTree_H125.root","signal_2e2mu");
  generateTemplates("HZZ4lTree_ZZTo*.root","background_qqZZ_2e2mu");

  setInputDir("root://lxcms02//data/Higgs/rootuplesOut/131012/JHU_8TeV/4mu/");
  setChan("4mu");
  generateTemplates("HZZ4lTree_jhuPseH125.root","signal_PS_4mu");
  setInputDir("root://lxcms02//data/Higgs/rootuplesOut/131012/JHU_8TeV/4e/");
  setChan("4e");
  generateTemplates("HZZ4lTree_jhuPseH125.root","signal_PS_4e");
  setInputDir("root://lxcms02//data/Higgs/rootuplesOut/131012/JHU_8TeV/2mu2e/");
  setChan("2e2mu");
  generateTemplates("HZZ4lTree_jhuPseH125.root","signal_PS_2e2mu");
}


