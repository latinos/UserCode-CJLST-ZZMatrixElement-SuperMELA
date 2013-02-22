void create_dataobs(){
  const float iniLumi=5.261; const string strIniLumi="5.261";
  const float finLumi=12.21; const string strFinLumi="19.63";

  const bool is7TeV=false;
  string str_sqrt = is7TeV ? "7TeV" : "8TeV";
  string str_lumi =  is7TeV ? "5.051" : "19.63";

  string chans[3]={"DoubleOr","DoubleEle","DoubleMu"};
  string pathToOrig= is7TeV ? "root://lxcms02//data/Higgs/rootuplesOut/130203/PRODFSR/data/" : "root://lxcms02//data/Higgs/rootuplesOut/130205/PRODFSR_"+str_sqrt+"/data/";
  for(int ich=0;ich<3;ich++){

    TFile *inFile=TFile::Open((pathToOrig+"HZZ4lTree_"+chans[ich]+".root").c_str());
    TTree *inTree=(TTree*)inFile->Get("SelectedTree");
    float mzz,mela,mzzErr;
    float p0plus,p0plusN, p0minus,p2minus, psigM4l, pbkgM4l,pbkg;
 float p0plusVA,p0hplusVA,p0minusVA,p1plusVA,p1minusVA,p2minimalVA,p2minimalVA_qq,pbkgVA;

    inTree->SetBranchAddress("ZZMass",&mzz);
    // inTree->SetBranchAddress("MC_weight_noxsec",&w);
   
    inTree->SetBranchAddress("p0plus_m4l",&psigM4l);
    inTree->SetBranchAddress("bkg_m4l",&pbkgM4l);
    
    inTree->SetBranchAddress("p0plus_melaNorm",&p0plusN);
    inTree->SetBranchAddress("p0plus_mela",&p0plus);
    inTree->SetBranchAddress("p0minus_mela",&p0minus);
    inTree->SetBranchAddress("p2_mela",&p2minus);
    inTree->SetBranchAddress("bkg_mela",&pbkg);

  inTree->SetBranchAddress("p0plus_VAJHU",&p0plusVA);
  inTree->SetBranchAddress("p0minus_VAJHU",&p0minusVA);
  inTree->SetBranchAddress("p0hplus_VAJHU",&p0hplusVA);
  inTree->SetBranchAddress("p1plus_VAJHU",&p1plusVA);
  inTree->SetBranchAddress("p1_VAJHU",&p1minusVA);
 inTree->SetBranchAddress("p2_VAJHU",&p2minimalVA);
 inTree->SetBranchAddress("p2qqb_VAJHU",&p2minimalVA_qq);
  inTree->SetBranchAddress("bkg_VAMCFMNorm",&pbkgVA);

  // double m4l,mela,m4lErr;
  //  inTree->SetBranchAddress("CMS_zz4l_mass",&m4l);
  //  inTree->SetBranchAddress("melaLD",&mela);
    //   inTree->SetBranchAddress("CMS_zz4l_massErr",&mzzErr);

    TFile *outF=new TFile(("./hzz"+chans[ich]+"_"+str_lumi+".20130205VA.root").c_str(),"RECREATE");
    TTree *outTree=new TTree("data_obs","data_obs");
    double KD, sKD,pseudoKD,graviKD,p0hKD,p1plusKD,p1minusKD,qqgraviKD;
    double m4l,m4lErr;
    outTree->Branch("CMS_zz4l_mass",&m4l,"CMS_zz4l_mass/D");
    outTree->Branch("melaLD",&KD,"melaLD/D");
    outTree->Branch("CMS_zz4l_massErr",&m4lErr,"CMS_zz4l_massErr/D");
    outTree->Branch("CMS_zz4l_KD",&KD,"CMS_zz4l_KD/D");
    outTree->Branch("CMS_zz4l_smd",&sKD,"CMS_zz4l_smd/D");
    outTree->Branch("CMS_zz4l_pseudoKD",&pseudoKD,"CMS_zz4l_pseudoKD/D");
    outTree->Branch("CMS_zz4l_graviKD",&graviKD,"CMS_zz4l_graviKD/D");
    outTree->Branch("CMS_zz4l_p0hplusKD",&p0hKD,"CMS_zz4l_p0hKD/D");
    outTree->Branch("CMS_zz4l_p1plusKD",&p1plusKD,"CMS_zz4l_p1plusKD/D");
    outTree->Branch("CMS_zz4l_p1minusKD",&p1minusKD,"CMS_zz4l_p1minusKD/D");
    outTree->Branch("CMS_zz4l_p2qqKD",&qqgraviKD,"CMS_zz4l_p2qqKD/D");

	
    for(int i=0;i<inTree->GetEntries();i++){
      inTree->GetEntry(i);
      //calculate discriminants from individual probabilities
      KD= p0plusN / (p0plusN+pbkg);
      //  sKD = p0plusN*psigM4l / (p0plusN*psigM4l + pbkg*pbkgM4l)  ;
      // float pseudoKD = p0plus / (p0plus + p0minus);

      // // using VA \\ \\
      sKD = double(p0plusVA)*double(psigM4l) / (double(p0plusVA)*double(psigM4l) + double(pbkgVA)*double(pbkgM4l))  ;
      pseudoKD = double(p0plusVA) / (double(p0plusVA) + double(p0minusVA));
      graviKD =  double(p0plusVA) / ( double(p0plusVA) +  double(p2minimalVA) );
      p0hKD = p0plusVA / (p0plusVA   + p0hplusVA);
      p1plusKD = p0plusVA / (p0plusVA   + p1plusVA);
      p1minusKD = p0plusVA / (p0plusVA   + p1minusVA);
      qqgraviKD =  p0plusVA   / (p0plusVA + p2minimalVA_qq);
      m4l=mzz;
      m4lErr=0.0;//mzzErr;

      if(p0plus<0.0 || p0minus<0.0 || mzz <106.0 || mzz>141.0 ){
	sKD=-1.0;
	//KD=-1.0;
      }

      outTree->Fill();

    }//end loop on entries of inTree

    outF->cd();
    outTree->Write();
    delete outF;
    delete inTree;
    delete inFile;
  }//end loop on chans

}//main
