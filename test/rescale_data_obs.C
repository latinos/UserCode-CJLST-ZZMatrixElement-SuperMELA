void rescale_data_obs(){

  const float iniLumi=5.261; const string strIniLumi="5.261";
  const float finLumi=12.0; const string strFinLumi="12.0";

  string chans[3]={"2e2mu","4e","4mu"};

  TRandom3 *myr=new TRandom3(36567);

  for(int ich=0;ich<3;ich++){
    TFile *inFile=new TFile(("hzz"+chans[ich]+"_"+strIniLumi+".root").c_str());
    TTree *inTree=(TTree*)inFile->Get("data_obs");
    
    double m4l,mela,m4lErr;
    inTree->SetBranchAddress("CMS_zz4l_mass",&m4l);
    inTree->SetBranchAddress("melaLD",&mela);
    inTree->SetBranchAddress("CMS_zz4l_massErr",&m4lErr);

    TFile *outF=new TFile(("hzz"+chans[ich]+"_"+strFinLumi+".TMP.root").c_str(),"RECREATE");
    TTree *outTree=new TTree("data_obs","data_obs");
    outTree->Branch("CMS_zz4l_mass",&m4l,"CMS_zz4l_mass/D");
    outTree->Branch("melaLD",&mela,"melaLD/D");
    outTree->Branch("CMS_zz4l_massErr",&m4lErr,"CMS_zz4l_massErr/D");
 
    const float scaleF=finLumi / iniLumi;
    float diff = scaleF;
    while(diff>0){

      if(diff>1.0){//just copy
	
	for(int i=0;i<inTree->GetEntries();i++){
	  inTree->GetEntry(i);
	  outTree->Fill();
	}//end loop on entries of inTree

      }
      else{
	for(int i=0;i<inTree->GetEntries();i++){
	  inTree->GetEntry(i);
	double rnum=myr->Rndm();
	if(rnum<diff)outTree->Fill();
	   }
      }
      diff--;      
    }//end while loop

    outF->cd();
    outTree->Write();
    delete outF;
    delete inTree;
    delete inFile;
  }//end loop on chans

}//main
