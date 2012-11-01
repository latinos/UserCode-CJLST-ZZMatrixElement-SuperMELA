void reduce_data_obs(string fileIn){

  TFile *fIn=new TFile(fileIn.c_str(),"UPDATE");
  TTree *tIn=(TTree*)fIn->Get("data_obs");
  double smd,psmela,m4l;
  float gravimela;
  tIn->SetBranchAddress("CMS_zz4l_mass",&m4l);
  tIn->SetBranchAddress("CMS_zz4l_smd",&smd);
  tIn->SetBranchAddress("CMS_zz4l_KD",&psmela);
  tIn->SetBranchAddress("graviLD",&gravimela);

  double gravimelaD;
  TTree *tOut =  new TTree("data_obs_red","data_obs_red");
  tOut->Branch("CMS_zz4l_smd",&smd,"CMS_zz4l_smd/D");
  tOut->Branch("CMS_zz4l_KD",&psmela,"CMS_zz4l_KD/D");
  tOut->Branch("CMS_zz4l_mass",&m4l,"CMS_zz4l_mass/D");
  tOut->Branch("graviLD",&gravimelaD,"graviLD/D");

  for(int i=0;i<tIn->GetEntries();i++){
    tIn->GetEntry(i);
    gravimelaD=gravimela;
    if(m4l>=105&&m4l<=140) tOut->Fill();
  }

  tOut->Write();
  delete fIn;
}
