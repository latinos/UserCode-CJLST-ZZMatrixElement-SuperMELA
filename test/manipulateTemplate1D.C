void manipulateBkgTemplate(string channel){


  TFile *f=new TFile(("Dbackground_qqZZ_"+channel+".root").c_str(),"UPDATE");
  TH1F *h1d=(TH1F*)f->Get("h_superD");
  TH1F *hnew=(TH1F*)h1d->Clone("h_superD_mod");

  float bin1=hnew->GetBinContent(1);
  float bin2=hnew->GetBinContent(2);
  float bin3=hnew->GetBinContent(3);
  float bin4=hnew->GetBinContent(4);
  hnew->SetBinContent(1,bin1*0.7);
  hnew->SetBinContent(2,bin2*1.2);
  hnew->SetBinContent(3,bin3*1.2);
  hnew->SetBinContent(4,bin4*1.15);

 
 //smear the other bins
  TH1F *hnew2=(TH1F*)hnew->Clone("h_superD_mod_rndm");
  TRandom3 *myrandom=new TRandom3(9978);
  for(int i=5;i<=hnew2->GetNbinsX();i++){
    hnew2->SetBinContent(i,hnew2->GetBinContent(i) * myrandom->Gaus(1.0,0.10));
  }

  hnew->Scale(1.0/hnew->Integral());
  hnew2->Scale(1.0/hnew->Integral());

  for(int i=1;i<=hnew->GetNbinsX();i++){
    hnew->SetBinError(i,0.);
  }

  hnew->Write();
  hnew2->Write();
  delete myrandom;
  delete f;
}

void manipulateSigTemplate(string channel){


  TFile *f=new TFile(("Dsignal_"+channel+".root").c_str(),"UPDATE");
  TH1F *h1d=(TH1F*)f->Get("h_superD");
  TH1F *hnew=(TH1F*)h1d->Clone("h_superD_mod");


  int nb=hnew->GetNbinsX();
  float bin1=hnew->GetBinContent(nb);
  float bin2=hnew->GetBinContent(nb-1);
  float bin3=hnew->GetBinContent(nb-2);
  float bin4=hnew->GetBinContent(nb-3);
  hnew->SetBinContent(nb,bin1*0.9);
  hnew->SetBinContent(nb-1,bin2*0.95);
  hnew->SetBinContent(nb-2,bin3*1.05);
  hnew->SetBinContent(nb-3,bin4*1.15);

  //smear the other bins
  TH1F *hnew2=(TH1F*)hnew->Clone("h_superD_mod_rndm");
  TRandom3 *myrandom=new TRandom3(9978);
  for(int i=nb-4;i>=1;i--){
    hnew2->SetBinContent(i,hnew2->GetBinContent(i) * myrandom->Gaus(1.0,0.05));
  }


  hnew->Scale(1.0/hnew->Integral());
  hnew2->Scale(1.0/hnew->Integral());

  for(int i=1;i<=hnew->GetNbinsX();i++){
    hnew->SetBinError(i,0.);
  }

  hnew->Write();
  hnew2->Write();
  delete myrandom;
  delete f;
}


void ratio1DSigTemplate(string channel){


  TFile *f=new TFile(("Dsignal_"+channel+".root").c_str(),"UPDATE");
  TH1F *h1d=(TH1F*)f->Get("h_superD");
  TH1F *hnew=(TH1F*)f->Get("h_superDfromProjX");
  TH1F *hmod=(TH1F*)f->Get("h_superD_mod");
  TH1F *hmod2=(TH1F*)f->Get("h_superD_mod_rndm");

  cout<<"Area histo#1: "<<h1d->Integral()<<"   hist#2:"<<hnew->Integral()<<endl;

  TCanvas *c1=new TCanvas("can1D","CAN - 1D Signal templates",800,800);
  c1->cd();
  h1d->SetMarkerColor(kBlue);
  hnew->SetMarkerColor(kMagenta);
  hmod->SetMarkerColor(kOrange);
  h1d->Draw("P");
  hnew->Draw("Psames");
  hmod->Draw("Psames");
  TLegend *l=new TLegend(0.33,0.5,0.5,0.8);
  l->AddEntry(h1d,"ORIG 1D template","P");
  l->AddEntry(hmod,"MODIFIED 1D template","P");
  l->AddEntry(hmod2,"MODIFIED 1D template + randomized","P");
  l->AddEntry(hnew,"Projected 2D template","P");
  l->SetFillColor(kWhite);
  l->Draw();
  c1->SaveAs(("can_plotSig_2DprojX_1D_"+channel+".root").c_str());
  delete c1;
  cout<<"ratio Sig 1"<<endl;
  TH1F *ratio=(TH1F*)hnew->Clone(("ratio_"+channel).c_str());
  // ratio->Sumw2();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kMagenta);
  ratio->Divide(h1d);
  ratio->SetMaximum(1.5);
  ratio->SetMinimum(0.5);
  cout<<"ratio Sig 2"<<endl;
 TH1F *ratiomod=(TH1F*)hmod->Clone(("ratio_mod_"+channel).c_str());
  // ratio->Sumw2();
  ratiomod->SetMarkerStyle(21);
  ratiomod->SetMarkerColor(kOrange);
  ratiomod->Divide(h1d);
  ratiomod->SetMaximum(1.5);
  ratiomod->SetMinimum(0.5);
  cout<<"ratio Sig 3"<<endl;
 TH1F *ratiomod2=(TH1F*)hmod2->Clone(("ratio_mod_rndm_"+channel).c_str());
  // ratio->Sumw2();
  ratiomod2->SetMarkerStyle(21);
  ratiomod2->SetMarkerColor(kGreen+3);
  ratiomod2->Divide(h1d);
  ratiomod2->SetMaximum(1.5);
  ratiomod2->SetMinimum(0.5);


  for(int i=1;i<=ratio->GetNbinsX();i++){
    ratio->SetBinError(i,0.0);
    ratiomod->SetBinError(i,0.0);
    ratiomod2->SetBinError(i,0.0);
  }

  TCanvas *cr=new TCanvas("canRatio","CAN - RATIO - SIG TEMPLATES",800,800);
  cr->cd();
  ratio->Draw("P");
  ratiomod->Draw("Psames");
  ratiomod2->Draw("Psames");
 l->Draw();
  cr->SaveAs(("can_ratioSig_2DprojX_1D_"+channel+".root").c_str());
  delete cr;
  delete f;
}


void ratio1DBkgTemplate(string channel){

  cout<<"Plot bkg "<<channel.c_str()<<endl;
  TFile *f=new TFile(("Dbackground_qqZZ_"+channel+".root").c_str(),"UPDATE");
  TH1F *h1d=(TH1F*)f->Get("h_superD");
  TH1F *hnew=(TH1F*)f->Get("h_superDfromProjX");
  TH1F *hmod=(TH1F*)f->Get("h_superD_mod");
  TH1F *hmod2=(TH1F*)f->Get("h_superD_mod_rndm");

  cout<<"Area histo#1: "<<h1d->Integral()<<"   hist#2:"<<hnew->Integral()<<endl;

  TCanvas *c1=new TCanvas("can1D","CAN - 1D Signal templates",800,800);
  c1->cd();
  h1d->SetMarkerColor(kBlue);
  hnew->SetMarkerColor(kMagenta);
  hmod->SetMarkerColor(kOrange);
  h1d->Draw("P");
  hnew->Draw("Psames");
  hmod->Draw("Psames");
  TLegend *l=new TLegend(0.33,0.5,0.5,0.8);
  l->AddEntry(h1d,"ORIG 1D template","P");
  l->AddEntry(hmod,"MODIFIED 1D template","P");
  l->AddEntry(hmod2,"MODIFIED 1D template + randomized","P");
  l->AddEntry(hnew,"Projected 2D template","P");
  l->SetFillColor(kWhite);
  l->Draw();
  c1->SaveAs(("can_plotBkg_2DprojX_1D_"+channel+".root").c_str());
  delete c1;

  TH1F *ratio=(TH1F*)hnew->Clone(("ratio_"+channel).c_str());
  // ratio->Sumw2();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kMagenta);
  ratio->Divide(h1d);
  ratio->SetMaximum(1.5);
  ratio->SetMinimum(0.5);

 TH1F *ratiomod=(TH1F*)hmod->Clone(("ratio_mod_"+channel).c_str());
  // ratio->Sumw2();
  ratiomod->SetMarkerStyle(21);
  ratiomod->SetMarkerColor(kOrange);
  ratiomod->Divide(h1d);
  ratiomod->SetMaximum(1.5);
  ratiomod->SetMinimum(0.5);

 TH1F *ratiomod2=(TH1F*)hmod2->Clone(("ratio_mod_rndm_"+channel).c_str());
  // ratio->Sumw2();
  ratiomod2->SetMarkerStyle(21);
  ratiomod2->SetMarkerColor(kGreen+3);
  ratiomod2->Divide(h1d);
  ratiomod2->SetMaximum(1.5);
  ratiomod2->SetMinimum(0.5);


  for(int i=1;i<=ratio->GetNbinsX();i++){
    ratio->SetBinError(i,0.0);
    ratiomod->SetBinError(i,0.0);
    ratiomod2->SetBinError(i,0.0);
  }

  TCanvas *cr=new TCanvas("canRatio","CAN - RATIO - SIG TEMPLATES",800,800);
  cr->cd();
  ratio->Draw("P");
  ratiomod->Draw("Psames");
  ratiomod2->Draw("Psames");
 l->Draw();
  cr->SaveAs(("can_ratioBkg_2DprojX_1D_"+channel+".root").c_str());
  delete cr;
  delete f;
}

void manipulateTemplate1D(){
  manipulateBkgTemplate("4mu");
  manipulateBkgTemplate("4e");
  manipulateBkgTemplate("2e2mu");

  manipulateSigTemplate("4mu");
  manipulateSigTemplate("4e");
  manipulateSigTemplate("2e2mu");

  ratio1DSigTemplate("4mu");
  ratio1DSigTemplate("2e2mu");
  ratio1DSigTemplate("4e");

  ratio1DBkgTemplate("4mu");
  ratio1DBkgTemplate("2e2mu");
  ratio1DBkgTemplate("4e");
}
