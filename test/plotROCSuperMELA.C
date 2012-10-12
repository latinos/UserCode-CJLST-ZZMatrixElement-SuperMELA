TGraph* plotROC(TH1F *hSig,TH1F *hBkg, string label){

 //ROC curve
   float binw=hBkg->GetBinWidth(1);
   float tot_exp_bkg=hBkg->Integral();
   float tot_exp_sig=hSig->Integral();
   float exp_sig=0.0,exp_bkg=0.0;// tot_exp_sig, exp_bkg=tot_exp_bkg;
   TH1F *heff_sig=(TH1F*)hSig->Clone();
   heff_sig->Reset();
   heff_sig->SetName("hEFF_sig");
   heff_sig->SetTitle("EFF curve for SM Higg");
   TH1F *heff_bkg=(TH1F*)hBkg->Clone();
   heff_bkg->Reset();
   heff_bkg->SetName("hEFF_bkg");
   heff_bkg->SetTitle("EFF curve for SM qqZZ");
   int nbtmp=hBkg->GetNbinsX();

   const int nb=nbtmp;
   float roc_sig[nb],roc_bkg[nb];

   const float binwidth=hSig->GetBinWidth(1);
   bool flag1=false;
   for(int i=nb;i>=1;i--){
     exp_sig=hSig->Integral(i,nb);
     exp_bkg=hBkg->Integral(i,nb);

     heff_sig->SetBinContent(i,exp_sig/tot_exp_sig);
     heff_bkg->SetBinContent(i,exp_bkg/tot_exp_bkg);

     if(exp_sig<0.1&&!flag1){flag1=true;cout<<"Sig eff < 0.1 at i="<<i <<" (SMD>="<<float(i*binwidth)<<endl;}
     roc_sig[nb-i]=exp_sig/tot_exp_sig;
     roc_bkg[nb-i]=exp_bkg/tot_exp_bkg;
   }

   TGraph *grROC=new TGraph(nb+1,roc_sig,roc_bkg);   
   grROC->SetName(("grROC_"+label).c_str());
   grROC->SetMarkerStyle(20);
   grROC->SetMarkerColor(kBlue);
   TCanvas *c_SD4=new TCanvas("c4","CANVAS - ROC SuperMELA",800,800);
   c_SD4->cd();
   grROC->GetXaxis()->SetTitle("#varepsilon_{SIG}");
   grROC->GetYaxis()->SetTitle("#varepsilon_{BKG}");
   grROC->Draw("AP");
   gPad->SetGrid();
   c_SD4->SaveAs(("can_SuperMELA_ROC"+label+".root").c_str());

   return grROC;


}



void plotROCSuperMELA(){

  TFile *fSigOLD=new TFile("histos_superMELAOLD_H125To4mu.root");
  TFile *fBkgOLD=new TFile("histos_superMELAOLD_ZZTo4mu.root");
  TFile *fSig=new TFile("histos_superMELA_H125To4mu.root");
  TFile *fBkg=new TFile("histos_superMELA_ZZTo4mu.root");

  TH1F *hSigOLD=(TH1F*)fSigOLD->Get("hsmd");
  hSigOLD->SetName("hsmd_SigOLD");
  TH1F *hBkgOLD=(TH1F*)fBkgOLD->Get("hsmd");
  hBkgOLD->SetName("hsmd_BkgOLD");
  TH1F *hSig=(TH1F*)fSig->Get("hsmd");
  hSig->SetName("hsmd_Sig");
  TH1F *hBkg=(TH1F*)fBkg->Get("hsmd");
  hBkg->SetName("hsmd_Bkg");


  TGraph *grROCOLD=(TGraph*)plotROC(hSigOLD,hBkgOLD,"OLD");
  TGraph *grROCNEW=(TGraph*)plotROC(hSig,hBkg,"NEW");

  TCanvas *cROC=new TCanvas("canROC1","CANVAS ROC SuperMELA",900,900);
  cROC->cd();
  grROCOLD->SetMarkerStyle(20);
  grROCOLD->SetMarkerColor(kBlue);
  grROCNEW->SetMarkerStyle(21);
  grROCNEW->SetMarkerColor(kRed);
  grROCOLD->Draw("AP");
  grROCNEW->Draw("P");
  gPad->SetGrid();
  TLegend *l1=new TLegend(0.3,0.5,0.53,0.65);
  l1->AddEntry(grROCOLD,"OLD SuperMELA","P");
  l1->AddEntry(grROCNEW,"NEW SuperMELA","P");
  l1->SetFillColor(kWhite);
  l1->Draw();
  cROC->SaveAs("can_SuperMELA_ROC_OLDandNEW.root");
  cROC->SaveAs("can_SuperMELA_ROC_OLDandNEW.eps");
}
