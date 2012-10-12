#include <boost/algorithm/string.hpp>
#include <sstream>
#include "../interface/SuperMELA.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/include/HiggsCSandWidth.cc"
using namespace RooFit;

SuperMELA::SuperMELA(double mH,string channel,int LHCsqrts){

  verbose_=true;
  sqrts_=LHCsqrts;
  mHVal_=mH;
  mH_rrv_=new RooRealVar("mH","mH",mHVal_,0.0,sqrts_);
  strChan_=channel;
 
  mela_psig_=0.0;
  mela_pbkg_=0.0;

  recalculateMELA_=false;
  melaProd_=0;

  pathToCards_="../../../HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/SM_inputs_8TeV/";
  //  init();

}

SuperMELA::~SuperMELA(){

  delete sig_CB_;
  delete qqZZ_pdf_;
  if(melaProd_!=0)delete melaProd_;

}

void SuperMELA::SetDecayChannel(string myChan){
  strChan_=myChan;
  bool newChanOK=checkChannel();
  if(verbose_)std::cout<<"Setting decay channel of SuperMELA to "<<strChan_.c_str()<<" , re-initializing..."<<std::endl;
  init();
  if(verbose_&&newChanOK)std::cout<<"Decay channel set successfully to "<<strChan_.c_str()<<std::endl;
}//end SetDecayChannel

void SuperMELA::computeKD(double m4l,double PSigMelaIn,double PBkgMelaIn,double &superMELA,double &MELA,double &Psig,double &Pbkg){
  
  mela_psig_=PSigMelaIn;
  mela_pbkg_=PBkgMelaIn;
  float melaTmp=mela_psig_ / (mela_psig_+mela_pbkg_);

  if(recalculateMELA_){
    double melaValIn= MELA;

    float costheta1, costheta2, phi, costhetastar, phistar1;
    bool withPt = false;
    bool withY = false;
    
    //    melaProd_.computeKD(pL11, id11, pL12, id12, pL21, id21, pL22, id22,
    melaProd_->computeKD(Z1_lept1_,Z1_lept1Id_,Z1_lept2_,Z1_lept2Id_, Z2_lept1_,Z2_lept1Id_,Z2_lept2_, Z2_lept2Id_,
			 costhetastar, costheta1, costheta2, phi, phistar1, melaTmp, mela_psig_, mela_pbkg_,
			withPt, withY);
   
    if(verbose_)std::cout<<"MELA recalculated. MELA value set from "<<melaValIn<<" to "<<melaTmp<<std::endl;
  }//end if recalculateMELA
  
  //
  MELA=double(melaTmp);
  if(m4l<80 || m4l>180){
    if(verbose_)    std::cout<<"WARNING from void SuperMELA::computeKD ! m4l outside range [80, 180]: "<<m4l<<" . Setting SuperMELA to dummy values."<<std::endl;
    Psig =0.0;
    Pbkg =0.0;
    superMELA=-1.0;
    return;
  }

  pair<double,double> supP =  superMelaLikelihoodDiscriminant ( m4l,mela_psig_, mela_pbkg_);
  //


  Psig =supP.first;
  Pbkg =supP.second;
  superMELA=Psig/(Psig+Pbkg);
  if(verbose_)std::cout<<"SuperMELA="<<superMELA<<"   PSig="<<Psig<<"   PBkg="<<Pbkg<<std::endl;

}//end computeKD


void SuperMELA::computeKD(double m4l,bool use4vectors,double &superMELA,double &MELA,double &Psig,double &Pbkg){
  try{
    if(!recalculateMELA_){
      throw 101;
    }

    if(use4vectors&&(Z1_lept1_.P()<=0.0 ||Z1_lept2_.P()<=0.0||Z2_lept1_.P()<=0.0||Z2_lept2_.P()<=0.0) ){
      throw 102;
    }
    if(!use4vectors && (m1_<=0.0||m2_<=0.0||fabs(hs_)>1.0)){
      throw 103;
    }
  }
  catch(int e){

    if(e==101) std::cerr<<"Exception in void SuperMELA::computeKD(double m4l,bool use4vectors,double &superMELA,double &MELA,double &Psig,double &Pbkg): you cannot call this method if you did not turn on the flag for recalculating MELA on the fly."<<std::endl;
    if(e==102)std::cerr<<"Exception in void SuperMELA::computeKD(double m4l,bool use4vectors,double &superMELA,double &MELA,double &Psig,double &Pbkg): wrong value of input 4-momenta. Magnitudes of 3-momenta are "<<Z1_lept1_.P()<<"  "<<Z1_lept2_.P()<< "  "<<Z2_lept1_.P()<<"   "<<Z2_lept2_.P()<<std::endl;
    if(e==103)std::cerr<<"Exception in void SuperMELA::computeKD(double m4l,bool use4vectors,double &superMELA,double &MELA,double &Psig,double &Pbkg): wrong value of input hel angles. M1 "<<m1_<<"  M2="<<m2_<< "  CosThetaStar="<<hs_<<"  CosTheta1="<<h1_<<std::endl;

  }


  float melaTmp=-1.0;
  
  float costheta1, costheta2, phi, costhetastar, phistar1;
  bool withPt = false;
  bool withY = false;
  
  if(use4vectors){
    melaProd_->computeKD(Z1_lept1_,Z1_lept1Id_,Z1_lept2_,Z1_lept2Id_, Z2_lept1_,Z2_lept1Id_,Z2_lept2_, Z2_lept2Id_,
			 costhetastar, costheta1, costheta2, phi, phistar1, melaTmp, mela_psig_, mela_pbkg_,
			 withPt, withY);
  }
  else{
    float dumpt=0.0,dumy=0.0;
    melaProd_->computeKD(m4l,m1_,m2_, hs_, h1_, h2_, phi_, phistar1_, melaTmp, mela_psig_, mela_pbkg_,
			 withPt,dumpt, withY,dumy);
  }
  MELA=double(melaTmp);
  if(verbose_)std::cout<<"MELA recalculated. MELA value set to "<<MELA<<std::endl;


  //
  if(m4l<80 || m4l>180){
    if(verbose_)    std::cout<<"WARNING from void SuperMELA::computeKD ! m4l outside range [80, 180]: "<<m4l<<" . Setting SuperMELA to dummy values."<<std::endl;
    Psig =0.0;
    Pbkg =0.0;
    superMELA=-1.0;
    return;
  }

  pair<double,double> supP =  superMelaLikelihoodDiscriminant ( m4l,mela_psig_, mela_pbkg_);
  //

  Psig =supP.first;
  Pbkg =supP.second;
  superMELA=Psig/(Psig+Pbkg);
  if(verbose_)std::cout<<"SuperMELA="<<superMELA<<"   PSig="<<Psig<<"   PBkg="<<Pbkg<<std::endl;
}//end computeKD


void SuperMELA::init(){
  if(verbose_)std::cout<<"INITIALIZING..."<<std::endl;
  //calculate m4l ranges for the given mH, set range of rrv
  calc_mZZ_range(mHVal_,lowMH_,highMH_);
  if(verbose_)cout<<"Range width="<<highMH_ - lowMH_<<endl;
  m4l_rrv_=new RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",mHVal_,lowMH_,highMH_);//this one is the variable!
  m4l_rrv_->setBins(2000,"fft") ;
  m4l_rrv_->setRange("shape",lowMH_,highMH_);


  //set parameters for signal m4l shape and calculate normalization
  string str_n_CB ;
  string   str_alpha_CB ;
  string   str_mean_CB ;
  string   str_sigma_CB ;
  //  double   mean_BW_d = mHVal_;
  if(verbose_)std::cout<<"Reading signal shape formulas"<<std::endl;
  readSigParsFromFile( str_mean_CB,str_sigma_CB ,str_n_CB ,str_alpha_CB);
  if(verbose_){
    std::cout<<"Read from input card the following formulas: "<<std::endl;
    std::cout<<"Mean RooFormula (string): "<<str_mean_CB.c_str()<<std::endl;
    std::cout<<"Sigma RooFormula (string): "<<str_sigma_CB.c_str()<<std::endl;
  }
  RooRealVar  dummyOne("one","one",1.0);
  dummyOne.setConstant(true);
  char rrvName[96];
  sprintf(rrvName,"CMS_zz4l_n_sig_%s_%d",strChan_.c_str(),int(sqrts_));
  n_CB_=new RooFormulaVar(rrvName,str_n_CB.c_str() ,RooArgList(*mH_rrv_));
  sprintf(rrvName,"CMS_zz4l_alpha_sig_%s_%d",strChan_.c_str(),int(sqrts_));
  alpha_CB_=new RooFormulaVar(rrvName,str_alpha_CB.c_str() ,RooArgList(dummyOne));

  RooRealVar corr_mean_sig("CMS_zz4l_mean_sig_corrMH","CMS_zz4l_mean_sig_corrMH",0.0,-10.0,10.0);
  RooRealVar corr_sigma_sig("CMS_zz4l_sigma_sig_corrMH","CMS_zz4l_sigma_sig_corrMH",0.0,-10.0,10.0);
  mean_CB_=new RooFormulaVar("CMS_zz4l_mean_m_sig",("("+str_mean_CB+")+@0*@1").c_str(),RooArgList(*mH_rrv_,corr_mean_sig));//this is normalized by  mHVal_
  meanTOT_CB_=new RooFormulaVar("CMS_zz4l_mean_sig","(@0+@1)",RooArgList(*mH_rrv_,*mean_CB_));

  if(verbose_){    std::cout<<"Signal Mean vals -> Correction: "<<corr_mean_sig.getVal()<<"  Mean: "<<mean_CB_->getVal()<<"  Total: "<<meanTOT_CB_->getVal()<<std::endl;}

  sigma_CB_=new RooFormulaVar("CMS_zz4l_sigma_m_sig",("("+str_sigma_CB+")*(1+@1)").c_str(),RooArgList(*mH_rrv_,corr_sigma_sig ));
   
  if(verbose_){
    std::cout<<"Signal shape parameter values: "<<std::endl;
    std::cout<<"Mean (formula value) = "<<meanTOT_CB_->getVal()<<std::endl;
    std::cout<<"Sigma (formula value) = "<<sigma_CB_->getVal()<<std::endl;
    std::cout<<"n (formula value) = "<<n_CB_->getVal()<<std::endl;
    std::cout<<"alpha (formula value) = "<<alpha_CB_->getVal()<<std::endl;
  }

  /***
  //////dummy for testing:
  mean_dummy_=new RooRealVar("dummy_mean","DUMMY MEAN",124.8,0.0,200.0);
  sigma_dummy_=new RooRealVar("dummy_mean","DUMMY MEAN",1.55,0.0,200.0);
  alpha_dummy_=new RooRealVar("dummy_mean","DUMMY MEAN",1.48,0.0,200.0);
  n_dummy_=new RooRealVar("dummy_mean","DUMMY MEAN",2.49,0.0,200.0);
  mean_dummy_->setConstant(kTRUE);
  sigma_dummy_->setConstant(kTRUE);
  alpha_dummy_->setConstant(kTRUE);
  n_dummy_->setConstant(kTRUE);
 sig_CB_=new RooCBShape("signalCB_ggH","signalCB_ggH",*m4l_rrv_,*mean_dummy_,*sigma_dummy_,*alpha_dummy_,*n_dummy_);
  ***/

  sig_CB_=new RooCBShape("signalCB_ggH","signalCB_ggH",*m4l_rrv_,*meanTOT_CB_,*sigma_CB_,*alpha_CB_,*n_CB_);
 
  if(verbose_)std::cout<<"Value of signal m4l shape is "<<sig_CB_->getVal()<<std::endl;
  norm_sig_CB_=sig_CB_->createIntegral( RooArgSet(*m4l_rrv_), RooFit::Range("shape"))->getVal();
  if(verbose_)std::cout<<"Normalization of signal m4l shape is "<<norm_sig_CB_<<std::endl;
  //set parameters for background m4l shape and calculate normalization
  if(verbose_)std::cout<<"Reading background shape parameters"<<std::endl;
  std::vector<double> v_apars;
  readBkgParsFromFile(v_apars);
  
  if(verbose_){
    cout<<"Size of vector with bkg shape pars is "<<v_apars.size()<<endl;
    std::cout<<"Param [0]="<<v_apars.at(0)<<" [13]="<<v_apars.at(13)<<std::endl;
  }

  a0_qqZZ_=new RooRealVar("CMS_zz4l_a0_qqZZ","CMS_zz4l_a0_qqZZ",v_apars.at(0),0.,200.);
  a1_qqZZ_=new RooRealVar("CMS_zz4l_a1_qqZZ","CMS_zz4l_a1_qqZZ",v_apars.at(1),0.,200.);
  a2_qqZZ_=new RooRealVar("CMS_zz4l_a2_qqZZ","CMS_zz4l_a2_qqZZ",v_apars.at(2),0.,200.);
  a3_qqZZ_=new RooRealVar("CMS_zz4l_a3_qqZZ","CMS_zz4l_a3_qqZZ",v_apars.at(3),0.,1.);
  a4_qqZZ_=new RooRealVar("CMS_zz4l_a4_qqZZ","CMS_zz4l_a4_qqZZ",v_apars.at(4),0.,200.);
  a5_qqZZ_=new RooRealVar("CMS_zz4l_a5_qqZZ","CMS_zz4l_a5_qqZZ",v_apars.at(5),0.,200.);
  a6_qqZZ_=new RooRealVar("CMS_zz4l_a6_qqZZ","CMS_zz4l_a6_qqZZ",v_apars.at(6),0.,100.);
  a7_qqZZ_=new RooRealVar("CMS_zz4l_a7_qqZZ","CMS_zz4l_a7_qqZZ",v_apars.at(7),0.,1.);
  a8_qqZZ_=new RooRealVar("CMS_zz4l_a8_qqZZ","CMS_zz4l_a8_qqZZ",v_apars.at(8),0.,200.);
  a9_qqZZ_=new RooRealVar("CMS_zz4l_a9_qqZZ","CMS_zz4l_a9_qqZZ",v_apars.at(9),0.,1.);
  a10_qqZZ_=new RooRealVar("CMS_zz4l_a10_qqZZ","CMS_zz4l_a10_qqZZ",v_apars.at(10),0.,200.);
  a11_qqZZ_=new RooRealVar("CMS_zz4l_a11_qqZZ","CMS_zz4l_a11_qqZZ",v_apars.at(11),-100.,100.);
  a12_qqZZ_=new RooRealVar("CMS_zz4l_a12_qqZZ","CMS_zz4l_a12_qqZZ",v_apars.at(12),0.,10000.);
  a13_qqZZ_=new RooRealVar("CMS_zz4l_a13_qqZZ","CMS_zz4l_a13_qqZZ",v_apars.at(13),0.,1.);
  a0_qqZZ_->setConstant(kTRUE);
  a1_qqZZ_->setConstant(kTRUE);
  a2_qqZZ_->setConstant(kTRUE);
  a3_qqZZ_->setConstant(kTRUE);
  a4_qqZZ_->setConstant(kTRUE);
  a5_qqZZ_->setConstant(kTRUE);
  a6_qqZZ_->setConstant(kTRUE);
  a7_qqZZ_->setConstant(kTRUE);
  a8_qqZZ_->setConstant(kTRUE);
  a9_qqZZ_->setConstant(kTRUE);
  a10_qqZZ_->setConstant(kTRUE);
  a11_qqZZ_->setConstant(kTRUE);
  a12_qqZZ_->setConstant(kTRUE);
  a13_qqZZ_->setConstant(kTRUE);


  qqZZ_pdf_ = new RooqqZZPdf_v2("bkg_qqzz","bkg_qqzz",*m4l_rrv_,*a0_qqZZ_, *a1_qqZZ_, *a2_qqZZ_, *a3_qqZZ_,             
				*a4_qqZZ_, *a5_qqZZ_, *a6_qqZZ_, *a7_qqZZ_,
				*a8_qqZZ_, *a9_qqZZ_, *a10_qqZZ_,*a11_qqZZ_,
				*a12_qqZZ_,   *a13_qqZZ_ );

  norm_bkg_qqZZ_=qqZZ_pdf_->createIntegral( RooArgSet(*m4l_rrv_), RooFit::Range("shape"))->getVal();

  if(recalculateMELA_){//instantiate the MELA producer
    if(verbose_)std::cout<<"Creating Mela producer."<<std::endl;
    melaProd_=new Mela(false,sqrts_);//use analytical bkgd parametrization

  }

}//end init

void  SuperMELA::readSigParsFromFile(string &str_mean_CB,string &str_sigma_CB ,string &str_n_CB ,string &str_alpha_CB){

  bool meanOK=false,sigmaOK=false,nOK=false,alphaOK=false;
  //open text file
  string fCardName=pathToCards_+"inputs_"+strChan_+".txt";
  if(verbose_)std::cout<<"Parsing input card "<<fCardName.c_str()<<std::endl;
  ifstream card(fCardName.c_str(),ios::in);
  string line;
  while(card.good()){
    getline(card, line);
    std::vector<string> fields;
    split( fields, line, boost::is_any_of( " " ), boost::token_compress_on );
    if(fields[0]!="signalShape")continue;
    //ok, we found somethign interesting
    if(fields.size()!=3){
      std::cout<<"Error in SuperMELA::readSigParsFromFile! Incorrect format of line "<<line.c_str()<<std::endl; 
      break;
    }
    if(fields[1]=="n_CB"){str_n_CB=fields[2];nOK=true;}
    if(fields[1]=="alpha_CB"){str_alpha_CB=fields[2];alphaOK=true;}
    if(fields[1]=="mean_CB"){str_mean_CB=fields[2];meanOK=true;}
    if(fields[1]=="sigma_CB"){str_sigma_CB=fields[2];sigmaOK=true;}
    
    if(meanOK&&sigmaOK&&alphaOK&&nOK)break;
  }//end while loop on lines

  try{
    if(!(meanOK&&sigmaOK&&alphaOK&&nOK)){
      throw 20;
    }
  }
  catch(int e){
    std::cout<<"Exception "<<e <<" in SuperMELA::readSigParsFromFile! Not all signal shape formulas were read "<<meanOK<<" "<<sigmaOK<<"  "<<alphaOK<<"  "<<nOK<<std::endl;
  }

  card.close();

}//end readSigParsFromFile


void SuperMELA::readBkgParsFromFile(std::vector<double> &apars ){

  const int nPars=14;
  int nFound=0;
  apars.resize(14);
  string fCardName=pathToCards_+"inputs_"+strChan_+".txt";
  if(verbose_)std::cout<<"Parsing input card "<<fCardName.c_str()<<std::endl;
  ifstream card(fCardName.c_str(),ios::in);
  string line;
 
  while(card.good()){
    getline(card, line);
    std::vector<string> fields;
    split( fields, line, boost::is_any_of( " " ), boost::token_compress_on );
    if(fields[0]!="qqZZshape")continue;
    //  cout<<"Bkg Line selected: "<<line.c_str()<<std::endl;
    if(fields.size()!=3){
      std::cout<<"Error in SuperMELA::readSigParsFromFile! Incorrect format of line "<<line.c_str()<<std::endl; 
      break;
    }

    stringstream ssip;
    ssip<<nFound;
    if(fields[1]=="a"+ssip.str()+"_bkgd"){
      //    cout<<"ip="<<nFound<<" Field[1] matched to "<<("a"+ssip.str()+"_bkgd").c_str()<<endl;
      // cout<<"converting to int "<<fields[2].c_str()<<endl;
      apars[nFound]=atof(fields[2].c_str());
      nFound++;
    }
  }//end loop on lines



  try{
    if(nFound!=nPars){
      throw 30;
    }
  }
  catch (int e){
    if(e==30)std::cerr<<"Exception from void SuperMELA::readBkgParsFromFile(std::vector<double> apars ). Mismatched number of params of qqZZ shape read from card file "<<fCardName.c_str()<<" ---> "<<nFound<<" (it should be "<<nPars<<std::endl;
  }
  card.close();
}//end readBkgParsFromFile

void SuperMELA::calc_mZZ_range(const double mHVal,double &low_M,double &high_M){
  HiggsCSandWidth *myCSW = new HiggsCSandWidth();
  double widthHVal =  myCSW->HiggsWidth(1,mHVal);
  if(verbose_)cout<<"Width="<<widthHVal<<endl;
  double windowVal = max( widthHVal, 1. );
  double lowside = 100.;
  if (mHVal >= 275){ lowside = 180.; }
  else { lowside = 100.; }
  low_M = max( (mHVal - 20.*windowVal), lowside) ;
  high_M = min( (mHVal + 15.*windowVal), 999.5) ;

  if(mHVal==125){low_M=105;high_M=140.0;}
  if(mHVal==126){low_M=106;high_M=141.0;}

  delete myCSW;

}//end calc_mZZ_range


std::pair<double,double>  SuperMELA::superMelaLikelihoodDiscriminant (double m4l,double melaPsig,double melaPbkg){

  m4l_rrv_->setVal(m4l);
  if(verbose_)std::cout<<"In SuperMELA::superMelaLikelihoodDiscriminant,  m4l="<<m4l<<"  MELA-psig="<<melaPsig<<"  MELA-Pbkg="<<melaPbkg<<std::endl;
  //calculate value of signal m4l pdf (normalize pdf to 1) 
  double m4lPsig=sig_CB_->getVal() / norm_sig_CB_;
  if(verbose_)std::cout<<"  m4lPsig="<<m4lPsig<<std::flush;
  //calculate value of background m4l pdf  (normalize pdf to 1) 
  double m4lPbkg=qqZZ_pdf_->getVal() / norm_bkg_qqZZ_;
  if(verbose_)std::cout<<"  m4lPbkg="<<m4lPbkg<<std::endl;

  //the angular probs given back by the Mela producer are already normalized to 1
  double Psig=melaPsig*m4lPsig;
  double Pbkg=melaPbkg*m4lPbkg;

  return make_pair(Psig,Pbkg);
}//end superMelaLikelihoodDiscriminant

bool SuperMELA::checkChannel(){
  try{
    if(strChan_!="4mu" && strChan_!="4e" && strChan_!="2e2mu"  && strChan_!="2mu2e" ){
      throw 10;
    }
  }
  catch (int e){
    std::cerr<<"Exception "<<e<<" from SuperMELA::SetDecayChannel(string myChan). Unrecognized string for decay channel: "<<strChan_.c_str()<<std::endl;
    return false;
  }

  if(strChan_=="4mu")ch_=0;
  else if(strChan_=="4e")ch_=1;
  else ch_=2;

  if(strChan_=="2mu2e")strChan_="2e2mu";
  return true;

}
