#include "HiggsCSandWidth.h"
#include "TChain.h"
#include "TH1F.h"

TString inputDir[2]={"root://lxcms02//data/Higgs/rootuplesOut/130720c/JHU","root://lxcms02//data/Higgs/rootuplesOut/130720c/JHU_8TeV"};
TString SMinputDir[2]={"root://lxcms02//data/Higgs/rootuplesOut/130720c/PRODFSR","root://lxcms02//data/Higgs/rootuplesOut/130720c/PRODFSR_8TeV"};

double lumi[2] = {5.015,19.79};

enum model {k0mplus,k0minus,k0hplus,k1minus,k1plus,k2mplus_gg,k2mplus_qqb,k2hplus,k2hminus,k2bplus,kNUM_MODEL};
enum channel {k4e,k4mu,k2e2mu,kNUM_CHAN};

TString chanDir[3]={"4e","4mu","2mu2e"};

//double eventsGen[10]={300000,300000,300000,200000,200000,300000,200000,300000,300000,300000};

//2e2mu branching ratiodouble 
double branchRatio[3][10] =
	{{ .2592, .2382, .2458, .2395, .2466, .2368, .2368, .2453, .2426, .2340 },
	 { .2592, .2382, .2458, .2395, .2466, .2368, .2368, .2453, .2426, .2340 },
	 { .4816, .5236, .5084, .5210, .5068, .5265, .5265, .5094, .5148, .5319 }};

TString sampleName7TeV[10]={"powheg15jhuGenV3H126",
			    "powheg15jhuGenV3PseH126",
			    "powheg15jhuGenV3ScaHH126",
			    "jhuGenV3Vec1MH126",
			    "jhuGenV3Vec1PH126",
			    "jhuGenV3Grav2PMH126",
			    "jhuGenV3qqGravH126",
			    "jhuGenV3Grav2PHH126",
			    "jhuGenV3Grav2MHH126",
			    "jhuGenV3Grav2PBH126"};

TString sampleName8TeV[10]={"powheg15jhuGenV3H126",
			    "powheg15jhuGenV3PseHH126",
			    "powheg15jhuGenV3ScaHH126",
			    "powheg15jhuGenV3Vec1MH126",
			    "powheg15jhuGenV3Vec1PH126",
			    "powheg15jhuGenV3GravH126",
			    "powheg15jhuGenV3qqGravH126",
			    "powheg15jhuGenV3Grav2PHH126",
			    "powheg15jhuGenV3Grav2MHH126",
			    "powheg15jhuGenV3Grav2PBH126"};


TString sampleLabel[10] =
{
	"$0^{+}_{m}$",
	"$0^{-}$",
	"$0^{+}_{h}$",
	"$1^{-}$",
	"$1^{+}$",
	"$2^{+}_{m} (gg)$",
	"$2^{+}_{m} (q\\bar{q})$",
	"$2^{+}_{h}$",
	"$2^{-}_{h}$",
	"$2^{+}_{b}$"
};

double sqrts[2] = {7.,8.};
TString  Estring[2] = { "7TeV", "8TeV"};


double calcAcceptance(model myModel, channel myChan, int e, bool useWeights=true){

  TChain* t = new TChain("SelectedTree");
  if(myModel == k0mplus)
    {
      if(e == 0)
	{
	  t->Add(SMinputDir[e]+"/"+chanDir[myChan]+"/HZZ4lTree_"+sampleName7TeV[myModel]+".root");
	}
      else
	{
	  t->Add(SMinputDir[e]+"/"+chanDir[myChan]+"/HZZ4lTree_"+sampleName8TeV[myModel]+".root");
	}
    }
  else
    {
      if(e == 0)
	{
	  t->Add(inputDir[e]+"/"+chanDir[myChan]+"/HZZ4lTree_"+sampleName7TeV[myModel]+".root");
	}
      else
	{
	  t->Add(inputDir[e]+"/"+chanDir[myChan]+"/HZZ4lTree_"+sampleName8TeV[myModel]+".root");
	}
    }
  if(!useWeights){
    t->Draw("ZZMass>>temp","MC_weight_norm*(ZZMass>106. && ZZMass<141.)");
    return (double)temp->Integral();
  }else{
    t->Draw("ZZMass>>temp","MC_weight*(ZZMass>106. && ZZMass<141.)");
    //cout << "MC_weight*(ZZMass>106. && ZZMass <141.) " << temp->Integral() << endl;;
    return (double)temp->Integral();
  }

}
	    
void tabulateAcceptance(double mass = 126.){

  HiggsCSandWidth* myCSW = new HiggsCSandWidth("../../CombinationPy/CreateDatacards/include/txtFiles");

  double corr_total, corr_norm;
      
  double NexpJP[3][2],NrecoJP[3][2];
  double NexpH[3][2],NrecoH[3][2];
  double NnormJP[3][2];

  cout << " \\multicolumn{9}{|c|}{" << sampleLabel[k0mplus] << "} \\\\ \\hline " << endl;
  cout << "channel & $\\sqrt{s}$ & $f_{i}^{J^P}$ & $\\alpha_{\\rm ideal} (i)$ & $\\epsilon_{\\rm reco} (i)$ & $\\alpha_{\\rm exp} (i)$ & $N^{J^P}_{\\rm exp} (i)$ & $\\alpha_{\\rm norm} (i)$ & $N^{J^P}_{\\rm norm} (i)$\\\\ \\hline " << endl; 
   
    for (int e = 0; e < 2; e++)
      {   
	
	double mass_XsecRatio = myCSW->HiggsCS(1,mass,sqrts[e])/myCSW->HiggsCS(1,126.,sqrts[e]);
	//cout << endl << "XsecRatio: " << mass_XsecRatio << endl;

	for (channel i=k4e; i != kNUM_CHAN; i = channel(i+1)){
   

	
	  NexpH[i][e]=lumi[e]*calcAcceptance(k0mplus,i,e)*mass_XsecRatio;
	  NrecoH[i][e]=NexpH[i][e]*mass_XsecRatio;

	
	  // channel column
	  cout << chanDir[i];
	  
	  //energy column
	  cout << " & " << Estring[e];

	  // BR column
	  cout << " & " << branchRatio[i][k0mplus] ;

	  //\alpha_{ideal}
	  cout << " &  1.0 ";

	  // \epsilon_{reco} column
	  cout << " & " << calcAcceptance(k0mplus,i,e,false)*mass_XsecRatio;

	  //\alpha_{exp}
	  cout << " &  1.0 ";
	  
	  // N_reco column
	  cout << " & " << NrecoH[i][e] ;
	   
	  
	  // \alpha_norm
	  cout << " &  1.0 ";

	  // N_norm column
	  cout << " & " << NrecoH[i][e] ;

	  cout << " \\\\ \\hline ";
	}
	cout << "\\hline ";
	
      }
      
  /////////////////////////////////////
  // tabulate for other JP models
  /////////////////////////////////////
    for(model j=k0minus; j != kNUM_MODEL; j = model(j+1)){
    cout << " \\multicolumn{9}{|c|}{" << sampleLabel[j] << "} \\\\ \\hline " << endl;
    cout << "channel & $\\sqrt{s}$ & $f_{i}^{J^P}$ & $\\alpha_{\\rm ideal} (i)$ & $\\epsilon_{\\rm reco} (i)$ & $\\alpha_{\\rm exp} (i)$ & $N^{J^P}_{\\rm exp} (i)$ & $\\alpha_{\\rm norm} (i)$ & $N^{J^P}_{\\rm norm} (i)$\\\\ \\hline " << endl; 	
    
    // get variable for formulas, reduce redundant calculations
    for(int e = 0; e < 2; e++)
      {
	double mass_XsecRatio = myCSW->HiggsCS(1,mass,sqrts[e])/myCSW->HiggsCS(1,126.,sqrts[e]);
	//cout << endl << "XsecRatio: " << mass_XsecRatio << endl;

	for (channel i=k4e; i != kNUM_CHAN; i = channel(i+1)){
	  
	  NrecoJP[i][e] = lumi[e] *calcAcceptance(j,i,e)*(branchRatio[k2e2mu][k0mplus]/branchRatio[k2e2mu][j])*mass_XsecRatio;
	  NexpJP[i][e] = NexpH[i][e] * (NrecoJP[i][e])/(NrecoH[i][e]);
	}
	
      }
    
    // print calculations for table
    for(int e = 0; e < 2; e++)
      {
	for (channel i=k4e; i != kNUM_CHAN; i = channel(i+1)){
	  
	  // channel column
	  cout << chanDir[i];
	  
	  //energy column
	  cout << " & " << Estring[e];

	  // BR column
	  cout << " & " << branchRatio[i][j] ;
	  
	  // corr_{BR} column
	  if(i==2)
	    cout << " & 1.0 " ;
	  else
	    cout << " & " << branchRatio[k2e2mu][k0mplus]*(1-branchRatio[k2e2mu][j])/(1-branchRatio[k2e2mu][k0mplus])/branchRatio[k2e2mu][j] ;

	  // \epsilon_{reco} column
	  cout << " & " << calcAcceptance(j,i,e,false)*mass_XsecRatio;

	  // corr_{exp}
	  cout << " & " << NexpJP[i][e]/NexpH[i][e] << endl; 
	  
	  // N_reco column
	  cout << " & " << NrecoJP[i][e] ;
	  
	  // corr_{norm}
	  cout << "% N^{Higg} = " << NexpH[0][0]+NexpH[1][0]+NexpH[2][0]+NexpH[0][1]+NexpH[1][1]+NexpH[2][1] << endl;
	  cout << "% N_{exp}^{JP} = " << NexpJP[0][0] << " + " << NexpJP[1][0] << " + " << NexpJP[2][0] << " + " << NexpJP[0][1] << " + " << NexpJP[1][1] << " + " << NexpJP[2][1] << " = " << (NexpJP[0][0]+NexpJP[1][0]+NexpJP[2][0]+NexpJP[0][1]+NexpJP[1][1]+NexpJP[2][1]) << endl;
	  corr_norm = NexpJP[i][e]/NexpH[i][e] * (NexpH[0][0]+NexpH[1][0]+NexpH[2][0]+NexpH[0][1]+NexpH[1][1]+NexpH[2][1])/(NexpJP[0][0]+NexpJP[1][0]+NexpJP[2][0]+NexpJP[0][1]+NexpJP[1][1]+NexpJP[2][1]) ;
	  
	  cout << " & " << corr_norm;

	  //N_{norm}
	  NnormJP[i][e] = corr_norm * NexpH[i][e];
	  cout << " & " << NnormJP[i][e];
	  
	  
	  cout << " \\\\ \\hline ";
	}

	cout << "\\hline "<< endl;;
    }
    cout << "% Sum N_{norm}^{JP} = " << NnormJP[0][0] << " + " << NnormJP[1][0] << " + " << NnormJP[2][0] << " + " << NnormJP[0][1] << " + " << NnormJP[1][1] << " + " << NnormJP[2][1] << " = " << (NnormJP[0][0]+NnormJP[1][0]+NnormJP[2][0]+NnormJP[0][1]+NnormJP[1][1]+NnormJP[2][1]) << endl;
  }
}
