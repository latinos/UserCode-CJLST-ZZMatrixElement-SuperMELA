//
// MELA root package loader - see testKD.C for instructions
//
{

  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->Load("libZZMatrixElementMELA.so");
  gSystem->Load("libZZMatrixElementSuperMELA.so");
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/ ");  
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h+");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/SuperMELA/interface/SuperMELA.h+");
  //  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/SuperMELA/src/SuperMela_proto.h+");
  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/PseudoMELA.h+");
}
