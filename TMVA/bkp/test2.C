#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "test/TMVAGui.C"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"


void test2(){
  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  TString outfileName( "trainingBDT_tZq.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  
  
  
  TFile *input_sig      = TFile::Open( "../TreeReader/outputroot/histofile_tZq.root" );
  TFile *input_wz       = TFile::Open( "../TreeReader/outputroot/histofile_WZ.root" );
  
  
  TTree *signal            = (TTree*)input_sig->Get("Ttree_tZq");
  TTree *background     = (TTree*)input_wz->Get("Ttree_WZ");
  
  factory->AddSignalTree    ( signal,	   1.);
  factory->AddBackgroundTree( background,  1.);
  
  
  std::vector<TString > varList;
  varList.push_back("tree_cosThetaStar");;
  varList.push_back("tree_topMass");     
  varList.push_back("tree_totMass");     
  varList.push_back("tree_deltaPhilb");  
  varList.push_back("tree_deltaRlb");    
  varList.push_back("tree_deltaRTopZ");  
  varList.push_back("tree_asym");        
  varList.push_back("tree_Zpt");         
  varList.push_back("tree_ZEta");        
  varList.push_back("tree_topPt");       
  varList.push_back("tree_topEta");      
  varList.push_back("tree_NJets");       
  varList.push_back("tree_NBJets");	 
  varList.push_back("tree_deltaRZl");	 
  varList.push_back("tree_deltaPhiZmet");
  varList.push_back("tree_btagDiscri");  
  
  varList.push_back("tree_totPt");	
  varList.push_back("tree_totEta");	
  
  
  varList.push_back("tree_leptWPt");	 
  varList.push_back("tree_leptWEta");	 
  varList.push_back("tree_leadJetPt");   
  varList.push_back("tree_leadJetEta");  
  varList.push_back("tree_deltaRZleptW");
  varList.push_back("tree_deltaPhiZleptW");
  
  
  varList.push_back("tree_met" );
  varList.push_back("tree_mTW" );
  
  
  for(unsigned int i=0; i< varList.size() ; i++) factory->AddVariable( varList[i].Data(),    'F');
  
  factory->SetSignalWeightExpression    ("tree_EvtWeight");
  factory->SetBackgroundWeightExpression("tree_EvtWeight");
   
  
  // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   
   //factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
//   factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );
   factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

 


   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );

  
  
  
}
