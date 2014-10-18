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


void doBDT_tZq_SM(){


   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

  TString outfileName( "trainingBDT_tZq.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

 
   // global event weights per tree (see below for setting event-wise weights)
   //Double_t signalWeight     = 0.003582;
   //Double_t backgroundWeight = 0.0269;
   
   Double_t signalWeight     = 1;
   Double_t backgroundWeight = 1;
   //proof_mergedN_BjetEqu1_loose_IsoWcan0p03.root
   TFile *input_sig      = TFile::Open( "../TreeReader/outputroot/histofile_tZq.root" );
   TFile *input_wz       = TFile::Open( "../TreeReader/outputroot/histofile_WZ.root" );
   TFile *input_zjets    = TFile::Open( "../TreeReader/outputroot/histofile_Zjets.root" );
   
   
   
   TTree *signal            = (TTree*)input_sig->Get("Ttree_tZq");
   TTree *background     = (TTree*)input_wz->Get("Ttree_WZ");
   //TTree *background_Zjets  = (TTree*)input_zjets->Get("Ttree_Zjets");
   /*
   
    // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,            signalWeight     );
   factory->AddBackgroundTree( background_WZ,      1.);
   //factory->AddBackgroundTree( background_WZ,     0.93);
   //factory->AddBackgroundTree( background_Zjets,  0.07);
     */
   
   
   
   
   
   std::vector<TString > varList;
   varList.push_back("tree_topMass"    );
   varList.push_back("tree_totMass"    );  
   varList.push_back("tree_deltaPhilb" );
   varList.push_back("tree_deltaRlb"   );
   varList.push_back("tree_deltaRTopZ" );
   varList.push_back("tree_asym"       );
   varList.push_back("tree_Zpt"        );
   varList.push_back("tree_ZEta"       );
   varList.push_back("tree_topPt"      );
   varList.push_back("tree_topEta"     );
   varList.push_back("tree_NJets"      );
   varList.push_back("tree_NBJets"      );
   varList.push_back("tree_deltaRZl"  );   
   varList.push_back("tree_deltaPhiZmet");
   varList.push_back("tree_btagDiscri"	);	   
   varList.push_back("tree_leptWPt"	);		  
   varList.push_back("tree_leptWEta"	);		  
   varList.push_back("tree_leadJetPt"	);	    
   varList.push_back("tree_leadJetEta"	);	  
   varList.push_back("tree_deltaPhiZleptW");  
   varList.push_back("tree_met");
   varList.push_back("tree_mTW");  
   
   
   //-----------------------------------------------------------------------
   //--------------------- to include by event weight-----------------------
   //-----------------------------------------------------------------------
   
   // vector has size of number of input variables
   unsigned int varsize = varList.size();
   std::vector<double > vars(  varsize );
   float  treevars[varsize];
   for(unsigned int i=0; i<varsize; i++) treevars[i] = 0;
   float weight_signal=0; 
   float weight_background=0;
   signal->ResetBranchAddresses() ;
   signal->SetBranchAddress( "tree_EvtWeight", &weight_signal );
   for (unsigned int ivar=0; ivar<varsize; ivar++) signal->SetBranchAddress( varList[ivar].Data(), &(treevars[ivar]) );
   for (unsigned int ivar=0; ivar<varsize; ivar++) cout << signal->GetBranchStatus(varList[ivar].Data()) << endl;
   
   for(int i=0; i< signal->GetEntries(); i++){
     cout << __LINE__  << "_________________________" << endl;
     cout << i << endl;
     //if(i == 7978) continue;
     for (unsigned int ivar=0; ivar<varsize; ivar++) {
       //if( signal->GetBranchStatus(varList[ivar].Data()) != 1 ) cout << "bad adress for var " << varList[ivar]<< endl;
        cout << "bad adress for var " << varList[ivar] << "  " << signal->GetBranchStatus(varList[ivar].Data()) << endl;
     }
     
     
     cout << "get entry " << endl;
     signal->GetEntry(i);
     cout << "get entry ok " << endl;
     for (unsigned int ivar=0; ivar<varsize; ivar++) vars[ivar] = double(treevars[ivar]);
     //cout <<"the signal weight " <<  weight_signal << endl;
     //for(unsigned int j=0; j < vars.size(); j++) cout << varList[j] <<  " " << vars[j] << endl;
     
     if (i < signal->GetEntries()/2.0) {cout << "for training " << endl; factory->AddSignalTrainingEvent( vars,  weight_signal); cout << "ok training " << endl;}
     else     {cout << "for testing  " << endl;factory->AddSignalTestEvent    ( vars,  weight_signal); cout << "ok testing " << endl;}
     
   }
  
  /*
  
   cout << __LINE__ << endl;
   // Background (has event weights)
   background->SetBranchAddress( "tree_EvtWeight", &weight_background );
   for (UInt_t ivar=0; ivar<varsize; ivar++) background->SetBranchAddress( varList[ivar].Data() , &(treevars[ivar]) );
   for (UInt_t i=0; i<background->GetEntries(); i++) {
      background->GetEntry(i);
      for (UInt_t ivar=0; ivar<varsize; ivar++) vars[ivar] = treevars[ivar];
      // add training and test events; here: first half is training, second is testing
      // note that the weight can also be event-wise
      if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars,  weight_background);
      else				  factory->AddBackgroundTestEvent    ( vars,  weight_background);
   } */
 
 /*
 
   cout << __LINE__ << endl;
   
   
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
*/

}



void trainingBDT_tZq_SM(){


   doBDT_tZq_SM ();
   
   
   
}



