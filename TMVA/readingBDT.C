#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "test/TMVAGui.C"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"



void loopOnSample(){
  
  
  
  
}






void readingBDT(){

 
   // This loads the library
   TMVA::Tools::Instance();
   
   //create the BDT reader
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
   
   TFile *input         = new TFile( "../TreeReader/outputroot/histofile_FCNCzct.root", "read");
//   TFile *target= new TFile( "HistoBDToutput/TMVApp_tZqsm.root","RECREATE" );
   TFile *target= new TFile( "HistoBDToutput/TMVApp_FCNCzct.root","RECREATE" );
   
  std::vector<TString > varList;
//  varList.push_back("tree_cosThetaStar");;
  varList.push_back("tree_topMass");     
//  varList.push_back("tree_totMass");     
  varList.push_back("tree_deltaPhilb");  
//  varList.push_back("tree_deltaRlb");    
//  varList.push_back("tree_deltaRTopZ");  
  varList.push_back("tree_asym");        
  varList.push_back("tree_Zpt");         
  varList.push_back("tree_ZEta");        
  varList.push_back("tree_topPt");       
  varList.push_back("tree_topEta");      
  varList.push_back("tree_NJets");       
  varList.push_back("tree_NBJets");	 
//  varList.push_back("tree_deltaRZl");	 
  varList.push_back("tree_deltaPhiZmet");
  varList.push_back("tree_btagDiscri");  
//  varList.push_back("tree_totPt");	
//  varList.push_back("tree_totEta");	
  varList.push_back("tree_leptWPt");	 
  varList.push_back("tree_leptWEta");	 
  varList.push_back("tree_leadJetPt");   
  varList.push_back("tree_leadJetEta");  
  varList.push_back("tree_deltaRZleptW");
  varList.push_back("tree_deltaPhiZleptW");
  varList.push_back("tree_met" );
  varList.push_back("tree_mTW" );
  
  unsigned int varsize = varList.size();
  float  treevars[varsize];
  for(unsigned int i=0; i<varsize; i++) treevars[i] = 0;
  
  for(unsigned int i=0; i< varList.size() ; i++) reader->AddVariable( varList[i].Data(),  &(treevars[i])  );
  
  
  reader->BookMVA( "BDT", "weights/BDT_trainning_FCNCzct_BDT.weights.xml" );
//  reader->BookMVA( "BDT", "weights/BDT_trainning_tzq_BDT.weights.xml" );
  
  
  
  
   //-----------------------------------------------------
   //for tZq signal
   //-----------------------------------------------------
  input->cd();
  TTree* theTree = (TTree*)input->Get("Ttree_FCNCzct");
//  TTree* theTree = (TTree*)input->Get("Ttree_tZq");
  for (unsigned int ivar=0; ivar<varsize; ivar++) theTree->SetBranchAddress( varList[ivar].Data(), &(treevars[ivar]) );
  
  for(int i=0; i< theTree->GetEntries(); i++){
    theTree->GetEntry(i);
    cout << reader->EvaluateMVA( "BDT") << endl;
    
    
  }
  
  

}


