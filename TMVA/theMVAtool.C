#define theMVAtool_cxx
#include "theMVAtool.h"

//#include <boost/algorithm/string/replace.hpp>

theMVAtool::theMVAtool(bool doCtrlReg){
  
  //constructor
//  varList.push_back("mTW");
  varList.push_back("wQuark1Pt");
  varList.push_back("wQuark1Eta");
  varList.push_back("wQuark1Phi");
  varList.push_back("wQuark2Pt");
  varList.push_back("wQuark2Eta");
  varList.push_back("wQuark2Phi");
  varList.push_back("wPairMass");
  varList.push_back("wPairPt");
  varList.push_back("wPairEta");
  varList.push_back("wPairPhi");
  varList.push_back("met");
  varList.push_back("nJets");
  varList.push_back("leadJetPt");
  varList.push_back("leadJetPhi");
  varList.push_back("leadJetEta");
  varList.push_back("leadJetbTag");
  varList.push_back("secJetPt");
  varList.push_back("secJetPhi");
  varList.push_back("secJetEta");
  varList.push_back("secJetbTag");
  varList.push_back("thirdJetPt");
  varList.push_back("thirdJetPhi");
  varList.push_back("thirdJetEta");
  varList.push_back("thirdJetbTag");
  varList.push_back("fourthJetPt");
  varList.push_back("fourthJetPhi");
  varList.push_back("fourthJetEta");
  varList.push_back("fourthJetbTag");
  varList.push_back("nBjets");
  varList.push_back("bTagDisc");
  varList.push_back("lep1Pt");
  varList.push_back("lep1Eta");
  varList.push_back("lep1Phi");
  varList.push_back("lep1RelIso");
  varList.push_back("lep1D0");
  varList.push_back("lep2Pt");
  varList.push_back("lep2Eta");
  varList.push_back("lep2Phi");
  varList.push_back("lep2RelIso");
  varList.push_back("lep2D0");
  varList.push_back("lepMass");
  varList.push_back("lepPt");
  varList.push_back("lepEta");
  varList.push_back("lepPhi");
  varList.push_back("zMass");
  varList.push_back("zPt");
  varList.push_back("zEta");
  varList.push_back("zPhi");
  varList.push_back("topMass");
  varList.push_back("topPt");
  varList.push_back("topEta");
  varList.push_back("topPhi");
  varList.push_back("jjdelR");
  varList.push_back("jjdelPhi");
  varList.push_back("wwdelR");
  varList.push_back("wwdelPhi");
  varList.push_back("zLepdelR");
  varList.push_back("zLepdelPhi");
  varList.push_back("zl1Quark1DelR");
  varList.push_back("zl1Quark1DelPhi");
  varList.push_back("zl1Quark2DelR");
  varList.push_back("zl1Quark2DelPhi");
  varList.push_back("zl2Quark1DelR");
  varList.push_back("zl2Quark1DelPhi");
  varList.push_back("zl2Quark2DelR");
  varList.push_back("zl2Quark2DelPhi");
  varList.push_back("zlb1DelR");
  varList.push_back("zlb1DelPhi");
  varList.push_back("zlb2DelR");
  varList.push_back("zlb2DelPhi");
  varList.push_back("lepHt");
  varList.push_back("wQuarkHt");
  varList.push_back("totPt");
  varList.push_back("totEta");
  varList.push_back("totPtVec");
  varList.push_back("totVecM");
//  varList.push_back("Channel");
  varList.push_back("totPt2Jet");
  varList.push_back("wzdelR");
  varList.push_back("wzdelPhi");
  varList.push_back("zQuark1DelR");
  varList.push_back("zQuark1DelPhi");
  varList.push_back("zQuark2DelR");
  varList.push_back("zQuark2DelPhi");
  varList.push_back("zTopDelR");
  varList.push_back("zTopDelPhi");
  varList.push_back("zl1TopDelR");
  varList.push_back("zl1TopDelPhi");
  varList.push_back("zl2TopDelR");
  varList.push_back("zl2TopDelPhi");
  varList.push_back("zjminR");
  varList.push_back("zjminPhi");
  varList.push_back("totHt");
  varList.push_back("jetHt");
  varList.push_back("totHtOverPt");

  if (doCtrlReg){
    regList.push_back("sig_");
    regList.push_back("ctrl_");
  }
  else{
    regList.push_back("");
  }

  for(unsigned int i=0; i< varList.size() ; i++) theVarMap[varList[i]] = i;
  
  samplelist.push_back("DataMu");
  samplelist.push_back("DataEG");
  samplelist.push_back("tZq");
//  samplelist.push_back("THQ"); // Empty tree - background totally cut
  samplelist.push_back("TTZ");
  samplelist.push_back("TTW");
  samplelist.push_back("TT");
//  samplelist.push_back("DYToLL_M10To50"); // Empty tree - background totally cut
  samplelist.push_back("DYToLL_M50");
//  samplelist.push_back("Wjets"); // Empty tree - background totally cut
//  samplelist.push_back("WW"); // Empty tree - background totally cut
//  samplelist.push_back("WZ"); // Empty tree - background totally cut for sig/ctrl, all has 1 event
  samplelist.push_back("ZZ");
//  samplelist.push_back("TsChan"); // Empty tree - background totally cut
  samplelist.push_back("TtChan"); // Empty tree - background totally cut
  samplelist.push_back("TbartChan");
  samplelist.push_back("TtW");
  samplelist.push_back("TbartW");
  
  systlist.push_back("");  
  systlist.push_back("__trig__plus");
  systlist.push_back("__trig__minus");
  systlist.push_back("__jer__plus"); 
  systlist.push_back("__jer__minus");
  systlist.push_back("__jes__plus");
  systlist.push_back("__jes__minus");
  systlist.push_back("__pileup__plus");
  systlist.push_back("__pileup__minus");
  systlist.push_back("__bTag__plus");
  systlist.push_back("__bTag__minus");
  systlist.push_back("__met__plus");
  systlist.push_back("__met__minus");
  systlist.push_back("__pdf__plus");
  systlist.push_back("__pdf__minus");
  systlist.push_back("__ME_PS__plus");
  systlist.push_back("__ME_PS__minus"); 
  
  } 


theMVAtool::theMVAtool(std::vector<TString > thevarlist, std::vector<TString > thesamplelist, std::vector<TString > thesystlist){
  for(unsigned int i=0; i<thevarlist.size(); i++) { varList.push_back(thevarlist[i]);theVarMap[varList[i]] = i;}
  for(unsigned int i=0; i<thesamplelist.size(); i++) { samplelist.push_back(thesamplelist[i]);}
  for(unsigned int i=0; i<thesystlist.size(); i++) { systlist.push_back(thesystlist[i]);}
} 




//---------------------------------------------------------------
 // to perform the training
//---------------------------------------------------------------

void theMVAtool::doTraining(TString channel, TString inDir, bool sigMode){

 //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  TString outfileName( ("trainingBDT_"+channel+"_tZq_.root").Data() );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  //TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_"+channel+"_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
  
  TFile *input_sig        = TFile::Open( inDir+"/histofile_tZq.root" );
  TFile *input_TTZ        = TFile::Open( inDir+"/histofile_TTZ.root" );
  TFile *input_TTW        = TFile::Open( inDir+"/histofile_TTW.root" );
//  TFile *input_THQ        = TFile::Open( inDir+"/histofile_THQ.root" ); // Empty tree - background totally cut

  TFile *input_TT         = TFile::Open( inDir+"/histofile_TT.root" );  
//  TFile *input_WW         = TFile::Open( inDir+"/histofile_WW.root" ); // Empty tree - background totally cut
//  TFile *input_WZ         = TFile::Open( inDir+"/histofile_WZ.root" );
  TFile *input_ZZ         = TFile::Open( inDir+"/histofile_ZZ.root" );

//  TFile *input_TtChan     = TFile::Open( inDir+"/histofile_TtChan.root" ); // Empty tree - background totally cut
  TFile *input_TbartChan  = TFile::Open( inDir+"/histofile_TbartChan.root" );
//  TFile *input_TsChan     = TFile::Open( inDir+"/histofile_TsChan.root" ); // Empty tree - background totally cut
  TFile *input_TtW        = TFile::Open( inDir+"/histofile_TtW.root" );
  TFile *input_TbartW     = TFile::Open( inDir+"/histofile_TbartW.root" );

  TFile *input_DY50	  = TFile::Open( inDir+"/histofile_DYToLL_M50.root" );
//  TFile *input_DY10To50	  = TFile::Open( inDir+"/histofile_DYToLL_M10To50.root" ); // Empty tree - background totally cut
  
  TString treePost = "";
  if (regList.size() > 1) { 
    if ( sigMode ) treePost = "sig_";
    else treePost = "ctrl_";
  }
  TTree *signal              = (TTree*)input_sig->Get("Ttree_"+treePost+"tZq");
  TTree *background_TTZ     = (TTree*)input_TTZ->Get("Ttree_"+treePost+"TTZ");
  TTree *background_TTW     = (TTree*)input_TTW->Get("Ttree_"+treePost+"TTW");
//  TTree *background_THQ     = (TTree*)input_TTW->Get("Ttree_"+treePost+"THQ");

  TTree *background_TT     = (TTree*)input_TT->Get("Ttree_"+treePost+"TT");
//  TTree *background_WW     = (TTree*)input_WW->Get("Ttree_"+treePost+"WW");
//  TTree *background_WZ     = (TTree*)input_WZ->Get("Ttree_"+treePost+"WZ");
  TTree *background_ZZ     = (TTree*)input_ZZ->Get("Ttree_"+treePost+"ZZ");
  
//  TTree *background_TtChan    = (TTree*)input_TtChan->Get("Ttree_"+treePost+"TtChan");
  TTree *background_TbartChan = (TTree*)input_TbartChan->Get("Ttree_"+treePost+"TbartChan");
//  TTree *background_TsChan    = (TTree*)input_TsChan->Get("Ttree_"+treePost+"TsChan");
  TTree *background_TtW       = (TTree*)input_TtW->Get("Ttree_"+treePost+"TtW");  
  TTree *background_TbartW    = (TTree*)input_TbartW->Get("Ttree_"+treePost+"TbartW");  

  TTree *background_DY50         = (TTree*)input_DY50->Get("Ttree_"+treePost+"DYToLL_M50");
//  TTree *background_DY10To50   = (TTree*)input_DY10To50->Get("Ttree_"+treePost+"DYToLL_M10To50");

  factory->AddSignalTree      ( signal, 1. );
  factory->AddBackgroundTree  ( background_TTZ, 1. );
  factory->AddBackgroundTree  ( background_TTW, 1. );

  factory->AddBackgroundTree  ( background_TT, 1. );
//  factory->AddBackgroundTree  ( background_WZ, 1. );
  factory->AddBackgroundTree  ( background_ZZ, 1. );

//  factory->AddBackgroundTree  ( background_TtChan, 1. );
  factory->AddBackgroundTree  ( background_TbartChan, 1. );
//  factory->AddBackgroundTree  ( background_TsChan, 1. );
  factory->AddBackgroundTree  ( background_TtW, 1. );
  factory->AddBackgroundTree  ( background_TbartW, 1. );

  factory->AddBackgroundTree  ( background_DY50, 1. );
//  factory->AddBackgroundTree  ( background_DY10To50, 1. );

  for(unsigned int i=0; i< varList.size() ; i++) factory->AddVariable( varList[i].Data(),    'F');
  
  factory->SetSignalWeightExpression    ("EvtWeight");
  factory->SetBackgroundWeightExpression("EvtWeight");
 
  // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
   
   if(channel !="all"){
     if     (channel == "mumu") mycuts = "Channel == 0";
     else if(channel == "ee"   ) mycuts = "Channel == 1";
     else cout << "WARNING wrong channel name while training " << endl;
   }
   
   
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   //for WZ
   //factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );   
   factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=5:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
   /*for (int nTrees = 25; nTrees < 125; nTrees+=25){
     for (int nCuts = 20; nCuts < 51; nCuts+=10){
       for (int mDepth = 2; mDepth < 5; mDepth++){
	 TString arg;
	 arg.Form("!H:!V:NTrees=%i:nEventsMin=100:MaxDepth=%i:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=%i:PruneMethod=NoPruning",nTrees,mDepth,nCuts);
	 TString name;
	 name.Form("BDT_%i_%i_%i",nTrees,mDepth,nCuts);
	 factory->BookMethod( TMVA::Types::kBDT,name,arg);
       }
     }
     }*/
   //   factory->BookMethod( TMVA::Types::kBDT, "BDT1", "!H:!V:NTrees=150:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=30:PruneMethod=NoPruning" );
   //   factory->BookMethod( TMVA::Types::kBDT, "BDT2", "!H:!V:NTrees=190:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=30:PruneMethod=NoPruning" );
   
 

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


//---------------------------------------------------------------
// to perform the reading
// also fills histograms of the varList used into the BDT
//---------------------------------------------------------------

void theMVAtool::doReading(float bdtcut, TString channel, TString inDir, TString outDir, bool usePseudoData){
  
  
  // This loads the library
  TMVA::Tools::Instance();
  
  //create the BDT reader
  reader = new TMVA::Reader( "!Color:!Silent" );  
  unsigned int varsize = varList.size();
  float  treevars[varsize];
  for(unsigned int i=0; i<varsize; i++) treevars[i] = 0;  
  for(unsigned int i=0; i< varList.size() ; i++) reader->AddVariable( varList[i].Data(),  &(treevars[i])  );
//  reader->AddVariable("mTW", &(treevars[varsize-1]));
   
  reader->BookMVA( "BDT", ("weights/BDT_trainning_"+channel+"_tzq_BDT.weights.xml").Data() ); 
  
  //cout << __LINE__ << "   SAMPLE SIZE " << samplelist.size() << endl;
  for(unsigned int i=0; i< samplelist.size(); i++){
    //cout << samplelist[i] << endl;
    TFile *input         = new TFile( (inDir+"histofile_"+samplelist[i]+".root").Data(), "read");
    TFile *theoutputfile = new TFile(   (outDir+"output_"+channel+"_"+samplelist[i]+".root").Data(), "recreate");
    theOutputFileMap[samplelist[i]] = theoutputfile;

    //---------------------
    //loop over systematics
    //---------------------
    cout <<  "processing sample " <<  samplelist[i]  << endl;
   
   
    for(unsigned int j=0; j<systlist.size(); j++){
      for (unsigned int k=0;k<regList.size(); k++){
	if( (samplelist[i] == "DataMu" || samplelist[i] == "DataEG") && systlist[j] != "") continue;
//	if( (samplelist[i] == "DataMuZenriched" || samplelist[i] == "DataEGZenriched") && !(systlist[j] == "" || systlist[j] == "__zPt__plus" || systlist[j] == "__zPt__minus") ) continue;
	createHisto(regList[k]+samplelist[i]+systlist[j], channel);
	//      cout << samplelist[i]+systlist[j] << channel << endl;
	loopInSample(input, regList[k]+samplelist[i]+systlist[j], treevars, bdtcut, channel);
	writeHisto(samplelist[i], systlist[j], regList[k]);
      }
    }
   
    delete theoutputfile;
  }
  
  
  
  
} 

void theMVAtool::loopInSample(TFile* input, TString sample, float *treevars, float bdtcut, TString channel){
  
  
  unsigned int varsize = varList.size();
  
  /*if (sample.Contains("ZZ")){
    fillHisto(sample, treevars, 0., 0);
    return;
    }*/

  input->cd();
  //  cout << "sample " << sample << endl;
  // cout << "tree name " << "Tree_"+sample  << endl;
  TTree* theTree = (TTree*)input->Get("Ttree_"+sample);
  for (unsigned int ivar=0; ivar<varsize; ivar++) theTree->SetBranchAddress( varList[ivar].Data(), &(treevars[ivar]) );
  float theweight = 0;
  theTree->SetBranchAddress( "EvtWeight", &theweight );
  int theChannel = -1;
  theTree->SetBranchAddress( "Channel", &theChannel );
  float mtwValue = -1;
  theTree->SetBranchAddress( "mTW", &mtwValue );
  int theSelChannel = -1;
  if(     channel == "mumu") theSelChannel = 0;
  else if(channel == "ee" ) theSelChannel = 1;
  
  if(theTree == 0) cout << "no TTree found with name " << "Ttree_"+sample << endl;
  for(int i=0; i< theTree->GetEntries(); i++){
    theTree->GetEntry(i);
    
    double sf_local = 1.;
        
/*    cout << "----------------------" << endl;
    cout << "channel " << channel << endl;
    cout << "theChannel    " << theChannel << endl;
    cout << "theSelChannel " <<  theSelChannel<< endl;*/
    if(channel != "all" && theChannel != theSelChannel ) continue;
    
    //cout << "pass sel " << endl;
/*    if(sample.Contains("WZ") ){
      //cout << "theChannel " << theChannel << endl;
      if(theChannel == 0) sf_local = sf_WZ[0]; 
      if(theChannel == 1) sf_local = sf_WZ[1]; 
      if(theChannel == 2) sf_local = sf_WZ[2]; 
      if(theChannel == 3) sf_local = sf_WZ[3]; 
    } 
    
    if(sample.Contains("Zjets") || sample.Contains("DYToLL_M10To50")){
      if(theChannel == 0) sf_local = sf_DY[0]; 
      if(theChannel == 1) sf_local = sf_DY[1]; 
      if(theChannel == 2) sf_local = sf_DY[2]; 
      if(theChannel == 3) sf_local = sf_DY[3]; 
    } */
    //    if(sample == "TTZ") cout << "weight 1 " << theweight << endl;
    theweight*=sf_local;
    //if(sample == "TTZ") cout << "weight 2 " << theweight << endl;
    
    double mvaValue = reader->EvaluateMVA( "BDT" );
    //cout << reader->EvaluateMVA( "BDT") << endl;
    if(mvaValue > bdtcut) continue;
    //cout << "the bdt " << bdtcut << endl;
    fillHisto(sample, treevars, mvaValue, mtwValue, theweight);
  }  
  
}



void theMVAtool::createHisto(TString sample, TString channel){
  std::vector<TH1F*> histovect;
  for(unsigned int j=0; j<varList.size(); j++){  
    
    int nbins=1;
    double xmin = -1000;
    double xmax = 1000;

    if(varList[j]=="mTW"	     ){     nbins = 10; xmin = 0;   xmax = 200;};	
    if(varList[j]=="wQuark1Pt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="wQuark1Eta"	     ){     nbins = 20; xmin = -5.; xmax = 5.;};
    if(varList[j]=="wQuark1Phi"	     ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="wQuark2Pt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="wQuark2Eta"	     ){     nbins = 20; xmin = -5.; xmax = 5.;};
    if(varList[j]=="wQuark2Phi"	     ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="wPairMass"	     ){     nbins = 10; xmin = 0;   xmax = 200;};
    if(varList[j]=="wPairPt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="wPairEta"	     ){     nbins = 20; xmin = -5.;   xmax = 5.;};
    if(varList[j]=="wPairPhi"	     ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="met"	     ){     nbins = 10; xmin = 0;   xmax = 200;};	   
    if(varList[j]=="nJets"           ){	    nbins = 8;  xmin = 0.5; xmax = 8.5;};      
    if(varList[j]=="leadJetPt"       ){     nbins = 20; xmin = 0;   xmax = 500;};	      
    if(varList[j]=="leadJetPhi"      ){     nbins = 20; xmin = -3.2; xmax = 3.2;};	       
    if(varList[j]=="leadJetEta"      ){     nbins = 10; xmin = -2.5; xmax = 2.5;};	       
    if(varList[j]=="leadJetbTag"     ){     nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="secJetPt"        ){     nbins = 15; xmin = 0;   xmax = 300;};
    if(varList[j]=="secJetPhi"       ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="secJetEta"       ){     nbins = 20; xmin = -5.;   xmax = 5.;};
    if(varList[j]=="secJetbTag"      ){     nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="thirdJetPt"        ){     nbins = 15; xmin = 0;   xmax = 300;};
    if(varList[j]=="thirdJetPhi"       ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="thirdJetEta"       ){     nbins = 20; xmin = -5.;   xmax = 5.;};
    if(varList[j]=="thirdJetbTag"      ){     nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="fourthJetPt"        ){     nbins = 15; xmin = 0;   xmax = 300;};
    if(varList[j]=="fourthJetPhi"       ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="fourthJetEta"       ){     nbins = 20; xmin = -5.;   xmax = 5.;};
    if(varList[j]=="fourthJetbTag"      ){     nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="nBjets"          ){     nbins = 4;  xmin = -0.5;xmax = 3.5;};  
    if(varList[j]=="bTagDisc"        ){     nbins = 10; xmin = 0.6;   xmax = 1;};
    if(varList[j]=="lep1Pt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="lep1Eta"	     ){     nbins = 20; xmin = -5.; xmax = 5.;};
    if(varList[j]=="lep1Phi"	     ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="lep1RelIso"	     ){     nbins = 50; xmin = 0.;   xmax = 0.2;};
    if(varList[j]=="lep1D0"	     ){     nbins = 10; xmin = 0.;   xmax = 0.1;};
    if(varList[j]=="lep2Pt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="lep2Eta"	     ){     nbins = 20; xmin = -5.; xmax = 5.;};
    if(varList[j]=="lep2Phi"	     ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="lep2RelIso"	     ){     nbins = 50; xmin = 0.;   xmax = 0.2;};
    if(varList[j]=="lep2D0"	     ){     nbins = 10; xmin = 0.;   xmax = 0.1;};
    if(varList[j]=="lepMass"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="lepPt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="lepEta"	     ){     nbins = 20; xmin = -5.; xmax = 5.;};
    if(varList[j]=="lepPhi"	     ){     nbins = 20; xmin = -3.2;   xmax = 3.2;};
    if(varList[j]=="zMass"	     ){	    nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="zPt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="zEta"	     ){	    nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="zPhi"	     ){	    nbins = 20; xmin = -3.2;xmax = 3.2;};
    if(varList[j]=="topMass"	     ){	    nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="topPt"	     ){	    nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="topEta"	     ){	    nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="topPhi"	     ){	    nbins = 20; xmin = -3.2;xmax = 3.2;};
    if(varList[j]=="jjdelR"          ){     nbins = 20;xmin = 0; xmax = 6;};
    if(varList[j]=="jjdelPhi"        ){     nbins = 10;xmin = 0; xmax = 3.2;};
    if(varList[j]=="wwdelR"          ){     nbins = 20;xmin = 0; xmax = 6;};
    if(varList[j]=="wwdelPhi"        ){     nbins = 10;xmin = 0; xmax = 3.2;};
    if(varList[j]=="zLepdelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zLepdelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zl1Quark1DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zl1Quark1DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zl1Quark2DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zl1Quark2DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zl2Quark1DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zl2Quark1DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zl2Quark2DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zl2Quark2DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zlb1DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zlb1DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zlb2DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zlb2DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="lepHt"     	     ){     nbins = 15;xmin = 0; xmax = 4;};
    if(varList[j]=="wQuarkHt" 	     ){     nbins = 15;xmin = 0; xmax = 150;};
    if(varList[j]=="totPt"	     ){     nbins = 10; xmin = 0;   xmax = 250;};
    if(varList[j]=="totEta"	     ){	    nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="totPtVec"        ){     nbins = 20;xmin = 0; xmax = 200;};
    if(varList[j]=="totVecM"         ){     nbins = 20;xmin = 0; xmax = 300;};
    if(varList[j]=="totPt2Jet"       ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="wzdelR"          ){     nbins = 20;xmin = 0; xmax = 7;};
    if(varList[j]=="wzdelPhi"        ){     nbins = 20;xmin = 0; xmax = 3.2;};
    if(varList[j]=="zQuark1DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zQuark1DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zQuark2DelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zQuark2DelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zTopDelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zTopDelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zl1TopDelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zl1TopDelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zl2TopDelR"        ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zl2TopDelPhi"      ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zjminR"          ){     nbins = 15;xmin = 0; xmax = 3.5;};
    if(varList[j]=="zjminPhi"        ){     nbins = 15;xmin = 0; xmax = 3.2;};
    if(varList[j]=="totHt"           ){     nbins = 20;xmin = 150; xmax = 1500;};
    if(varList[j]=="jetHt"           ){     nbins = 20;xmin = 0; xmax = 500;};
    if(varList[j]=="totHtOverPt"     ){     nbins = 20;xmin = 0; xmax = 4;};
    
    if(nbins==1)  cout << "warning : no TH1F definition for variable " << varList[j] << endl;
    TH1F * histo = new TH1F((varList[j]+"_"+channel+"__"+sample).Data(), (varList[j]+"_"+channel+"__"+sample).Data(), nbins, xmin, xmax);
    histo->Sumw2();
    histovect.push_back(histo);
  }

  TH1F * histomTW = new TH1F(("mTW_"+channel+"__"+sample).Data(), ("mTW_"+channel+"__"+sample).Data(), 15,0.,200.);
  histomTW->Sumw2();
  histovect.push_back(histomTW);

  TH1F * histo = new TH1F( ("MVA_"+channel+"__"+sample).Data(), ("MVA_"+channel+"__"+sample).Data(),  20, -1., 1.);
  histo->Sumw2();
  histovect.push_back(histo);
  theHistoMap[sample] = histovect;
  
}


void theMVAtool::fillHisto(TString sample, float* theVar, double mva, double mtw, double weight){
   
   //  cout << "sample " << sample << endl;
   std::vector<TH1F*> histovect = theHistoMap[sample];
   //cout << " histovect.size() " <<  histovect.size()  << endl;
   for(unsigned int i=0; i<  varList.size(); i++) histovect[i]->Fill(theVar[i], weight);
   histovect[varList.size()]->Fill(mtw,weight);
   histovect.back()->Fill(mva, weight);
   
}

void theMVAtool::scaleHisto(TString sample,double thescale){
  
  std::vector<TH1F*> histovect = theHistoMap[sample];
  for(unsigned int i=0; i< (varList.size()+1); i++) histovect[i]->Scale(thescale);
  
}



void theMVAtool::writeHisto(TString sample, TString syst, TString reg){
  if(theOutputFileMap[sample] == 0 ) cout << "no output file for sample " << sample << endl;
  theOutputFileMap[sample]->cd();
  std::vector<TH1F*> histovect = theHistoMap[reg+sample+syst];
  for(unsigned int i=0; i< (histovect.size()); i++) { 
    if(histovect[i] == 0) cout << "no histogram existing for " << reg+sample+syst << endl;
/*    if (reg != "")
    {
        std::string name = histovect[i]->GetName();
        boost::replace_last(name, reg.Data(), "");
        histovect[i]->SetName(name.c_str());
    }*/
    histovect[i]->Write();
  }
  
}

void theMVAtool::makePseudoData(TString inDir, TString channel, TString region) {

  TRandom3 therand(0); //Randomization

  TString pseudodata_output_name = "pseudoData/histoFile_"+region+"_"+channel+"_pseudoData.root";
  TFile *output_file = new TFile (pseudodata_output_name, "RECREATE");
  TTree *output_tree = new TTree("Ttree_pseudoData", "Ttree_pseudoData");

  std::cout << "\n--- GENERATION OF PSEUDODATA IN " << output_file->GetName() << " ! ---\n" << std::endl;

  TH1F *h_sum = 0;

  //file->ls(); //output the content of the file

  for(uint iVar=0; iVar<varList.size(); iVar++) {

    std::cout << " --- " << varList[iVar] << std::endl;

    std::vector<float> sum;

    // get the sum of all the samples to be used as expectation value of the poisson
    for(unsigned int i=0; i< samplelist.size(); i++){
      //cout << samplelist[i] << endl;

      if ( samplelist[i].Contains("Data") ) continue; //From MC only
      TFile *input         = new TFile( (inDir+"histofile_"+samplelist[i]+".root").Data(), "read");
      TTree *input_tree    = (TTree*)input->Get("Ttree_"+region+"_"+samplelist[i]);

      Float_t var;
      TString varName = varList[iVar];
      std::cout << __LINE__<< std::endl;
      input_tree->SetBranchAddress(varName, &var);
      std::cout<<__LINE__<<std::endl;
      for (Int_t nEntry = 0; nEntry < input_tree->GetEntries(); nEntry++) {
	input_tree->GetEntry(nEntry);
	sum[nEntry] += var;
      }
    }
  
    int nofbins = sum.size();
  
    for(int i=0; i<nofbins; i++)
    {
      int bin_content = sum[i]; // std::cout << "initial content = " << bin_content << std::endl;
      int new_bin_content = therand.Poisson(bin_content); //std::cout<<"new content = " << new_bin_content << std::endl;
      sum[i] = new_bin_content;
    }

    // Add new branch to tree
      std::cout << __LINE__<< std::endl;
    TBranch *newBranch = output_tree->Branch(varList[iVar], &sum[0], ( varList[iVar]+"/F") );
      std::cout << __LINE__<< std::endl;
    newBranch->Fill();
      std::cout << __LINE__<< std::endl;
  }

  output_file->cd();
  output_tree->Write();
  output_file->Close();

  std::cout << "--- Done with generation of pseudo-data" << std::endl;  
}
