#define theMVAtool_cxx
#include "theMVAtool.h"

theMVAtool::theMVAtool(bool doCtrlReg){
  
  //constructor
  varList.push_back("mTW");
  varList.push_back("wQuark1Pt");
  varList.push_back("wQuark1Eta");
  varList.push_back("wQuark1Phi");
  varList.push_back("wQuark2Pt");
  varList.push_back("wQuark2Eta");
  varList.push_back("wQuark2Phi");
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
  varList.push_back("zLepdelR");
  varList.push_back("zLepdelPhi");
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
  varList.push_back("chan");
  varList.push_back("totPt2Jet");
  varList.push_back("wzdelR");
  varList.push_back("wzdelPhi");
  varList.push_back("zminR");
  varList.push_back("zminPhi");
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
  samplelist.push_back("THQ");
  samplelist.push_back("TTZ");
  samplelist.push_back("TTW");
  samplelist.push_back("TT");
//  samplelist.push_back("DYToLL_M10-50"); // Empty tree - background totally cut
  samplelist.push_back("DYToLL_M50");
//  samplelist.push_back("Wjets"); // Empty tree - background totally cut
//  samplelist.push_back("WW"); // Empty tree - background totally cut
  samplelist.push_back("WZ");
  samplelist.push_back("ZZ");
//  samplelist.push_back("TsChan"); // Empty tree - background totally cut
  samplelist.push_back("TtChan");
  samplelist.push_back("TbartChan");
  samplelist.push_back("TtW");
  samplelist.push_back("TbartW");
  
  systlist.push_back("");  
  systlist.push_back("__pileup__plus");
  systlist.push_back("__pileup__minus");
  systlist.push_back("__trig__plus");
  systlist.push_back("__trig__minus");
  systlist.push_back("__bTag__plus");
  systlist.push_back("__bTag__minus");
  systlist.push_back("__pdf__plus");
  systlist.push_back("__pdf__minus");
  systlist.push_back("__jes__plus");
  systlist.push_back("__jes__minus");
  systlist.push_back("__jer__plus"); 
  systlist.push_back("__jer__minus");
  systlist.push_back("__ME_PS__plus");
  systlist.push_back("__ME_PS__minus"); 
  
  
  
// Trilepton sf's
/*  
  sf_DY.push_back(1.); sf_DY_err.push_back(1.);
  sf_DY.push_back(1.); sf_DY_err.push_back(1.);
  sf_DY.push_back(1.); sf_DY_err.push_back(1.);
  sf_DY.push_back(1.); sf_DY_err.push_back(1.);
  
  sf_WZ.push_back(1.); sf_WZ_err.push_back(1.);
  sf_WZ.push_back(1.); sf_WZ_err.push_back(1.);
  sf_WZ.push_back(1.); sf_WZ_err.push_back(1.);
  sf_WZ.push_back(1.); sf_WZ_err.push_back(1.);
*/
  /* These are from the lepton selection stage calculations
  sf_WZ.push_back(1.204); sf_WZ_err.push_back(0.04);
  sf_WZ.push_back(1.267); sf_WZ_err.push_back(0.06);
  sf_WZ.push_back(1.051); sf_WZ_err.push_back(0.05);
  sf_WZ.push_back(1.167); sf_WZ_err.push_back(0.04);*/
  /*    sf_DY.push_back(0.39); sf_DY_err.push_back(0.63);
    sf_DY.push_back(1.13); sf_DY_err.push_back(0.3); 
    sf_DY.push_back(3.36); sf_DY_err.push_back(1.39);
    sf_DY.push_back(0.754); sf_DY_err.push_back(0.03); 
  
  
    sf_WZ.push_back(1.01); sf_WZ_err.push_back(0.05);
    sf_WZ.push_back(1.16); sf_WZ_err.push_back(0.09); 
    sf_WZ.push_back(0.94); sf_WZ_err.push_back(0.07);
    sf_WZ.push_back(1.651); sf_WZ_err.push_back(0.03); 
  */
} 


theMVAtool::theMVAtool(std::vector<TString > thevarlist, std::vector<TString > thesamplelist, std::vector<TString > thesystlist){
  for(unsigned int i=0; i<thevarlist.size(); i++) { varList.push_back(thevarlist[i]);theVarMap[varList[i]] = i;}
  for(unsigned int i=0; i<thesamplelist.size(); i++) { samplelist.push_back(thesamplelist[i]);}
  for(unsigned int i=0; i<thesystlist.size(); i++) { systlist.push_back(thesystlist[i]);}
} 




//---------------------------------------------------------------
 // to perform the training
//---------------------------------------------------------------

void theMVAtool::doTraining(TString channel, TString inDir){

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
  TFile *input_THQ        = TFile::Open( inDir+"/histofile_THQ.root" );

  TFile *input_TT         = TFile::Open( inDir+"/histofile_TT.root" );  
//  TFile *input_WW         = TFile::Open( inDir+"/histofile_WW.root" ); // Empty tree - background totally cut
  TFile *input_WZ         = TFile::Open( inDir+"/histofile_WZ.root" );
  TFile *input_ZZ         = TFile::Open( inDir+"/histofile_ZZ.root" );

  TFile *input_TtChan     = TFile::Open( inDir+"/histofile_TtChan.root" );
  TFile *input_TbartChan  = TFile::Open( inDir+"/histofile_TbartChan.root" );
//  TFile *input_TsChan     = TFile::Open( inDir+"/histofile_TsChan.root" ); // Empty tree - background totally cut
  TFile *input_TtW        = TFile::Open( inDir+"/histofile_TtW.root" );
  TFile *input_TbartW     = TFile::Open( inDir+"/histofile_TbartW.root" );

  TFile *input_DY50	  = TFile::Open( inDir+"/histofile_DYToLL_M-50.root" );
//  TFile *input_DY10to50	  = TFile::Open( inDir+"/histofile_DYToLL_M10-50.root" ); // Empty tree - background totally cut
  
  TTree *signal              = (TTree*)input_sig->Get("Ttree_tZq");

  TTree *background_TTZ     = (TTree*)input_TTZ->Get("Ttree_TTZ");
  TTree *background_TTW     = (TTree*)input_TTW->Get("Ttree_TTW");
  TTree *background_THQ     = (TTree*)input_TTW->Get("Ttree_THQ");

  TTree *background_TT     = (TTree*)input_TT->Get("Ttree_TT");
//  TTree *background_WW     = (TTree*)input_WW->Get("Ttree_WW");
  TTree *background_WZ     = (TTree*)input_WZ->Get("Ttree_WZ");
  TTree *background_ZZ     = (TTree*)input_ZZ->Get("Ttree_ZZ");
  
  TTree *background_TtChan    = (TTree*)input_TtChan->Get("Ttree_TtChan");
  TTree *background_TbartChan = (TTree*)input_TbartChan->Get("Ttree_TbartChan");
//  TTree *background_TsChan    = (TTree*)input_TsChan->Get("Ttree_TsChan");
  TTree *background_TtW       = (TTree*)input_TtW->Get("Ttree_TtW");  
  TTree *background_TbartW    = (TTree*)input_TbartW->Get("Ttree_TbartW");  

  TTree *background_DY50         = (TTree*)input_DY50->Get("Ttree_DYToLL_M-50");
//  TTree *background_DY10to50   = (TTree*)input_DY10to50->Get("Tree_DYToLL_M10-50");

  factory->AddSignalTree      ( signal, 1. );
  factory->AddBackgroundTree  ( background_TTZ, 1. );
  factory->AddBackgroundTree  ( background_TTW, 1. );
//  factory->AddBackgroundTree  ( background_THQ, 1. );

  factory->AddBackgroundTree  ( background_TT, 1. );
//  factory->AddBackgroundTree  ( background_WW, 1. );
  factory->AddBackgroundTree  ( background_WZ, 1. );
  factory->AddBackgroundTree  ( background_ZZ, 1. );

  factory->AddBackgroundTree  ( background_TtChan, 1. );
  factory->AddBackgroundTree  ( background_TbartChan, 1. );
//  factory->AddBackgroundTree  ( background_TsChan, 1. );
  factory->AddBackgroundTree  ( background_TtW, 1. );
  factory->AddBackgroundTree  ( background_TbartW, 1. );

  factory->AddBackgroundTree  ( background_DY50, 1. );
//  factory->AddBackgroundTree  ( background_DY10to50, 1. );

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
   factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=75:nEventsMin=100:MaxDepth=2:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
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

void theMVAtool::doReading(float bdtcut, TString channel, TString inDir, TString outDir){
  
  
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
  
 // cout << __LINE__ << "   SAMPLE SIZE " << samplelist.size() << endl;
  for(unsigned int i=0; i< samplelist.size(); i++){
    //cout << samplelist[i] << endl;
    TFile *input         = new TFile( (inDir+"histofile_"+samplelist[i]+".root").Data(), "read");
    TFile * theoutputfile = new TFile(   (outDir+"output_"+channel+"_"+samplelist[i]+".root").Data(), "recreate");
    theOutputFileMap[samplelist[i]] = theoutputfile;
    
    
    //---------------------
    //loop over systematics
    //---------------------
    cout <<  "processing sample " <<  samplelist[i]  << endl;
   
   
    for(unsigned int j=0; j<systlist.size(); j++){
     // cout << __LINE__  << "  " << samplelist[i]+systlist[j] << endl;
      for (unsigned int k=0; k<regList.size(); k++){
	if( (samplelist[i] == "DataMu" || samplelist[i] == "DataEG" || samplelist[i] == "DataMuEG") && systlist[j] != "") continue;
	if( (samplelist[i] == "DataMuZenriched" || samplelist[i] == "DataEGZenriched" || samplelist[i] == "DataMuEGZenriched") && !(systlist[j] == "" || systlist[j] == "__zPt__plus" || systlist[j] == "__zPt__minus") ) continue;
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
  // cout << "tree name " << "Ttree_"+sample  << endl;
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
    
    
    /*cout << "----------------------" << endl;
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
    
    if(sample.Contains("Zjets") || sample.Contains("DYToLL_M10-50")){
      if(theChannel == 0) sf_local = sf_DY[0]; 
      if(theChannel == 1) sf_local = sf_DY[1]; 
      if(theChannel == 2) sf_local = sf_DY[2]; 
      if(theChannel == 3) sf_local = sf_DY[3]; 
    } */
    //    if(sample == "TTZ") cout << "weight 1 " << theweight << endl;
    theweight*=sf_local;
    //if(sample == "TTZ") cout << "weight 2 " << theweight << endl;
    
    double mvaValue = reader->EvaluateMVA( "BDT");
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
    if(varList[j]=="cosThetaStar")  {nbins = 20; xmin = -1;  xmax = 1;};
    if(varList[j]=="topMass")       {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="totMass")       {nbins = 20; xmin = 100; xmax = 700;};
    if(varList[j]=="deltaPhilb")    {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="deltaRlb")      {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="deltaRTopZ")    {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="asym")	   {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="Zpt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="Zeta")	   {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="topPt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="totPt")	   {nbins = 10; xmin = 0;   xmax = 250;};
    if(varList[j]=="topEta")        {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="totEta")        {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="NJets")	   {nbins = 6;  xmin = 0.5; xmax = 6.5;};      
    if(varList[j]=="NBJets")        {nbins = 4;  xmin = -0.5;xmax = 3.5;};  
    if(varList[j]=="deltaRZl")      {nbins = 10; xmin = 0;   xmax = 3.2;}; 
    if(varList[j]=="deltaPhiZmet")  {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="btagDiscri")    {nbins = 10; xmin = 0.6;   xmax = 1;};
    if(varList[j]=="leptWPt")       {nbins = 15; xmin = 0;   xmax = 300;};
    if(varList[j]=="leptWEta")      {nbins = 10; xmin = -2.5;xmax = 2.5;};	     
    if(varList[j]=="leadJetPt")     {nbins = 20; xmin = 0;   xmax = 500;};	      
    if(varList[j]=="leadJetEta")    {nbins = 10; xmin = -2.5; xmax = 2.5;};	       
    if(varList[j]=="deltaPhiZleptW"){nbins = 10; xmin = 0;   xmax = 3.2;};	  	       
    if(varList[j]=="deltaRZleptW"){nbins = 10; xmin = 0;   xmax = 3.2;};	   
    if(varList[j]=="met"	     ){     nbins = 10; xmin = 0;   xmax = 200;};	   
    if(varList[j]=="mTW"	     ){     nbins = 10; xmin = 0;   xmax = 200;};
    if(varList[j]=="lepPt"	     ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="lepMetPt"   ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="totPt2Jet"  ){     nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="secJetEta"  ){     nbins = 20; xmin = -5.;   xmax = 5.;};
    if(varList[j]=="secJetPt"  ){     nbins = 15; xmin = 0;   xmax = 300;};
    if(varList[j]=="secJetbTag"  ){     nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="leadJetbTag"  ){     nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="wzdelR"     ){     nbins = 20;xmin = 0; xmax = 7;};
    if(varList[j]=="jjdelR"     ){     nbins = 20;xmin = 0; xmax = 6;};
    if(varList[j]=="zjminR"     ){     nbins = 15;xmin = 0; xmax = 3.5;};
    if(varList[j]=="ZlepWdelPhi"     ){     nbins = 20;xmin = -3.5; xmax = 3.5;};
    if(varList[j]=="ZmetdelPhi"     ){     nbins = 20;xmin = -3.5; xmax = 3.5;};
    if(varList[j]=="ZlepWdelR"     ){     nbins = 20;xmin = 0; xmax = 6;};
    if(varList[j]=="lbDelR"     ){     nbins = 20;xmin = 0; xmax = 6;};
    if(varList[j]=="totHt"     ){     nbins = 20;xmin = 150; xmax = 1500;};
    if(varList[j]=="totVecM"     ){     nbins = 20;xmin = 0; xmax = 300;};
    if(varList[j]=="zlb1DelR"     ){     nbins = 15;xmin = 0; xmax = 6;};
    if(varList[j]=="zlb1DelPhi"     ){     nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="zlb2DelPhi"     ){     nbins = 20;xmin = -4; xmax = 4;};

    if(varList[j]=="zlb2DelR"     ){nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="lbDelPhi"     ){nbins = 20;xmin = -4; xmax = 4;};
    if(varList[j]=="lepHt"     ){nbins = 15;xmin = 0; xmax = 4;};
    if(varList[j]=="jetHt"     ){nbins = 20;xmin = 0; xmax = 500;};
    if(varList[j]=="lepMetHt"     ){nbins = 20;xmin = 0; xmax = 500;};
    if(varList[j]=="totHtOverPt"     ){nbins = 20;xmin = 0; xmax = 4;};
    if(varList[j]=="totPtVec"     ){nbins = 20;xmin = 0; xmax = 200;};
    if(varList[j]=="wLepjminR"     ){     nbins = 15;xmin = 0; xmax = 3.5;};

    
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
   
  //   cout << "sample " << sample << endl;
   std::vector<TH1F*> histovect = theHistoMap[sample];
   // cout << " histovect.size() " <<  histovect.size()  << endl;
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
    histovect[i]->Write();
  }
  
}


