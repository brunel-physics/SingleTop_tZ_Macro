#define theMVAtool_cxx
#include "theMVAtool.h"

theMVAtool::theMVAtool(){
  
  //constructor

  varList.push_back("tree_topMass");     
  varList.push_back("tree_Zpt");         
  varList.push_back("tree_Zeta");        
  varList.push_back("tree_topPt");       
  varList.push_back("tree_topEta");      
  varList.push_back("tree_NJets");       
  varList.push_back("tree_NBJets");	 
  varList.push_back("tree_btagDiscri");
  //  varList.push_back("tree_totPt");	
  //  varList.push_back("tree_totEta");	
  varList.push_back("tree_leptWPt");
  varList.push_back("tree_leptWEta");
  varList.push_back("tree_leadJetPt");   
  varList.push_back("tree_leadJetEta");  
  //  varList.push_back("tree_met" );
  //  varList.push_back("tree_mTW" );

  //  User Added
  //  varList.push_back("tree_totPtVec");
  //  varList.push_back("tree_totVecM");
  //  varList.push_back("tree_lepPt");
  varList.push_back("tree_lepMetPt");
  //  varList.push_back("tree_totPt2Jet");
  //  varList.push_back("tree_leadJetbTag");
  varList.push_back("tree_secJetbTag");
  varList.push_back("tree_secJetPt");
  varList.push_back("tree_secJetEta");
  /*  varList.push_back("tree_wzdelR");
  varList.push_back("tree_jjdelR");
  varList.push_back("tree_zjminR");
  varList.push_back("tree_ZlepWdelPhi");
  varList.push_back("tree_ZmetdelPhi");
  varList.push_back("tree_ZlepWdelR");
  varList.push_back("tree_lbDelR");
  varList.push_back("tree_lbDelPhi");
  varList.push_back("tree_zlb1DelR");
  varList.push_back("tree_zlb1DelPhi");
  varList.push_back("tree_zlb2DelR");
  varList.push_back("tree_zlb2DelPhi");
  varList.push_back("tree_totHt");
  varList.push_back("tree_lepHt");
  varList.push_back("tree_jetHt");
  varList.push_back("tree_lepMetHt");
  varList.push_back("tree_totHtOverPt");
  varList.push_back("tree_zMass");*/

  //  End user added varList
 
  for(unsigned int i=0; i< varList.size() ; i++) theVarMap[varList[i]] = i;
  
  samplelist.push_back("DataMu");
  samplelist.push_back("DataMuEG");
  samplelist.push_back("DataEG");
  samplelist.push_back("DataEGZenriched");
  samplelist.push_back("DataMuEGZenriched");
  samplelist.push_back("DataMuZenriched");
  samplelist.push_back("tZq");
  samplelist.push_back("TTZ");
  samplelist.push_back("TTW");
  samplelist.push_back("TT");
  samplelist.push_back("DYToLL_M10-50");
  samplelist.push_back("Zjets");
  samplelist.push_back("Wjets");
  samplelist.push_back("WZ");
  samplelist.push_back("ZZ");
  samplelist.push_back("TbarsChan");
   samplelist.push_back("TsChan");
   samplelist.push_back("TtChan");
  samplelist.push_back("TbartChan");
  samplelist.push_back("TtW");
  samplelist.push_back("TbartW");
  //  samplelist.push_back("WW");
  
  systlist.push_back("");
  systlist.push_back("__lept__plus");
  systlist.push_back("__lept__minus");
  systlist.push_back("__trig__plus");
  systlist.push_back("__trig__minus");
  systlist.push_back("__PDF__plus");
  systlist.push_back("__PDF__minus");
  systlist.push_back("__jes__plus");
  systlist.push_back("__jes__minus");
  systlist.push_back("__jer__plus"); 
  systlist.push_back("__jer__minus");
  systlist.push_back("__metuncls__plus");
  systlist.push_back("__metuncls__minus");
			 
  sf_DY.push_back(0.057); sf_DY_err.push_back(0.00057);
  sf_DY.push_back(0.612); sf_DY_err.push_back(0.0051);
  sf_DY.push_back(0.0637); sf_DY_err.push_back(0.0005);
  sf_DY.push_back(9.216); sf_DY_err.push_back(0.0011);
  
  sf_WZ.push_back(1.153); sf_WZ_err.push_back(0.04);
  sf_WZ.push_back(1.267); sf_WZ_err.push_back(0.06);
  sf_WZ.push_back(1.928); sf_WZ_err.push_back(0.20);
  sf_WZ.push_back(1.135); sf_WZ_err.push_back(0.04);

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

void theMVAtool::doTraining(TString channel){

 //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  TString outfileName( ("trainingBDT_"+channel+"_tZq_.root").Data() );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_"+channel+"_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
  
  
  
  TFile *input_sig        = TFile::Open( "inputroot/histofile_tZq.root" );

  TFile *input_wz         = TFile::Open( "inputroot/histofile_WZ.root" );
  TFile *input_wJets      = TFile::Open( "inputroot/histofile_Wjets.root");
  TFile *input_ww         = TFile::Open( "inputroot/histofile_WW.root");
  
  TFile *input_TTW        = TFile::Open( "inputroot/histofile_TTW.root" );
  TFile *input_TTZ        = TFile::Open( "inputroot/histofile_TTZ.root" );
  TFile *input_TT         = TFile::Open( "inputroot/histofile_TT.root" );
  /* 
  TFile *input_TbarsChan  = TFile::Open( "inputroot/histofile_TbarsChan.root");
  TFile *input_TbartChan  = TFile::Open( "inputroot/histofile_TbartChan.root");
  TFile *input_TbartW     = TFile::Open( "inputroot/histofile_TbartW.root");
  TFile *input_TsChan     = TFile::Open( "inputroot/histofile_TsChan.root");
  TFile *input_TtChan     = TFile::Open( "inputroot/histofile_TtChan.root");
  */
  TFile *input_Zjets      = TFile::Open( "inputroot/histofile_Zjets.root");
  TFile *input_ZZ         = TFile::Open( "inputroot/histofile_ZZ.root");
  
  TFile *input_DY1        = TFile::Open( "inputroot/histofile_DataEGZenriched.root" );
  TFile *input_DY2        = TFile::Open( "inputroot/histofile_DataMuEGZenriched.root" );
  TFile *input_DY3        = TFile::Open( "inputroot/histofile_DataMuZenriched.root" );
  //  TFile *input_DY4        = TFile::Open( "inputroot/histofile_DYToLL_M10-50.root" );

  
  TTree *signal              = (TTree*)input_sig->Get("Ttree_tZq");

  TTree *background_WZ       = (TTree*)input_wz->Get("Ttree_WZ");
  TTree *background_wJets      = (TTree*)input_wJets->Get("Ttree_Wjets.root");
  //  TTree *background_ww         = (TTree*)input_ww->Get("Ttree_WW.root");

  TTree *background_TT      = (TTree*)input_TT->Get("Ttree_TT");
  TTree *background_TTW     = (TTree*)input_TTW->Get("Ttree_TTW");
  TTree *background_TTZ     = (TTree*)input_TTZ->Get("Ttree_TTZ");
  /*
  TTree *background_TbarsChan  = (TTree*)input_TbarsChan->Get("Ttree_TbarsChan.root");
  TTree *background_TbartChan  = (TTree*)input_TbartChan->Get("Ttree_TbartChan.root");
  TTree *background_TbartW     = (TTree*)input_TbartW->Get("Ttree_TbartW.root");
  TTree *background_TsChan     = (TTree*)input_TsChan->Get("Ttree_TsChan.root");
  TTree *background_TtChan     = (TTree*)input_TtChan->Get("Ttree_TtChan.root");
  */
  TTree *background_Zjets = (TTree*)input_Zjets->Get("Ttree_Zjets");
  TTree *background_ZZ    = (TTree*)input_ZZ->Get("Ttree_ZZ");
  
  TTree *background_DY1     = (TTree*)input_DY1->Get("Ttree_DataEGZenriched");
  TTree *background_DY2     = (TTree*)input_DY2->Get("Ttree_DataMuEGZenriched");
  TTree *background_DY3     = (TTree*)input_DY3->Get("Ttree_DataMuZenriched");
  //  TTree *background_DY4     = (TTree*)input_DY4->Get("Ttree_DYToLL_M10-50");

  factory->AddSignalTree (signal, 1.0);

  //  factory->AddBackgroundTree (background_WZ,  0.102); 
  factory->AddBackgroundTree (background_WZ,  1.0); 

  //  factory->AddBackgroundTree (background_wJets, 1.2);
//  factory->AddBackgroundTree (background_ww, 0.102);
  
//  factory->AddBackgroundTree (background_TT,  0.35);
//  factory->AddBackgroundTree (background_TTW, 0.39); 
  factory->AddBackgroundTree (background_TTZ, 1.0); 
  //  factory->AddBackgroundTree (background_TTZ, 0.26);     

  //  factory->AddBackgroundTree (background_TbarsChan, 0.2);
  //  factory->AddBackgroundTree (background_TbartChan, 0.3);
  //  factory->AddBackgroundTree (background_TbartW, 0.4);
  //  factory->AddBackgroundTree (background_TsChan, 0.30);
  //  factory->AddBackgroundTree (background_TtChan, 0.30);
  
  //  factory->AddBackgroundTree (background_Zjets, 2.5);
  //  factory->AddBackgroundTree (background_ZZ, 0.0);
  
  //  factory->AddBackgroundTree (background_DY1,  0.01);
  //  factory->AddBackgroundTree (background_DY2,  0.10);
  //  factory->AddBackgroundTree (background_DY3,  0.20);
  //  factory->AddBackgroundTree (background_DY4, 2.50);
  
 
  for(unsigned int i=0; i< varList.size() ; i++) factory->AddVariable( varList[i].Data(),    'F');
  
  factory->SetSignalWeightExpression    ("tree_EvtWeight");
  factory->SetBackgroundWeightExpression("tree_EvtWeight");
 
  // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";
   
   if(channel !="all"){
     if     (channel == "mumumu") mycuts = "tree_Channel == 0";
     else if(channel == "mumue" ) mycuts = "tree_Channel == 1";
     else if(channel == "eemu"  ) mycuts = "tree_Channel == 2";
     else if(channel == "eee"   ) mycuts = "tree_Channel == 3";
     else cout << "WARNING wrong channel name while training " << endl;
   }
   
   
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   //for WZ
   //factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
   factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   //Default: NTrees= 100,nEventsMin=100,MaxDepth=3
   //Duncan's copy: NTrees= 800,nEventsMin=100,MaxDepth=2

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

void theMVAtool::doReading(float bdtcut, TString channel){
     
  // This loads the library
  TMVA::Tools::Instance();
  
  //create the BDT reader
  reader = new TMVA::Reader( "!Color:!Silent" );  
  unsigned int varsize = varList.size(); 
  float  treevars[varsize];
  for(unsigned int i=0; i<varsize; i++) treevars[i] = 0;  
  for(unsigned int i=0; i< varList.size() ; i++) reader->AddVariable( varList[i].Data(),  &(treevars[i])  );
   
  //  reader->BookMVA( "BDT", ("weights/BDT_trainning_"+channel+"_FCNCzct_BDT.weights.xml").Data() );
  reader->BookMVA( "BDT", ("weights/BDT_trainning_"+channel+"_tzq_BDT.weights.xml").Data() );
  
  // cout << __LINE__ << "   SAMPLE SIZE " << samplelist.size() << endl;
  for(unsigned int i=0; i< samplelist.size(); i++){
    // cout << "samplelist[" << i << "]: " << samplelist[i] << endl;
    TFile *input         = new TFile( ("inputroot/histofile_"+samplelist[i]+".root").Data(), "read");
    TFile * theoutputfile = new TFile(   ("outputroot/output_"+channel+"_"+samplelist[i]+".root").Data(), "recreate");
    theOutputFileMap[samplelist[i]] = theoutputfile;
    

    //---------------------
    //loop over systematics
    //---------------------
    cout <<  "processing sample " <<  samplelist[i]  << endl;
   
    for(unsigned int j=0; j<systlist.size(); j++){
      //    cout << __LINE__  << "  " << samplelist[i]+systlist[j] << endl;
       if(
         (samplelist[i] == "DataMu" || samplelist[i] == "DataEG" || samplelist[i] == "DataMuEG" ||
         samplelist[i] == "DataMuZenriched" || samplelist[i] == "DataEGZenriched" || samplelist[i] == "DataMuEGZenriched") 
         && systlist[j] != "" ) continue;
      createHisto(samplelist[i]+systlist[j], channel);
      loopInSample(input, samplelist[i]+systlist[j], treevars, bdtcut, channel);
      writeHisto(samplelist[i], systlist[j]);

    }
   
    delete theoutputfile;
  }
  
  
  
  
} 

void theMVAtool::loopInSample(TFile* input, TString sample, float *treevars, float bdtcut, TString channel){
  
  
  
  unsigned int varsize = varList.size();
  
  input->cd();
  // cout << "sample " << sample << endl;
  // cout << "tree name " << "Ttree_"+sample  << endl;
  TTree* theTree = (TTree*)input->Get("Ttree_"+sample);
  for (unsigned int ivar=0; ivar<varsize; ivar++) theTree->SetBranchAddress( varList[ivar].Data(), &(treevars[ivar]) );
  float theweight = 0;
  theTree->SetBranchAddress( "tree_EvtWeight", &theweight );
  int theChannel = -1;
  theTree->SetBranchAddress( "tree_Channel", &theChannel );
  
  int theSelChannel = -1;
  if(     channel == "mumumu") theSelChannel = 0;
  else if(channel == "mumue" ) theSelChannel = 1;
  else if(channel == "eemu"  ) theSelChannel = 2;
  else if(channel == "eee"   ) theSelChannel = 3;
  

  if(theTree == 0) cout << "no TTree found with name " << "Ttree_"+sample << endl;
  for(int i=0; i< theTree->GetEntries(); i++){
    theTree->GetEntry(i);
    
    double sf_local = 1.;
    
    
    //cout << "----------------------" << endl;
    //cout << "channel " << channel << endl;
    //cout << "theChannel    " << theChannel << endl;
    //cout << "theSelChannel " <<  theSelChannel<< endl;
    if(channel != "all" && theChannel != theSelChannel ) continue;
    
    //cout << "pass sel " << endl;
  
    if(sample == "WZ" ){
      //      cout << "theChannel " << theChannel << endl;
      if(theChannel == 0) sf_local = sf_WZ[0]; 
      if(theChannel == 1) sf_local = sf_WZ[1]; 
      if(theChannel == 2) sf_local = sf_WZ[2]; 
      if(theChannel == 3) sf_local = sf_WZ[3]; 
    } 
    
    if(sample == "Zjets" || sample ==  "DYToLL_M10-50"){
      if(theChannel == 0) sf_local = sf_DY[0]; 
      if(theChannel == 1) sf_local = sf_DY[1]; 
      if(theChannel == 2) sf_local = sf_DY[2]; 
      if(theChannel == 3) sf_local = sf_DY[3]; 
    } 
    //if(sample == "WZ") cout << "weight 1 " << theweight << endl;
    theweight*=sf_local;
    //if(sample == "WZ") cout << "weight 2 " << theweight << endl;
    
    double mvaValue = reader->EvaluateMVA( "BDT");
    //    cout << reader->EvaluateMVA( "BDT") << endl;
    if(mvaValue > bdtcut) continue;
    //    cout << "the bdt " << bdtcut << endl;
    fillHisto(sample, treevars, mvaValue, theweight);
  }
  
  
  
  
}






void theMVAtool::createHisto(TString sample, TString channel){
  std::vector<TH1F*> histovect;
  for(unsigned int j=0; j<varList.size(); j++){  
    
    int nbins=1;
    double xmin = -1000;
    double xmax = 1000;
    if(varList[j]=="tree_topMass")       {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_Zpt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_Zeta")	   {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_topPt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_topEta")        {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_btagDiscri")    {nbins = 10; xmin = 0;   xmax = 1;};
    //    if(varList[j]=="tree_totPt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    //    if(varList[j]=="tree_totEta")        {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_leptWPt")       {nbins = 20; xmin = 0;   xmax = 350;};
    if(varList[j]=="tree_leptWEta")      {nbins = 10; xmin = -2.5;xmax = 2.5;};    
    if(varList[j]=="tree_leadJetPt")     {nbins = 20; xmin = 0;   xmax = 500;};	      
    if(varList[j]=="tree_leadJetEta")    {nbins = 20; xmin = -2.5; xmax = 2.5;};	       
    //    if(varList[j]=="tree_met"){     nbins = 10; xmin = 0;   xmax = 200;};	   
    //    if(varList[j]=="tree_mTW"){     nbins = 10; xmin = 0;   xmax = 200;};
    if(varList[j]=="tree_NJets"){     nbins = 10; xmin = 0;   xmax = 200;};
    if(varList[j]=="tree_NBJets"){     nbins = 10; xmin = 0;   xmax = 200;};
    //    if(varList[j]=="tree_totPtVec"){nbins = 20; xmin = 0;   xmax = 500;};
    //    if(varList[j]=="tree_totVecM"){nbins = 20; xmin = 0;   xmax = 500;};
    //    if(varList[j]=="tree_lepPt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_lepMetPt") {nbins = 20; xmin = 0;   xmax = 500;};
    //    if(varList[j]=="tree_totPt2Jet") {nbins = 20; xmin = 0;   xmax = 500;};
    //    if(varList[j]=="tree_leadJetbTag") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_secJetbTag") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_secJetPt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_secJetEta") {nbins = 20; xmin = 0;   xmax = 500;};
    /*    if(varList[j]=="tree_wzdelR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_jjdelR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_zjminR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_ZlepWdelPhi") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_ZmetdelPhi") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_ZlepWdelR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_lbDelR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_lbDelPhi") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_zlb1DelR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_zlb1DelPhi") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_zlb2DelR") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_zlb2DelPhi") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_totHt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_lepHt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_jetHt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_lepMetHt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_totHtOverPt") {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_zMass") {nbins = 20; xmin = 0;   xmax = 500;};*/
    
    if(nbins==1)  cout << "warning : no TH1F definition for variable " << varList[j] << endl;
    TH1F * histo = new TH1F((varList[j]+"_"+channel+"__"+sample).Data(), (varList[j]+"_"+channel+"__"+sample).Data(), nbins, xmin, xmax);
    histo->Sumw2();
    histovect.push_back(histo);
  }
  
  
  TH1F * histo = new TH1F( ("MVA_"+channel+"__"+sample).Data(), ("MVA_"+channel+"__"+sample).Data(),  20, -1, 1);
  histo->Sumw2();
  histovect.push_back(histo);
  theHistoMap[sample] = histovect;
  
}


void theMVAtool::fillHisto(TString sample, float* theVar, double mva, double weight){
   
   //cout << "sample " << sample << endl;
   std::vector<TH1F*> histovect = theHistoMap[sample];
   //cout << " histovect.size() " <<  histovect.size()  << endl;
   for(unsigned int i=0; i<  varList.size(); i++) histovect[i]->Fill(theVar[i], weight);
   histovect.back()->Fill(mva, weight);
   
}

void theMVAtool::scaleHisto(TString sample,double thescale){
  
  std::vector<TH1F*> histovect = theHistoMap[sample];
  for(unsigned int i=0; i< (varList.size()+1); i++) histovect[i]->Scale(thescale);
  
}



void theMVAtool::writeHisto(TString sample, TString syst){
  if(theOutputFileMap[sample] == 0 ) cout << "no output file for sample " << sample << endl;
  theOutputFileMap[sample]->cd();
  std::vector<TH1F*> histovect = theHistoMap[sample+syst];
  for(unsigned int i=0; i< (varList.size()+1); i++) { 
    if(histovect[i] == 0) cout << "no histogram existing for " << sample+syst << endl;
    histovect[i]->Write();
  }
  
}


