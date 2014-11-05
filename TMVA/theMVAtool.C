#define theMVAtool_cxx
#include "theMVAtool.h"

theMVAtool::theMVAtool(){
  
  //constructor
  ////varList.push_back("tree_cosThetaStar");
  varList.push_back("tree_topMass");     
  ////varList.push_back("tree_totMass");     
  ////varList.push_back("tree_deltaPhilb");  
  ////varList.push_back("tree_deltaRlb");    
  //varList.push_back("tree_deltaRTopZ");  
  varList.push_back("tree_asym");        
  //varList.push_back("tree_Zpt");         
  //varList.push_back("tree_ZEta");        
  //varList.push_back("tree_topPt");       
  //varList.push_back("tree_topEta");      
  //varList.push_back("tree_NJets");       
  varList.push_back("tree_NBJets");	 
  //varList.push_back("tree_deltaRZl");	 
  //varList.push_back("tree_deltaPhiZmet");
  varList.push_back("tree_btagDiscri");  
  varList.push_back("tree_totPt");	
  //varList.push_back("tree_totEta");	
  varList.push_back("tree_leptWPt");	 
  //varList.push_back("tree_leptWEta");	 
  varList.push_back("tree_leadJetPt");   
  //varList.push_back("tree_leadJetEta");  
  //varList.push_back("tree_deltaRZleptW");
  //varList.push_back("tree_deltaPhiZleptW");
  //varList.push_back("tree_met" );
  //varList.push_back("tree_mTW" );
  
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
  samplelist.push_back("WZHF");
  samplelist.push_back("ZZ");
  samplelist.push_back("TbarsChan");
  samplelist.push_back("TsChan");
  samplelist.push_back("TtChan");
  samplelist.push_back("TbartChan");
  samplelist.push_back("TtW");
  samplelist.push_back("TbartW");
  //samplelist.push_back("WW");
  
  systlist.push_back("");  
  systlist.push_back("__lept__plus");
  systlist.push_back("__lept__minus");
  //systlist.push_back("__trig__plus");
  //systlist.push_back("__trig__minus");
  //systlist.push_back("__PDF__plus");
  //systlist.push_back("__PDF__minus");*/
  systlist.push_back("__jes__plus");
  systlist.push_back("__jes__minus");
  systlist.push_back("__jer__plus"); 
  systlist.push_back("__jer__minus");
  systlist.push_back("__metuncls__plus");
  systlist.push_back("__metuncls__minus"); 
  
  
  
  sf_DY.push_back(0.77); sf_DY_err.push_back(0.64);
  sf_DY.push_back(0.97); sf_DY_err.push_back(0.21); 
  sf_DY.push_back(2.89); sf_DY_err.push_back(1.27);
  sf_DY.push_back(1.39); sf_DY_err.push_back(0.23); 
  
  
  sf_WZ.push_back(1.01); sf_WZ_err.push_back(0.04);
  sf_WZ.push_back(0.87); sf_WZ_err.push_back(0.05); 
  sf_WZ.push_back(0.67); sf_WZ_err.push_back(0.04);
  sf_WZ.push_back(1.14); sf_WZ_err.push_back(0.07);
  

  
  
  
  
} 


theMVAtool::theMVAtool(std::vector<TString > thevarlist, std::vector<TString > thesamplelist, std::vector<TString > thesystlist){
  for(unsigned int i=0; i<thevarlist.size(); i++) { varList.push_back(thevarlist[i]);theVarMap[varList[i]] = i;}
  for(unsigned int i=0; i<thesamplelist.size(); i++) { samplelist.push_back(thesamplelist[i]);}
  for(unsigned int i=0; i<thesystlist.size(); i++) { systlist.push_back(thesystlist[i]);}
} 




//---------------------------------------------------------------
 // to perform the training
//---------------------------------------------------------------

void theMVAtool::doTraining(){

 //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  TString outfileName( "trainingBDT_tZq.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  //TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  TMVA::Factory *factory = new TMVA::Factory( "BDT_trainning_tzq", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
  
  
  
  TFile *input_sig        = TFile::Open( "../TreeReader/outputroot/histofile_tZq.root" );
  TFile *input_wz         = TFile::Open( "../TreeReader/outputroot/histofile_WZ.root" );
  TFile *input_wzHF       = TFile::Open( "../TreeReader/outputroot/histofile_WZHF.root" );
  
  
  /*TFile *input_TTW        = TFile::Open( "../TreeReader/outputroot/histofile_TTW.root" );
  TFile *input_TTZ        = TFile::Open( "../TreeReader/outputroot/histofile_TTZ.root" );
  TFile *input_TT        = TFile::Open( "../TreeReader/outputroot/histofile_TT.root" );
  TFile *input_DY1        = TFile::Open( "../TreeReader/outputroot/histofile_DataMuZenriched.root" );
  TFile *input_DY2        = TFile::Open( "../TreeReader/outputroot/histofile_DataEGZenriched.root" );
  TFile *input_DY3        = TFile::Open( "../TreeReader/outputroot/histofile_DataMuZenriched.root" );*/
  
  
  TTree *signal              = (TTree*)input_sig->Get("Ttree_tZq");
  TTree *background_WZ       = (TTree*)input_wz->Get("Ttree_WZ");
  TTree *background_WZHF     = (TTree*)input_wzHF->Get("Ttree_WZHF");
  /*TTree *background_TT      = (TTree*)input_TT->Get("Ttree_TT");
  TTree *background_TTW     = (TTree*)input_TTW->Get("Ttree_TTW");
  TTree *background_TTZ     = (TTree*)input_TTZ->Get("Ttree_TTZ");*/
  
  
  
  //TTree *background_DY1     = (TTree*)input_DY1->Get("Ttree_DataMuZenriched");
  //TTree *background_DY2     = (TTree*)input_DY2->Get("Ttree_DataEGZenriched");
  //TTree *background_DY3     = (TTree*)input_DY3->Get("Ttree_DataMuZenriched");
  
  factory->AddSignalTree    ( signal,	   1.);
  factory->AddBackgroundTree( background_WZ,  0.102);
  factory->AddBackgroundTree( background_WZHF,  0.898);
  /*factory->AddBackgroundTree    (background_TT,        0.35);
  factory->AddBackgroundTree    (background_TTW,       0.39);
  factory->AddBackgroundTree    (background_TTZ,       0.26);*/
  
  
  //factory->AddBackgroundTree( background_DY1,  0.01);
  //factory->AddBackgroundTree( background_DY2,  0.10);
  //factory->AddBackgroundTree( background_DY3,  0.20);
  
 
  for(unsigned int i=0; i< varList.size() ; i++) factory->AddVariable( varList[i].Data(),    'F');
  
  factory->SetSignalWeightExpression    ("tree_EvtWeight");
  factory->SetBackgroundWeightExpression("tree_EvtWeight");
 
  // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   //for WZ
   //factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
   factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=100:nEventsMin=100:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

 

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

void theMVAtool::doReading(float bdtcut){
  
  
   
  // This loads the library
  TMVA::Tools::Instance();
  
  //create the BDT reader
  reader = new TMVA::Reader( "!Color:!Silent" );  
  unsigned int varsize = varList.size();
  float  treevars[varsize];
  for(unsigned int i=0; i<varsize; i++) treevars[i] = 0;  
  
  for(unsigned int i=0; i< varList.size() ; i++) reader->AddVariable( varList[i].Data(),  &(treevars[i])  );
  
   
  reader->BookMVA( "BDT", "weights/BDT_trainning_tzq_BDT.weights.xml" ); 
  
 // cout << __LINE__ << "   SAMPLE SIZE " << samplelist.size() << endl;
  for(unsigned int i=0; i< samplelist.size(); i++){
    //cout << samplelist[i] << endl;
    TFile *input         = new TFile( ("../TreeReader/outputroot/histofile_"+samplelist[i]+".root").Data(), "read");
    TFile * theoutputfile = new TFile(   ("outputroot/output_"+samplelist[i]+".root").Data(), "recreate");
    theOutputFileMap[samplelist[i]] = theoutputfile;
    
    
    //---------------------
    //loop over systematics
    //---------------------
   cout <<  "processing sample " <<  samplelist[i]  << endl;
   
   
    for(unsigned int j=0; j<systlist.size(); j++){
     // cout << __LINE__  << "  " << samplelist[i]+systlist[j] << endl;
       if(
         (samplelist[i] == "DataMu" || samplelist[i] == "DataEG" || samplelist[i] == "DataMuEG" ||
         samplelist[i] == "DataMuZenriched" || samplelist[i] == "DataEGZenriched" || samplelist[i] == "DataMuEGZenriched") 
         && systlist[j] != "" ) continue;
      createHisto(samplelist[i]+systlist[j]);
      loopInSample(input, samplelist[i]+systlist[j], treevars, bdtcut);
      writeHisto(samplelist[i], systlist[j]);
    }
   
    delete theoutputfile;
  }
  
  
  
  
}

void theMVAtool::loopInSample(TFile* input, TString sample, float *treevars, float bdtcut){
  
  
  
  unsigned int varsize = varList.size();
  
  
  input->cd();
  //cout << "sample " << sample << endl;
  TTree* theTree = (TTree*)input->Get("Ttree_"+sample);
  for (unsigned int ivar=0; ivar<varsize; ivar++) theTree->SetBranchAddress( varList[ivar].Data(), &(treevars[ivar]) );
  float theweight = 0;
  theTree->SetBranchAddress( "tree_EvtWeight", &theweight );
  int theChannel = -1;
  theTree->SetBranchAddress( "tree_Channel", &theChannel );
  
  if(theTree == 0) cout << "no TTree found with name " << "Ttree_"+sample << endl;
  for(int i=0; i< theTree->GetEntries(); i++){
    theTree->GetEntry(i);
    
    double sf_local = 1.;
  
    if(sample == "WZ" ){
      //cout << "theChannel " << theChannel << endl;
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
    //cout << reader->EvaluateMVA( "BDT") << endl;
    if(mvaValue > bdtcut) continue;
    //cout << "the bdt " << bdtcut << endl;
    fillHisto(sample, treevars, mvaValue, theweight);
  }
  
  
  
  
}






void theMVAtool::createHisto(TString sample){
  std::vector<TH1F*> histovect;
  for(unsigned int j=0; j<varList.size(); j++){  
    
    int nbins=1;
    double xmin = -1000;
    double xmax = 1000;
    if(varList[j]=="tree_cosThetaStar")  {nbins = 20; xmin = -1;  xmax = 1;};
    if(varList[j]=="tree_topMass")       {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_totMass")       {nbins = 20; xmin = 100; xmax = 700;};
    if(varList[j]=="tree_deltaPhilb")    {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="tree_deltaRlb")      {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="tree_deltaRTopZ")    {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="tree_asym")	   {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_Zpt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_ZEta")	   {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_topPt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_totPt")	   {nbins = 20; xmin = 0;   xmax = 500;};
    if(varList[j]=="tree_topEta")        {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_totEta")        {nbins = 10; xmin = -2.5;xmax = 2.5;};
    if(varList[j]=="tree_NJets")	   {nbins = 4;  xmin = 0.5; xmax = 4.5;};      
    if(varList[j]=="tree_NBJets")        {nbins = 2;  xmin = -0.5;xmax = 1.5;};  
    if(varList[j]=="tree_deltaRZl")      {nbins = 10; xmin = 0;   xmax = 3.2;}; 
    if(varList[j]=="tree_deltaPhiZmet")  {nbins = 10; xmin = 0;   xmax = 3.2;};
    if(varList[j]=="tree_btagDiscri")    {nbins = 10; xmin = 0;   xmax = 1;};
    if(varList[j]=="tree_leptWPt")       {nbins = 20; xmin = 0;   xmax = 350;};
    if(varList[j]=="tree_leptWEta")      {nbins = 10; xmin = -2.5;xmax = 2.5;};	     
    if(varList[j]=="tree_leadJetPt")     {nbins = 20; xmin = 0;   xmax = 500;};	      
    if(varList[j]=="tree_leadJetEta")    {nbins = 20; xmin = -2.5; xmax = 2.5;};	       
    if(varList[j]=="tree_deltaPhiZleptW"){nbins = 10; xmin = 0;   xmax = 3.2;};	  	       
    if(varList[j]=="tree_deltaRZleptW"){nbins = 10; xmin = 0;   xmax = 3.2;};	   
    if(varList[j]=="tree_met"	     ){     nbins = 10; xmin = 0;   xmax = 200;};	   
    if(varList[j]=="tree_mTW"	     ){     nbins = 10; xmin = 0;   xmax = 200;};
    
    
    if(nbins==1)  cout << "warning : no TH1F definition for variable " << varList[j] << endl;
    TH1F * histo = new TH1F((varList[j]+"__"+sample).Data(), (varList[j]+"__"+sample).Data(), nbins, xmin, xmax);
    histovect.push_back(histo);
  }
  
  
  TH1F * histo = new TH1F( ("MVA__"+sample).Data(), ("MVA__"+sample).Data(),  20, -1, 1);
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


