#ifndef theMVAtool_h
#define theMVAtool_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>

#include <TLorentzVector.h>

#include <iostream>

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
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"




class theMVAtool {
public :
   
   
   theMVAtool();
   theMVAtool(std::vector<TString > thevarlist, std::vector<TString > thesamplelist, std::vector<TString > thesystlist);
   ~theMVAtool(){};
   
   
   void doTraining();
   void doReading();
   
   void loopInSample(TFile* input, TString sample, float *treevars);
   
   //--------------------
   //deal with histograms
   //--------------------
   void createHisto(TString sample);
   void fillHisto(TString sample,  float* theVar, double mva, double weight=1);
   void writeHisto(TString sample, TString syst);
   void scaleHisto(TString sample, double thescale);
   
   std::vector<TString > varList;
   std::vector<TString > samplelist;
   std::vector<TString > systlist;
   
   
   std::vector<TH1F*> histovect;
   
   std::map<TString,std::vector<TH1F*> > theHistoMap;
   std::map<TString, int > theVarMap;
   
   std::map<TString, TFile* > theOutputFileMap;
   
   
   
   TMVA::Reader *reader; 
   
};



#endif
