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
   
   
   theMVAtool(bool doCtrlReg = false);
   theMVAtool(std::vector<TString > thevarlist, std::vector<TString > thesamplelist, std::vector<TString > thesystlist);
   ~theMVAtool(){};
   
   void doTraining(TString channel = "all",TString inDir = "inputroot/met0mtw0/", bool sigMode = false);
   void doReading(float bdtcut = 100., TString channel = "all", TString inDir = "inputroot/met0mtw0/", TString outDir = "outputroot/", bool usePseudoData = true);
   
   void loopInSample(TFile* input, TString sample, float *treevars, float bdtcut, TString channel);
   
   //--------------------
   //deal with histograms
   //--------------------
   void createHisto(TString sample, TString channel);
   void fillHisto(TString sample,  float* theVar, double mva, double mtw, double weight=1);
   void writeHisto(TString sample, TString syst, TString reg);
   void scaleHisto(TString sample, double thescale);
   
   std::vector<TString > varList;
   std::vector<TString > samplelist;
   std::vector<TString > systlist;
   std::vector<TString > regList;
   
   void makePseudoData  ( TString inDir, TString channel = "all", TString region = "" );
  
   std::vector<TH1F*> histovect;
   
   std::map<TString,std::vector<TH1F*> > theHistoMap;
   std::map<TString, int > theVarMap;
   
   std::map<TString, TFile* > theOutputFileMap;
   
   //to apply scale factor when reading
   std::vector<double > sf_DY, sf_DY_err;
   std::vector<double > sf_WZ, sf_WZ_err;
   
   TMVA::Reader *reader; 
   
};



#endif
