#include "TString.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include <iostream>



void Merge_Data(){

   TFile * filechannel_DataMu   = new TFile("outputroot/output_DataMu.root");
   cout << __LINE__ << endl;
   TFile * filechannel_DataEG   = new TFile("outputroot/output_DataEG.root");
   cout << __LINE__ << endl;
   TFile * filechannel_DataMuEG = new TFile("outputroot/output_DataMuEG.root");

   TH1F * histBdt_DataMu   = (TH1F*)filechannel_DataMu->Get( "MVA_DataMu")->Clone();
   cout << __LINE__ << endl;
   TH1F * histBdt_DataEG   = (TH1F*)filechannel_DataEG->Get( "MVA_DataEG")->Clone();
   cout << __LINE__ << endl;
   TH1F * histBdt_DataMuEG = (TH1F*)filechannel_DataMuEG->Get( "MVA_DataMuEG")->Clone();
   cout << __LINE__ << endl;
   
   cout << __LINE__ << endl;
   histBdt_DataMu->Add(histBdt_DataEG);
   cout << __LINE__ << endl;
   histBdt_DataMu->Add(histBdt_DataMuEG);
   cout << __LINE__ << endl;
   
   
   TFile * outputfile = new TFile ("outputroot/datafile.root", "recreate");
   cout << __LINE__ << endl;
   outputfile->cd ();
   
   cout << __LINE__ << endl;
   histBdt_DataMu->Write("MVA__data");
   cout << __LINE__ << endl;
}
