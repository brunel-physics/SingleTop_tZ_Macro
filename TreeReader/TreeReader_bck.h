//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 24 16:03:04 2014 by ROOT version 5.34/21
// from TTree SmallTree_WZ/SmallTree_WZ
// found on file: proof_WZ.root
//////////////////////////////////////////////////////////

#ifndef TreeReader_h
#define TreeReader_h

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

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           smalltree_nlepton;
   Float_t         smalltree_lept_pt[3];   //[smalltree_nlepton]
   Float_t         smalltree_lept_eta[3];   //[smalltree_nlepton]
   Float_t         smalltree_lept_phi[3];   //[smalltree_nlepton]
   Float_t         smalltree_lept_iso[3];   //[smalltree_nlepton]
   Int_t           smalltree_lept_flav[3];   //[smalltree_nlepton]
   Int_t           smalltree_njets;
   Float_t         smalltree_jet_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri[100];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri_up[100];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri_down[100];   //[smalltree_njets]
   Int_t           smalltree_jet_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jesup_njets;
   Float_t         smalltree_jet_jesup_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jesup_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jesdown_njets;
   Float_t         smalltree_jet_jesdown_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jesdown_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jerup_njets;
   Float_t         smalltree_jet_jerup_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jerup_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jerdown_njets;
   Float_t         smalltree_jet_jerdown_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jerdown_flav[100];   //[smalltree_njets]
   Float_t         smalltree_met_jesup_pt;
   Float_t         smalltree_met_jesup_phi;
   Float_t         smalltree_met_jesdown_pt;
   Float_t         smalltree_met_jesdown_phi;
   Float_t         smalltree_met_jerup_pt;
   Float_t         smalltree_met_jerup_phi;
   Float_t         smalltree_met_jerdown_pt;
   Float_t         smalltree_met_jerdown_phi;
   Float_t         smalltree_met_unclsup_pt;
   Float_t         smalltree_met_unclsup_phi;
   Float_t         smalltree_met_unclsdown_pt;
   Float_t         smalltree_met_unclsdown_phi;
   Float_t         smalltree_met_pt;
   Float_t         smalltree_met_phi;
   Float_t         smalltree_weight_trigup;
   Float_t         smalltree_weight_trigdown;

   // List of branches
   TBranch        *b_smalltree_nlepton;   //!
   TBranch        *b_smalltree_lept_pt;   //!
   TBranch        *b_smalltree_lept_eta;   //!
   TBranch        *b_smalltree_lept_phi;   //!
   TBranch        *b_smalltree_lept_iso;   //!
   TBranch        *b_smalltree_lept_flav;   //!
   TBranch        *b_smalltree_njets;   //!
   TBranch        *b_smalltree_jet_pt;   //!
   TBranch        *b_smalltree_jet_eta;   //!
   TBranch        *b_smalltree_jet_phi;   //!
   TBranch        *b_smalltree_jet_btagdiscri;   //!
   TBranch        *b_smalltree_jet_btagdiscri_up;   //!
   TBranch        *b_smalltree_jet_btagdiscri_down;   //!
   TBranch        *b_smalltree_jet_flav;   //!
   TBranch        *b_smalltree_jesup_njets;   //!
   TBranch        *b_smalltree_jet_jesup_pt;   //!
   TBranch        *b_smalltree_jet_jesup_eta;   //!
   TBranch        *b_smalltree_jet_jesup_phi;   //!
   TBranch        *b_smalltree_jet_jesup_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jesup_flav;   //!
   TBranch        *b_smalltree_jesdown_njets;   //!
   TBranch        *b_smalltree_jet_jesdown_pt;   //!
   TBranch        *b_smalltree_jet_jesdown_eta;   //!
   TBranch        *b_smalltree_jet_jesdown_phi;   //!
   TBranch        *b_smalltree_jet_jesdown_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jesdown_flav;   //!
   TBranch        *b_smalltree_jerup_njets;   //!
   TBranch        *b_smalltree_jet_jerup_pt;   //!
   TBranch        *b_smalltree_jet_jerup_eta;   //!
   TBranch        *b_smalltree_jet_jerup_phi;   //!
   TBranch        *b_smalltree_jet_jerup_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jerup_flav;   //!
   TBranch        *b_smalltree_jerdown_njets;   //!
   TBranch        *b_smalltree_jet_jerdown_pt;   //!
   TBranch        *b_smalltree_jet_jerdown_eta;   //!
   TBranch        *b_smalltree_jet_jerdown_phi;   //!
   TBranch        *b_smalltree_jet_jerdown_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jerdown_flav;   //!
   TBranch        *b_smalltree_met_jesup_pt;   //!
   TBranch        *b_smalltree_met_jesup_phi;   //!
   TBranch        *b_smalltree_met_jesdown_pt;   //!
   TBranch        *b_smalltree_met_jesdown_phi;   //!
   TBranch        *b_smalltree_met_jerup_pt;   //!
   TBranch        *b_smalltree_met_jerup_phi;   //!
   TBranch        *b_smalltree_met_jerdown_pt;   //!
   TBranch        *b_smalltree_met_jerdown_phi;   //!
   TBranch        *b_smalltree_met_unclsup_pt;   //!
   TBranch        *b_smalltree_met_unclsup_phi;   //!
   TBranch        *b_smalltree_met_unclsdown_pt;   //!
   TBranch        *b_smalltree_met_unclsdown_phi;   //!
   TBranch        *b_smalltree_met_pt;   //!
   TBranch        *b_smalltree_met_phi;   //!
   TBranch        *b_smalltree_weight_trigup;   //!
   TBranch        *b_smalltree_weight_trigdown;   //!

   TreeReader(TTree *tree=0);
   virtual ~TreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString sample);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   void initializeHisto(TString sample);
   void addHisto( TString var, TString selstep, TString sample, int nbins, float min, float max);
   void fillHisto(TString channel, TString var, TString selstep, TString sample, float val, float weight);
   
   void applyEventSel(TString channel, TString systtype, TString sample);
   
   
   TString determineChannel(int leptflav1, int leptflav2, int leptflav3);
   
   //void writeHisto();
   
   std::vector<TH1F*> histo_list_mmm;
   std::map<std::string, int> histo_map_mmm;
   
   std::vector<TH1F*> histo_list_mme;
   std::map<std::string, int> histo_map_mme;
   
   std::vector<TH1F*> histo_list_eem;
   std::map<std::string, int> histo_map_eem;
   
   std::vector<TH1F*> histo_list_eee;
   std::map<std::string, int> histo_map_eee;
   
   int numb_histo;
   
   
   
   
};

#endif

#ifdef TreeReader_cxx
TreeReader::TreeReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../RootFiles/proof_WZ.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../RootFiles/proof_WZ.root");
      }
      f->GetObject("SmallTree_WZ",tree);

   }
   Init(tree);
}

TreeReader::~TreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("smalltree_nlepton", &smalltree_nlepton, &b_smalltree_nlepton);
   fChain->SetBranchAddress("smalltree_lept_pt", smalltree_lept_pt, &b_smalltree_lept_pt);
   fChain->SetBranchAddress("smalltree_lept_eta", smalltree_lept_eta, &b_smalltree_lept_eta);
   fChain->SetBranchAddress("smalltree_lept_phi", smalltree_lept_phi, &b_smalltree_lept_phi);
   fChain->SetBranchAddress("smalltree_lept_iso", smalltree_lept_iso, &b_smalltree_lept_iso);
   fChain->SetBranchAddress("smalltree_lept_flav", smalltree_lept_flav, &b_smalltree_lept_flav);
   fChain->SetBranchAddress("smalltree_njets", &smalltree_njets, &b_smalltree_njets);
   fChain->SetBranchAddress("smalltree_jet_pt", smalltree_jet_pt, &b_smalltree_jet_pt);
   fChain->SetBranchAddress("smalltree_jet_eta", smalltree_jet_eta, &b_smalltree_jet_eta);
   fChain->SetBranchAddress("smalltree_jet_phi", smalltree_jet_phi, &b_smalltree_jet_phi);
   fChain->SetBranchAddress("smalltree_jet_btagdiscri", smalltree_jet_btagdiscri, &b_smalltree_jet_btagdiscri);
   fChain->SetBranchAddress("smalltree_jet_btagdiscri_up", smalltree_jet_btagdiscri_up, &b_smalltree_jet_btagdiscri_up);
   fChain->SetBranchAddress("smalltree_jet_btagdiscri_down", smalltree_jet_btagdiscri_down, &b_smalltree_jet_btagdiscri_down);
   fChain->SetBranchAddress("smalltree_jet_flav", smalltree_jet_flav, &b_smalltree_jet_flav);
   fChain->SetBranchAddress("smalltree_jesup_njets", &smalltree_jesup_njets, &b_smalltree_jesup_njets);
   fChain->SetBranchAddress("smalltree_jet_jesup_pt", smalltree_jet_jesup_pt, &b_smalltree_jet_jesup_pt);
   fChain->SetBranchAddress("smalltree_jet_jesup_eta", smalltree_jet_jesup_eta, &b_smalltree_jet_jesup_eta);
   fChain->SetBranchAddress("smalltree_jet_jesup_phi", smalltree_jet_jesup_phi, &b_smalltree_jet_jesup_phi);
   fChain->SetBranchAddress("smalltree_jet_jesup_btagdiscri", smalltree_jet_jesup_btagdiscri, &b_smalltree_jet_jesup_btagdiscri);
   fChain->SetBranchAddress("smalltree_jet_jesup_flav", smalltree_jet_jesup_flav, &b_smalltree_jet_jesup_flav);
   fChain->SetBranchAddress("smalltree_jesdown_njets", &smalltree_jesdown_njets, &b_smalltree_jesdown_njets);
   fChain->SetBranchAddress("smalltree_jet_jesdown_pt", smalltree_jet_jesdown_pt, &b_smalltree_jet_jesdown_pt);
   fChain->SetBranchAddress("smalltree_jet_jesdown_eta", smalltree_jet_jesdown_eta, &b_smalltree_jet_jesdown_eta);
   fChain->SetBranchAddress("smalltree_jet_jesdown_phi", smalltree_jet_jesdown_phi, &b_smalltree_jet_jesdown_phi);
   fChain->SetBranchAddress("smalltree_jet_jesdown_btagdiscri", smalltree_jet_jesdown_btagdiscri, &b_smalltree_jet_jesdown_btagdiscri);
   fChain->SetBranchAddress("smalltree_jet_jesdown_flav", smalltree_jet_jesdown_flav, &b_smalltree_jet_jesdown_flav);
   fChain->SetBranchAddress("smalltree_jerup_njets", &smalltree_jerup_njets, &b_smalltree_jerup_njets);
   fChain->SetBranchAddress("smalltree_jet_jerup_pt", smalltree_jet_jerup_pt, &b_smalltree_jet_jerup_pt);
   fChain->SetBranchAddress("smalltree_jet_jerup_eta", smalltree_jet_jerup_eta, &b_smalltree_jet_jerup_eta);
   fChain->SetBranchAddress("smalltree_jet_jerup_phi", smalltree_jet_jerup_phi, &b_smalltree_jet_jerup_phi);
   fChain->SetBranchAddress("smalltree_jet_jerup_btagdiscri", smalltree_jet_jerup_btagdiscri, &b_smalltree_jet_jerup_btagdiscri);
   fChain->SetBranchAddress("smalltree_jet_jerup_flav", smalltree_jet_jerup_flav, &b_smalltree_jet_jerup_flav);
   fChain->SetBranchAddress("smalltree_jerdown_njets", &smalltree_jerdown_njets, &b_smalltree_jerdown_njets);
   fChain->SetBranchAddress("smalltree_jet_jerdown_pt", smalltree_jet_jerdown_pt, &b_smalltree_jet_jerdown_pt);
   fChain->SetBranchAddress("smalltree_jet_jerdown_eta", smalltree_jet_jerdown_eta, &b_smalltree_jet_jerdown_eta);
   fChain->SetBranchAddress("smalltree_jet_jerdown_phi", smalltree_jet_jerdown_phi, &b_smalltree_jet_jerdown_phi);
   fChain->SetBranchAddress("smalltree_jet_jerdown_btagdiscri", smalltree_jet_jerdown_btagdiscri, &b_smalltree_jet_jerdown_btagdiscri);
   fChain->SetBranchAddress("smalltree_jet_jerdown_flav", smalltree_jet_jerdown_flav, &b_smalltree_jet_jerdown_flav);
   fChain->SetBranchAddress("smalltree_met_jesup_pt", &smalltree_met_jesup_pt, &b_smalltree_met_jesup_pt);
   fChain->SetBranchAddress("smalltree_met_jesup_phi", &smalltree_met_jesup_phi, &b_smalltree_met_jesup_phi);
   fChain->SetBranchAddress("smalltree_met_jesdown_pt", &smalltree_met_jesdown_pt, &b_smalltree_met_jesdown_pt);
   fChain->SetBranchAddress("smalltree_met_jesdown_phi", &smalltree_met_jesdown_phi, &b_smalltree_met_jesdown_phi);
   fChain->SetBranchAddress("smalltree_met_jerup_pt", &smalltree_met_jerup_pt, &b_smalltree_met_jerup_pt);
   fChain->SetBranchAddress("smalltree_met_jerup_phi", &smalltree_met_jerup_phi, &b_smalltree_met_jerup_phi);
   fChain->SetBranchAddress("smalltree_met_jerdown_pt", &smalltree_met_jerdown_pt, &b_smalltree_met_jerdown_pt);
   fChain->SetBranchAddress("smalltree_met_jerdown_phi", &smalltree_met_jerdown_phi, &b_smalltree_met_jerdown_phi);
   fChain->SetBranchAddress("smalltree_met_unclsup_pt", &smalltree_met_unclsup_pt, &b_smalltree_met_unclsup_pt);
   fChain->SetBranchAddress("smalltree_met_unclsup_phi", &smalltree_met_unclsup_phi, &b_smalltree_met_unclsup_phi);
   fChain->SetBranchAddress("smalltree_met_unclsdown_pt", &smalltree_met_unclsdown_pt, &b_smalltree_met_unclsdown_pt);
   fChain->SetBranchAddress("smalltree_met_unclsdown_phi", &smalltree_met_unclsdown_phi, &b_smalltree_met_unclsdown_phi);
   fChain->SetBranchAddress("smalltree_met_pt", &smalltree_met_pt, &b_smalltree_met_pt);
   fChain->SetBranchAddress("smalltree_met_phi", &smalltree_met_phi, &b_smalltree_met_phi);
   fChain->SetBranchAddress("smalltree_weight_trigup", &smalltree_weight_trigup, &b_smalltree_weight_trigup);
   fChain->SetBranchAddress("smalltree_weight_trigdown", &smalltree_weight_trigdown, &b_smalltree_weight_trigdown);
   Notify();
}

Bool_t TreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeReader_cxx
