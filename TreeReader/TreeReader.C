#define TreeReader_cxx
#include "TreeReader.h"

void TreeReader::Loop(TString sample)
{
   TString thesample  = sample;
   TFile * theoutputfile = new TFile( ("outputroot/histofile_"+thesample+".root").Data() , "recreate");
   
    
   for(unsigned int i=0; i< systlist.size(); i++){
     TString samplename = "";
     if( systlist[i]== "") samplename = thesample;
     else                  samplename = thesample+"__"+systlist[i];
     
     bool firstinit = false;
     if(i==0) firstinit = true;
     //cout << " iiii " << i << endl;
     //cout << samplename << endl;
     initializeHisto(samplename, firstinit);
    }

   isData = (sample== "DataEGZenriched" || sample== "DataMuEGZenriched" || sample== "DataMuZenriched" ||
	     sample== "DataEG" || sample== "DataMuEG" || sample== "DataMu");
   SetUpCSVreweighting();
   
    
   cout << "starting loops on events " << endl;
   
   
   
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(jentry%10000 == 0) cout << " processing " << sample << " event number " << jentry << endl;
      
      TString thechannel = determineChannel(smalltree_lept_flav[0], smalltree_lept_flav[1], smalltree_lept_flav[2]);
      
      
      //----------------------------------------------------------------------
      //apply event selection
      //second argument is for the systematics
      //"" means no sysetmatics
      //else, enter the systematic strings that appears in the histograms name
      //systematic names are :  
      //"leptup", "leptdown", "trigup", "trigdown", "PDFup", "PDFdown" 
      //"jesup", "jesdown", "jerup", "jerdown", "metunclsup", "metunclsdown"  
      //----------------------------------------------------------------------
      
      for(unsigned int isyst=0; isyst< systlist.size(); isyst++) applyEventSel(thechannel, systlist[isyst], thesample);
      
      
      
      
   }
   
   
   theoutputfile->Write();
   deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;
   
}



bool TreeReader::applyEventSel(TString thechannel, TString systtype, TString sample){
  
      bool applyCSV = 1;
      
      TLorentzVector leptZ1, leptZ2, leptW;
      
      leptZ1.SetPtEtaPhiM(smalltree_lept_pt[0], smalltree_lept_eta[0], smalltree_lept_phi[0], 0);
      leptZ2.SetPtEtaPhiM(smalltree_lept_pt[1], smalltree_lept_eta[1], smalltree_lept_phi[1], 0);
      leptW.SetPtEtaPhiM( smalltree_lept_pt[2], smalltree_lept_eta[2], smalltree_lept_phi[2], 0);
      
      
      double met_pt    = smalltree_met_pt;
      double met_phi   = smalltree_met_phi;
      double evtweight = -1;
      
      TString thesample = sample;
      if(systtype != "") thesample = thesample + "__"+systtype;
      //cout << thesample << thesample << endl;
      //cout << "smalltree_met_pt   " << smalltree_met_pt<< endl;
      int iter_jets          = 0;
      float * jet_pt         = 0;
      float * jet_eta        = 0;
      float * jet_phi        = 0;
      float * jet_btagdiscri = 0;
      int   * jet_flav       = 0;
      
      double wCSV = 1.;
      wCSV = GetCSVweight(0,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
   
      evtweight = smalltree_evtweight;
      iter_jets      = smalltree_njets;
      jet_pt	     = smalltree_jet_pt;
      jet_eta	     = smalltree_jet_eta;
      jet_phi	     = smalltree_jet_phi;
      jet_btagdiscri = smalltree_jet_btagdiscri;
      jet_flav       = smalltree_jet_flav;
      
      
      if(systtype == "" || 
         systtype == "lept__plus" ||   systtype == "lept__minus" || 
         systtype == "trig__plus" ||   systtype == "trig__minus" ||
         systtype == "PDF__plus"  ||   systtype == "PDF__minus"  
	 
      ){
	iter_jets      = smalltree_njets;
        jet_pt         = smalltree_jet_pt;
        jet_eta        = smalltree_jet_eta;
        jet_phi        = smalltree_jet_phi;
        jet_btagdiscri = smalltree_jet_btagdiscri;
        jet_flav       = smalltree_jet_flav;
	
	evtweight = smalltree_evtweight;
	//cout << "----------------" << endl;
	//if(     systtype == "lept__plus")  cout << "in lept up " << endl;
	//cout << "evtweight "  << evtweight << endl;
	if(     systtype == "lept__plus")    evtweight = smalltree_weight_leptup;
	else if(systtype == "lept__minus")   evtweight = smalltree_weight_leptdown;
	else if(systtype == "trig__plus")    evtweight = smalltree_weight_trigup;
	else if(systtype == "trig__minus")   evtweight = smalltree_weight_trigdown;
	else if(systtype == "PDF__plus")     evtweight = smalltree_weight_PDFup;
	else if(systtype == "PDF__minus")    evtweight = smalltree_weight_PDFdown;
	//cout << "evtweight "  << evtweight << endl;
	
        if(sample == "WZ" || sample == "WZHF" || sample == "tZq" || sample == "TTZ" || sample == "TTW" || sample == "ZZ"){
          if(thechannel == "mumumu" ) evtweight *= 0.9871;
          if(thechannel == "mumue"  ) evtweight *= 0.9001;
          if(thechannel == "eemu"   ) evtweight *= 0.9451;
          if(thechannel == "eee"    ) evtweight *= 0.9975;
        }
      
      }else if(systtype == "jes__plus"){
      
	met_pt    = smalltree_met_jesup_pt;
	met_phi   = smalltree_met_jesup_phi;
	iter_jets      = smalltree_jesup_njets;
        jet_pt         = smalltree_jet_jesup_pt;
        jet_eta        = smalltree_jet_jesup_eta;
        jet_phi        = smalltree_jet_jesup_phi;
        jet_btagdiscri = smalltree_jet_jesup_btagdiscri;
        jet_flav       = smalltree_jet_jesup_flav;
	
      }else if(systtype == "jes__minus"){
	met_pt    = smalltree_met_jesdown_pt;
	met_phi   = smalltree_met_jesdown_phi;
	iter_jets      = smalltree_jesdown_njets;
        jet_pt         = smalltree_jet_jesdown_pt;
        jet_eta        = smalltree_jet_jesdown_eta;
        jet_phi        = smalltree_jet_jesdown_phi;
        jet_btagdiscri = smalltree_jet_jesdown_btagdiscri;
        jet_flav       = smalltree_jet_jesdown_flav;
      }else if(systtype == "jer__plus"){
	met_pt    = smalltree_met_jerup_pt;
	met_phi   = smalltree_met_jerup_phi;
	iter_jets      = smalltree_jerup_njets;
        jet_pt         = smalltree_jet_jerup_pt;
        jet_eta        = smalltree_jet_jerup_eta;
        jet_phi        = smalltree_jet_jerup_phi;
        jet_btagdiscri = smalltree_jet_jerup_btagdiscri;
        jet_flav       = smalltree_jet_jerup_flav;
      }else if(systtype == "jer__minus"){
	met_pt    = smalltree_met_jerdown_pt;
	met_phi   = smalltree_met_jerdown_phi;
	iter_jets      = smalltree_jerdown_njets;
        jet_pt         = smalltree_jet_jerdown_pt;
        jet_eta        = smalltree_jet_jerdown_eta;
        jet_phi        = smalltree_jet_jerdown_phi;
        jet_btagdiscri = smalltree_jet_jerdown_btagdiscri;
        jet_flav       = smalltree_jet_jerdown_flav;
      }else if(systtype == "metuncls__plus"){
	met_pt    = smalltree_met_unclsup_pt;
	met_phi   = smalltree_met_unclsup_phi;
      }else if(systtype == "metuncls__minus"){
	met_pt    = smalltree_met_unclsdown_pt;
	met_phi   = smalltree_met_unclsdown_phi;
      }else if(systtype == "btag__JESup"){
	wCSV = GetCSVweight(1,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__JESdown"){
	wCSV = GetCSVweight(2,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVLFup"){
	wCSV = GetCSVweight(3,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVLFdown"){
	wCSV = GetCSVweight(4,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVHFStats1up"){
	wCSV = GetCSVweight(5,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVHFStats1down"){
	wCSV = GetCSVweight(6,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVHFStats2up"){
	wCSV = GetCSVweight(7,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVHFStats2down"){
	wCSV = GetCSVweight(8,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVCErr1up"){
	wCSV = GetCSVweight(9,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVCErr1down"){
	wCSV = GetCSVweight(10,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVCErr2up"){
	wCSV = GetCSVweight(11,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVCErr2down"){
	wCSV = GetCSVweight(12,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVHFup"){
	wCSV = GetCSVweight(13,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVHFdown"){
	wCSV = GetCSVweight(14,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVLFStats1up"){
	wCSV = GetCSVweight(15,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVLFStats1down"){
	wCSV = GetCSVweight(16,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVLFStats2up"){
	wCSV = GetCSVweight(17,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else if(systtype == "btag__CSVLFStats2down"){
	wCSV = GetCSVweight(18,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);
      }else{
        cout << "WARNING syst type not recognized !! " << endl;
	cout << "correct syst types are " << endl;
	cout << " \"\",  \"lept__plus\", \"lept__minus\", \"trig__plus\", \"trig__minus\", \"PDF__plus\", \"PDF__minus\"  " << endl;
	cout << "\"jes__plus\", \"jes__minus\", \"jer__plus\", \"jer__minus\", \"metuncls__plus\", \"metuncls__minus\" " << endl;      
      }
      
      if( applyCSV ) evtweight *= wCSV;
   
      TLorentzVector Zcand = leptZ1+leptZ2;
      
      
      double deltaR_LeptW_LeptZ1 = leptW.DeltaR(leptZ1);
      double deltaR_LeptW_LeptZ2 = leptW.DeltaR(leptZ2);
      
      
      //reconstruction of the W transverse mass
      //cout << "smalltree_lept_flav[2] " << smalltree_lept_flav[2]  << endl;
      
      if(abs(smalltree_lept_flav[2]) == 11  && smalltree_lept_iso[2] > isoEl &&  
      		 (sample != "DataEGZenriched" && sample!= "DataMuEGZenriched" && sample!= "DataMuZenriched" ) ) return false;
      if(abs(smalltree_lept_flav[2]) == 13  && smalltree_lept_iso[2] > isoMu &&  
      		 (sample!= "DataEGZenriched" && sample!= "DataMuEGZenriched" && sample!= "DataMuZenriched" ) ) return false;
      
      
      if(sample!= "DataEGZenriched" && sample!= "DataMuEGZenriched" && sample!= "DataMuZenriched"){
      
        if(abs(smalltree_lept_flav[0]) == 11  && smalltree_lept_iso[0] > 0.1) return false;
        if(abs(smalltree_lept_flav[1]) == 11  && smalltree_lept_iso[1] > 0.1) return false;
        if(abs(smalltree_lept_flav[2]) == 11  && smalltree_lept_iso[2] > 0.1) return false;
      }
      
      
      //reconstruction of the W transverse mass
      
      double mTW = pow(
         		2*leptW.Pt()*met_pt*(1-cos(leptW.Phi() -  met_phi))
         		,0.5);
      //cout << "mTW " <<  mTW<< endl;
      //reconstruction of the Z invarian mass
      double InvMass_ll = Zcand.M();
     
     //cout << "thesample " << thesample << endl;
     
     bool contain_HF = false;
     for(int ijet=0; ijet<iter_jets; ijet++){
       if( abs(jet_flav[ijet]) == 5 || abs(jet_flav[ijet]) == 4){
         contain_HF = true;
         //cout << jet_flav[ijet] << endl;
       } 
       
     }
     //if(contain_HF) cout << "contain_HF " << contain_HF << endl;
     //if( sample == "WZHF" ) return false;
     
     if(!contain_HF && sample == "WZHF" ) return false;
     //if( contain_HF && sample == "WZ"   ) return false;
     
     
     int njets=0;
     int nbjets = 0;
     std::vector<int > btagged_jet_idx;
     for(int ijet=0; ijet<iter_jets; ijet++){
       if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.5) continue;     
       njets++;
       fillHisto(thechannel, "JetPt",     "afterleptsel",  thesample,  jet_pt[ijet] , evtweight);
       fillHisto(thechannel, "JetEta",    "afterleptsel",  thesample,  jet_eta[ijet] , evtweight);
       if(jet_btagdiscri[ijet] > 0.244)  {nbjets++;   btagged_jet_idx.push_back(ijet); }  
     }
     
     
     fillHisto(thechannel, "NJet",       "afterleptsel",  thesample,  iter_jets , evtweight);
     fillHisto(thechannel, "NBJet",      "afterleptsel",  thesample,   nbjets, evtweight);
     
     
     fillHisto(thechannel, "mWT",        "afterleptsel",  thesample,   mTW,          evtweight);
     fillHisto(thechannel, "InvM_ll",    "afterleptsel",  thesample,   InvMass_ll,   evtweight);
     fillHisto(thechannel, "LeptPt",     "afterleptsel",  thesample,   leptZ1.Pt(),  evtweight);
     fillHisto(thechannel, "LeptEta",    "afterleptsel",  thesample,   leptZ1.Eta(), evtweight);
     fillHisto(thechannel, "LeptPt",     "afterleptsel",  thesample,   leptZ2.Pt(),  evtweight);
     fillHisto(thechannel, "LeptEta",    "afterleptsel",  thesample,   leptZ2.Eta(), evtweight);
     fillHisto(thechannel, "LeptPt",     "afterleptsel",  thesample,   leptW.Pt(),   evtweight);
     fillHisto(thechannel, "LeptEta",    "afterleptsel",  thesample,   leptW.Eta(),  evtweight);
     fillHisto(thechannel, "LeptPtZ1",	 "afterleptsel",  thesample,  leptZ1.Pt(),  evtweight);
     fillHisto(thechannel, "LeptEtaZ1",  "afterleptsel",  thesample,  leptZ1.Eta(), evtweight);
     fillHisto(thechannel, "LeptPtZ2",	 "afterleptsel",  thesample,  leptZ2.Pt(),  evtweight);
     fillHisto(thechannel, "LeptEtaZ2",  "afterleptsel",  thesample,  leptZ2.Eta(), evtweight);
     fillHisto(thechannel, "LeptPtW",	 "afterleptsel",  thesample,  leptW.Pt(),   evtweight);
     fillHisto(thechannel, "LeptEtaW",	 "afterleptsel",  thesample,  leptW.Eta(),  evtweight);
     
     
     fillHisto(thechannel, "CutFlow", "",  thesample, 0, evtweight);
      
      
      
      //---------------------------------
      //apply dilepton invariant mass cut
       if( fabs(InvMass_ll-91) < 15 && deltaR_LeptW_LeptZ1 > 0.1 && deltaR_LeptW_LeptZ2 > 0.1 ){
       
         for(int ijet=0; ijet<iter_jets; ijet++){
           if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.5) continue;     
           fillHisto(thechannel, "JetPt",     "afterZsel",  thesample,  jet_pt[ijet] , evtweight);
           fillHisto(thechannel, "JetEta",    "afterZsel",  thesample,  jet_eta[ijet] , evtweight);
         }
     
     
         fillHisto(thechannel, "NJet",       "afterZsel",  thesample,  iter_jets , evtweight);
         fillHisto(thechannel, "NBJet",      "afterZsel",  thesample,   nbjets, evtweight);
     
     
         fillHisto(thechannel, "mWT",        "afterZsel",  thesample,   mTW,          evtweight);
         fillHisto(thechannel, "InvM_ll",    "afterZsel",  thesample,   InvMass_ll,   evtweight);
         fillHisto(thechannel, "LeptPt",     "afterZsel",  thesample,   leptZ1.Pt(),  evtweight);
         fillHisto(thechannel, "LeptEta",    "afterZsel",  thesample,   leptZ1.Eta(), evtweight);
         fillHisto(thechannel, "LeptPt",     "afterZsel",  thesample,   leptZ2.Pt(),  evtweight);
         fillHisto(thechannel, "LeptEta",    "afterZsel",  thesample,   leptZ2.Eta(), evtweight);
         fillHisto(thechannel, "LeptPt",     "afterZsel",  thesample,   leptW.Pt(),   evtweight);
         fillHisto(thechannel, "LeptEta",    "afterZsel",  thesample,   leptW.Eta(),  evtweight);
         fillHisto(thechannel, "LeptPtZ1",   "afterZsel",  thesample,   leptZ1.Pt(),  evtweight);
         fillHisto(thechannel, "LeptEtaZ1",  "afterZsel",  thesample,   leptZ1.Eta(), evtweight);
         fillHisto(thechannel, "LeptPtZ2",   "afterZsel",  thesample,   leptZ2.Pt(),  evtweight);
         fillHisto(thechannel, "LeptEtaZ2",  "afterZsel",  thesample,   leptZ2.Eta(), evtweight);
         fillHisto(thechannel, "LeptPtW",    "afterZsel",  thesample,   leptW.Pt(),   evtweight);
         fillHisto(thechannel, "LeptEtaW",   "afterZsel",  thesample,   leptW.Eta(),  evtweight);
     
     
  
         fillHisto(thechannel, "CutFlow", "",  thesample, 1, evtweight);
	
	
	 //------------------------
	 //ask for at least one jet
	 if(njets>=1){

	   //store treeCSV channel
//	   if(thechannel == "mumumu") treeCSV_Channel = 0;
//	   if(thechannel == "mumue" ) treeCSV_Channel = 1;
//	   if(thechannel == "eemu"  ) treeCSV_Channel = 2;
//	   if(thechannel == "eee"   ) treeCSV_Channel = 3;
	   
//	   treeCSV_EvtWeight  =  evtweight;
//	   treeCSV_jetN    = iter_jets;
	   
//	   for(int ijet=0; ijet<iter_jets; ijet++)
//	      {
//		 treeCSV_jetPt[ijet]    = jet_pt[ijet];
//		 treeCSV_jetEta[ijet]    = jet_eta[ijet];
//		 treeCSV_jetCSV[ijet]    = jet_btagdiscri[ijet];
//		 treeCSV_jetFlav[ijet]    = jet_flav[ijet];
//	      }	    
	    
//	   if(theTreeCSV_map[thesample] != 0)  theTreeCSV_map[thesample]->Fill();
//	   else cout << "wrong sample name given to CSV TTree " << thesample <<endl; 
	    
           for(int ijet=0; ijet<iter_jets; ijet++){
             if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.5) continue;     
             fillHisto(thechannel, "JetPt",     "afterjetsel",  thesample,  jet_pt[ijet] , evtweight);
             fillHisto(thechannel, "JetEta",    "afterjetsel",  thesample,  jet_eta[ijet] , evtweight);
           }
     
     
           fillHisto(thechannel, "NJet",       "afterjetsel",  thesample,  iter_jets , evtweight);
           fillHisto(thechannel, "NBJet",      "afterjetsel",  thesample,   nbjets, evtweight);
     
     
           fillHisto(thechannel, "mWT",        "afterjetsel",  thesample,   mTW,          evtweight);
           fillHisto(thechannel, "InvM_ll",    "afterjetsel",  thesample,   InvMass_ll,   evtweight);
           fillHisto(thechannel, "LeptPt",     "afterjetsel",  thesample,   leptZ1.Pt(),  evtweight);
           fillHisto(thechannel, "LeptEta",    "afterjetsel",  thesample,   leptZ1.Eta(), evtweight);
           fillHisto(thechannel, "LeptPt",     "afterjetsel",  thesample,   leptZ2.Pt(),  evtweight);
           fillHisto(thechannel, "LeptEta",    "afterjetsel",  thesample,   leptZ2.Eta(), evtweight);
           fillHisto(thechannel, "LeptPt",     "afterjetsel",  thesample,   leptW.Pt(),   evtweight);
           fillHisto(thechannel, "LeptEta",    "afterjetsel",  thesample,   leptW.Eta(),  evtweight);
           fillHisto(thechannel, "LeptPtZ1",   "afterjetsel",  thesample,   leptZ1.Pt(),  evtweight);
           fillHisto(thechannel, "LeptEtaZ1",  "afterjetsel",  thesample,   leptZ1.Eta(), evtweight);
           fillHisto(thechannel, "LeptPtZ2",   "afterjetsel",  thesample,   leptZ2.Pt(),  evtweight);
           fillHisto(thechannel, "LeptEtaZ2",  "afterjetsel",  thesample,   leptZ2.Eta(), evtweight);
           fillHisto(thechannel, "LeptPtW",    "afterjetsel",  thesample,   leptW.Pt(),   evtweight);
           fillHisto(thechannel, "LeptEtaW",   "afterjetsel",  thesample,   leptW.Eta(),  evtweight);
     
     
	  
           fillHisto(thechannel, "CutFlow", "",  thesample, 2, evtweight);
	   
	   
	   
	   
	   
	   
	   //----------------------------
	   //no more than one btagged jet 
	   if(nbjets <=2 && njets<=2 && mTW> 40 ){
	  
	     float HT = 0;
	     float ST = 0;
	     
             for(int ijet=0; ijet<iter_jets; ijet++){
               if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.5) continue;     
               fillHisto(thechannel, "JetPt",     "afterbjetsel",  thesample,  jet_pt[ijet] , evtweight);
               fillHisto(thechannel, "JetEta",    "afterbjetsel",  thesample,  jet_eta[ijet] , evtweight);
	       HT+= jet_pt[ijet];
	       ST+= jet_pt[ijet];
             }
             
	     ST+= leptZ1.Pt()+leptZ1.Pt()+leptW.Pt();
	     
	     
	     
	     
	     TLorentzVector bjet, themet;
	     if(nbjets >  0) bjet.SetPtEtaPhiM(jet_pt[btagged_jet_idx[0]], jet_eta[btagged_jet_idx[0]], jet_phi[btagged_jet_idx[0]], 0);
	     else            bjet.SetPtEtaPhiM(jet_pt[0],                  jet_eta[0],                  jet_phi[0],                  0);
	     
	     themet.SetPtEtaPhiM(met_pt, 0, met_phi, 0);
	     
	     
	     
	      tree_EvtWeight  =  evtweight;
	     	
             // Top mass computation
	     //leptW + selJets[0].p4 + met.p4; 

	     double term1 = leptW.Pz()*(leptW.Px()*themet.Px() + leptW.Py()*themet.Py() + (80.399)*(80.399)/2.);
	     double det = leptW.Px()*themet.Px() + leptW.Py()*themet.Py() + (80.399)*(80.399)/2.
	        	    - themet.Pt()*themet.Pt()*(leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz());
	     if(det<0) det=0;
	     double term2 = leptW.E()*pow(det, 0.5);
	     double denom = leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz();

             
	     double sol1 = (term1 - term2)/denom;
	     //double sol2 = (term1 + term2)/denom;
  
	     double neutrE = pow( pow(themet.Px(),2) + pow(themet.Py(),2) + sol1*sol1, 0.5);//neglecting neut mass 

	     TLorentzVector neutrino;
	     neutrino.SetPxPyPzE( themet.Px(), themet.Py(), sol1, neutrE);
	  
	     TLorentzVector theWcand = neutrino + leptW;

	     TLorentzVector topCand = neutrino + leptW + bjet ;
	  

	  
             tree_topMass    = topCand.M();
             tree_totMass    = (topCand + (leptZ1+leptZ2)).M();
             tree_deltaPhilb = fabs(leptW.DeltaPhi(bjet));
             tree_deltaRlb   = leptW.DeltaR(  bjet);
             tree_deltaRTopZ = (leptZ1+leptZ2).DeltaR(topCand);
	     
	     int  leptW_Charge = 1;
	     if(smalltree_lept_flav[2] < 0 ) leptW_Charge = -1;
             tree_asym	  = leptW_Charge*fabs(leptW.Eta());
             tree_Zpt	  = (leptZ1+leptZ2).Pt();
             tree_ZEta	  = (leptZ1+leptZ2).Eta();
             tree_topPt	  = topCand.Pt();
             tree_topEta	  = topCand.Eta();
	  
	     //if(topCand.M() < 0) cout << "topCand.M() " << tree_topMass << endl;
	     //if(topCand.M() > 100000) cout << "topCand.M() " <<  tree_topMass<< endl;
	    
	     tree_totPt	    = (topCand + (leptZ1+leptZ2)).Pt();
	     tree_totEta	    = (topCand + (leptZ1+leptZ2)).Eta(); 

  	     tree_deltaRZl     = (leptZ1+leptZ2).DeltaR(leptW);
  	     tree_deltaPhiZmet = (leptZ1+leptZ2).DeltaPhi(themet);
  	     if(nbjets >  0) tree_btagDiscri   = jet_btagdiscri[btagged_jet_idx[0]];
	     else tree_btagDiscri   = jet_btagdiscri[0];

	     
	     
	     
  	     tree_NJets	    = float(njets);
	  
	
 
	     tree_NBJets	    = float(nbjets);
	  
  	     tree_leptWPt        = leptW.Pt(); 
	     tree_leptWEta       = leptW.Eta();
  	     tree_leadJetPt      = bjet.Pt(); 
	     tree_leadJetEta     = bjet.Eta();
             tree_deltaRZleptW   = (leptZ1+leptZ2).DeltaR(leptW); 
	     tree_deltaPhiZleptW = (leptZ1+leptZ2).DeltaPhi(leptW);
	  
	     tree_met = themet.Pt();
	     tree_mTW = mTW;
	  
	  
	     double Zpt = (leptZ1+leptZ2).Pt();
	     double deltaPhilj = fabs(leptW.DeltaPhi(bjet));
	     
	     TLorentzVector theWcand_topRF = theWcand;
	     TLorentzVector letp_WRF = leptW;  
	     theWcand_topRF.Boost( TVector3(-1.0*topCand.BoostVector().Px(), -1.0*topCand.BoostVector().Py() ,-1.0*topCand.BoostVector().Pz()));
	     letp_WRF.Boost(       TVector3(-1.0*theWcand.BoostVector().Px(),-1.0*theWcand.BoostVector().Py(),-1.0*theWcand.BoostVector().Pz()));
	  
             double cosThetaStar = cos(letp_WRF.Vect().Angle(theWcand_topRF.Vect()));
            
             fillHisto(thechannel, "CutFlow", "",  thesample, 3, evtweight);
	     
	     //cut on top mass
	     if(tree_topMass < 260 && tree_totPt> 30){ 
	       fillHisto(thechannel, "HT",       "afterbjetsel",  thesample,  HT , evtweight);
	       fillHisto(thechannel, "ST",       "afterbjetsel",  thesample,  ST,  evtweight);
	     
             
               fillHisto(thechannel, "NJet",       "afterbjetsel",  thesample,  iter_jets , evtweight);
               fillHisto(thechannel, "NBJet",      "afterbjetsel",  thesample,   nbjets, evtweight);
     
     
               fillHisto(thechannel, "mWT",        "afterbjetsel",  thesample,   mTW,          evtweight);
               fillHisto(thechannel, "InvM_ll",    "afterbjetsel",  thesample,   InvMass_ll,   evtweight);
               fillHisto(thechannel, "LeptPt",     "afterbjetsel",  thesample,   leptZ1.Pt(),  evtweight);
               fillHisto(thechannel, "LeptEta",    "afterbjetsel",  thesample,   leptZ1.Eta(), evtweight);
               fillHisto(thechannel, "LeptPt",     "afterbjetsel",  thesample,   leptZ2.Pt(),  evtweight);
               fillHisto(thechannel, "LeptEta",    "afterbjetsel",  thesample,   leptZ2.Eta(), evtweight);
               fillHisto(thechannel, "LeptPt",     "afterbjetsel",  thesample,   leptW.Pt(),   evtweight);
               fillHisto(thechannel, "LeptEta",    "afterbjetsel",  thesample,   leptW.Eta(),  evtweight);
               fillHisto(thechannel, "LeptPtZ1",   "afterbjetsel",  thesample,   leptZ1.Pt(),  evtweight);
               fillHisto(thechannel, "LeptEtaZ1",  "afterbjetsel",  thesample,   leptZ1.Eta(), evtweight);
               fillHisto(thechannel, "LeptPtZ2",   "afterbjetsel",  thesample,   leptZ2.Pt(),  evtweight);
               fillHisto(thechannel, "LeptEtaZ2",  "afterbjetsel",  thesample,   leptZ2.Eta(), evtweight);
               fillHisto(thechannel, "LeptPtW",    "afterbjetsel",  thesample,   leptW.Pt(),   evtweight);
               fillHisto(thechannel, "LeptEtaW",   "afterbjetsel",  thesample,   leptW.Eta(),  evtweight);
     
	       fillHisto(thechannel, "BJetCSV",   "afterbjetsel",  thesample,   tree_btagDiscri,  evtweight);
	         
               fillHisto(thechannel, "CutFlow", "",  thesample, 4, evtweight);
               fillHisto(thechannel, "toppT",     "afterbjetsel",  thesample, tree_topPt,  evtweight);
               fillHisto(thechannel, "ZpT",       "afterbjetsel",  thesample, 	   Zpt,  evtweight);
               fillHisto(thechannel, "WpT",       "afterbjetsel",  thesample, tree_leptWPt   ,  evtweight);
               fillHisto(thechannel, "MET",       "afterbjetsel",  thesample,  tree_met ,  evtweight);
               fillHisto(thechannel, "totpT",     "afterbjetsel",  thesample,  tree_totPt, evtweight);
	       fillHisto(thechannel, "DeltaPhiTopZ",     "afterbjetsel",  thesample, (leptZ1+leptZ2).DeltaPhi(topCand) , evtweight);
	     
               fillHisto(thechannel, "topMass",     "afterbjetsel",  thesample, tree_topMass , evtweight);
               fillHisto(thechannel, "ZMass",       "afterbjetsel",  thesample,  InvMass_ll, evtweight);
               fillHisto(thechannel, "tZMass",      "afterbjetsel",  thesample, tree_totMass , evtweight);
	     
	        tree_cosThetaStar = cosThetaStar;
	      
	        //store tree channel
	        if(thechannel == "mumumu") tree_Channel = 0;
	        if(thechannel == "mumue" ) tree_Channel = 1;
	        if(thechannel == "eemu"  ) tree_Channel = 2;
	        if(thechannel == "eee"   ) tree_Channel = 3;
	      
	      
	      
	        if(theTree_map[thesample] != 0)  theTree_map[thesample]->Fill();
	        else cout << "wrong sample name given to TTree " << thesample <<endl;
	      
	      //cout << thesample << endl;
	    }// end top mass selection
	  }//end btag selection	  
	}//end jet selection
      }//end dilepton M cut
      
}




//------------------------------------------------------
//initialize the historams for the analysis
//------------------------------------------------------


void TreeReader::initializeHisto(TString sample, bool isfirstset){
  
  
  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  cout << " initialize histograms               " << endl;
  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  
  
  if(isfirstset){
    numb_histo = 0;
    TH1F * first_emptyHisto = new TH1F("first_emptyHisto", "first_emptyHisto", 100, 0, 1000);
    histo_list_mmm.push_back(first_emptyHisto);
    histo_list_mme.push_back(first_emptyHisto);
    histo_list_eem.push_back(first_emptyHisto);
    histo_list_eee.push_back(first_emptyHisto);
  
    numb_histo++;
  }
  cout << numb_histo << endl;
  addHisto("CutFlow", "", sample.Data(),  15,-0.5,14.5);
  addHisto("NVtx",    "", sample.Data(),  60, 0, 60); 
  
  
  
  //after lepton selection
  addHisto("NJet",      "afterleptsel",  sample.Data(),  5,-0.5,4.5);
  addHisto("NBJet",     "afterleptsel",  sample.Data(),   5,-0.5,4.5);
  addHisto("mWT",       "afterleptsel",  sample.Data(),   100,0,200);
  addHisto("InvM_ll",   "afterleptsel",  sample.Data(),   100,0,200);
  addHisto("JetPt",     "afterleptsel",  sample.Data(),   100,0,300) ;
  addHisto("JetEta",    "afterleptsel",  sample.Data(),   26, -2.5, 2.5 ) ;
  addHisto("LeptPt",    "afterleptsel",  sample.Data(),   100,0.,200);
  addHisto("LeptEta",   "afterleptsel",  sample.Data(),   26, -2.5, 2.5 );
  
  addHisto("LeptPtZ1",    "afterleptsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaZ1",   "afterleptsel",  sample.Data(),  26, -2.5, 2.5 );
  addHisto("LeptPtZ2",    "afterleptsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaZ2",   "afterleptsel",  sample.Data(),  26, -2.5, 2.5 );
  addHisto("LeptPtW",     "afterleptsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaW",    "afterleptsel",  sample.Data(),  26, -2.5, 2.5 );
  
  
  
  
  //after Z selection selection
  addHisto("NJet",      "afterZsel",  sample.Data(),  5,-0.5,4.5);
  addHisto("NBJet",     "afterZsel",  sample.Data(),   5,-0.5,4.5);
  addHisto("mWT",       "afterZsel",  sample.Data(),   100,0,200);
  addHisto("InvM_ll",   "afterZsel",  sample.Data(),   100,0,200);
  addHisto("JetPt",     "afterZsel",  sample.Data(),   100,0,300) ;
  addHisto("JetEta",    "afterZsel",  sample.Data(),   26, -2.5, 2.5 ) ;
  addHisto("LeptPt",    "afterZsel",  sample.Data(),   100,0.,200);
  addHisto("LeptEta",   "afterZsel",  sample.Data(),   26, -2.5, 2.5 );
  
  addHisto("LeptPtZ1",    "afterZsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaZ1",   "afterZsel",  sample.Data(),  26, -2.5, 2.5 );
  addHisto("LeptPtZ2",    "afterZsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaZ2",   "afterZsel",  sample.Data(),  26, -2.5, 2.5 );
  addHisto("LeptPtW",     "afterZsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaW",    "afterZsel",  sample.Data(),  26, -2.5, 2.5 );
  
  
   
  //after jet selection
  addHisto("NJet",      "afterjetsel",  sample.Data(),  5,-0.5,4.5);
  addHisto("NBJet",     "afterjetsel",  sample.Data(),   5,-0.5,4.5);
  addHisto("mWT",       "afterjetsel",  sample.Data(),   100,0,200);
  addHisto("InvM_ll",   "afterjetsel",  sample.Data(),   100,0,200);
  addHisto("JetPt",     "afterjetsel",  sample.Data(),   100,0,300) ;
  addHisto("JetEta",    "afterjetsel",  sample.Data(),   26, -2.5, 2.5 ) ;
  addHisto("LeptPt",    "afterjetsel",  sample.Data(),   100,0.,200);
  addHisto("LeptEta",   "afterjetsel",  sample.Data(),   26, -2.5, 2.5 );
  
  addHisto("LeptPtZ1",    "afterjetsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaZ1",   "afterjetsel",  sample.Data(),  26, -2.5, 2.5 );
  addHisto("LeptPtZ2",    "afterjetsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaZ2",   "afterjetsel",  sample.Data(),  26, -2.5, 2.5 );
  addHisto("LeptPtW",     "afterjetsel",  sample.Data(),  100,0.,200);
  addHisto("LeptEtaW",    "afterjetsel",  sample.Data(),  26, -2.5, 2.5 );
    
  //after b-jet selection
  addHisto("NJet",      "afterbjetsel",  sample.Data(),  5,-0.5,4.5);
  addHisto("NBJet",     "afterbjetsel",  sample.Data(),   5,-0.5,4.5);
  addHisto("mWT",       "afterbjetsel",  sample.Data(),   100,0,200);
  addHisto("InvM_ll",   "afterbjetsel",  sample.Data(),   100,0,200);
  addHisto("JetPt",     "afterbjetsel",  sample.Data(),   100,0,300) ;
  addHisto("JetEta",    "afterbjetsel",  sample.Data(),   26, -2.5, 2.5 ) ;
  addHisto("LeptPt",    "afterbjetsel",  sample.Data(),   100,0.,200);
  addHisto("LeptEta",   "afterbjetsel",  sample.Data(),   26, -2.5, 2.5 );
  
  addHisto("LeptPtZ1",    "afterbjetsel",  sample.Data(),  25,0.,200);
  addHisto("LeptEtaZ1",   "afterbjetsel",  sample.Data(),  13, -2.5, 2.5 );
  addHisto("LeptPtZ2",    "afterbjetsel",  sample.Data(),  25,0.,200);
  addHisto("LeptEtaZ2",   "afterbjetsel",  sample.Data(),  13, -2.5, 2.5 );
  addHisto("LeptPtW",     "afterbjetsel",  sample.Data(),  25,0.,200);
  addHisto("LeptEtaW",    "afterbjetsel",  sample.Data(),  13, -2.5, 2.5 );
  
  addHisto("HT",    "afterbjetsel",  sample.Data(),  50,0.,1000);
  addHisto("ST",    "afterbjetsel",  sample.Data(),  50,0.,1000);
  
  addHisto("toppT",    "afterbjetsel",  sample.Data(),  50,0.,500);
  addHisto("ZpT",      "afterbjetsel",  sample.Data(),  50,0.,500);
  addHisto("WpT",      "afterbjetsel",  sample.Data(),  50,0.,500);
  addHisto("MET",      "afterbjetsel",  sample.Data(),  50,0.,500);
  addHisto("totpT",    "afterbjetsel",  sample.Data(),  50,0.,500);
  addHisto("DeltaPhiTopZ",    "afterbjetsel",  sample.Data(),  50,0.,3.2);
  addHisto("topMass",         "afterbjetsel",  sample.Data(),  50,0.,700);
  addHisto("ZMass",           "afterbjetsel",  sample.Data(),  50,0.,500);
  addHisto("tZMass",          "afterbjetsel",  sample.Data(),  70, 150,1000);
  
  
  
  
  
  addHisto("BJetCSV",    "afterbjetsel",  sample.Data(),  30,0.,1.);
   
  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  cout << " histograms  initialized             " << endl;
  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  
    
  //--------------------------------------//
  //   Output TTree 	
  //--------------------------------------//
  TString treename = "Ttree_"+sample;
  cout << treename << endl;
  TTree * TheTree = new TTree(treename.Data(),treename.Data());
  TheTree->Branch("tree_cosThetaStar",     &tree_cosThetaStar,     "tree_cosThetaStar/F"   );
  TheTree->Branch("tree_topMass",     &tree_topMass,     "tree_topMass/F"   );
  TheTree->Branch("tree_totMass",     &tree_totMass,     "tree_totMass/F"   );
  TheTree->Branch("tree_deltaPhilb",  &tree_deltaPhilb,  "tree_deltaPhilb/F");
  TheTree->Branch("tree_deltaRlb",    &tree_deltaRlb,    "tree_deltaRlb/F"  );
  TheTree->Branch("tree_deltaRTopZ",  &tree_deltaRTopZ,  "tree_deltaRTopZ/F");
  TheTree->Branch("tree_asym",        &tree_asym,        "tree_asym/F"      );
  TheTree->Branch("tree_Zpt",         &tree_Zpt,         "tree_Zpt/F"       );
  TheTree->Branch("tree_ZEta",        &tree_ZEta,        "tree_ZEta/F"      );
  TheTree->Branch("tree_topPt",       &tree_topPt,       "tree_topPt/F"     );
  TheTree->Branch("tree_topEta",      &tree_topEta,      "tree_topEta/F"    );
  TheTree->Branch("tree_NJets",       &tree_NJets,       "tree_NJets/F"     );
  TheTree->Branch("tree_NBJets",      &tree_NBJets,      "tree_NBJets/F"    );
  TheTree->Branch("tree_deltaRZl",    &tree_deltaRZl,    "tree_deltaRZl/F"     );
  TheTree->Branch("tree_deltaPhiZmet",&tree_deltaPhiZmet,"tree_deltaPhiZmet/F" );
  TheTree->Branch("tree_btagDiscri",  &tree_btagDiscri,  "tree_btagDiscri/F"   );
  
  TheTree->Branch("tree_totPt",      &tree_totPt,      "tree_totPt/F"   );
  TheTree->Branch("tree_totEta",     &tree_totEta,     "tree_totEta/F"   );
  
  
  TheTree->Branch("tree_leptWPt",        &tree_leptWPt        , "tree_leptWPt/F"        );
  TheTree->Branch("tree_leptWEta",       &tree_leptWEta       , "tree_leptWEta/F"       );
  TheTree->Branch("tree_leadJetPt",      &tree_leadJetPt      , "tree_leadJetPt/F"      ); 
  TheTree->Branch("tree_leadJetEta",     &tree_leadJetEta     , "tree_leadJetEta/F"     );
  TheTree->Branch("tree_deltaRZleptW",   &tree_deltaRZleptW   , "tree_deltaRZleptW/F"   );
  TheTree->Branch("tree_deltaPhiZleptW", &tree_deltaPhiZleptW , "tree_deltaPhiZleptW/F" );
  
  
  TheTree->Branch("tree_met", &tree_met , "tree_met/F" );
  TheTree->Branch("tree_mTW", &tree_mTW , "tree_mTW/F" );
  
  
  TheTree->Branch("tree_EvtWeight",   &tree_EvtWeight,   "tree_EvtWeight/F" );
  TheTree->Branch("tree_SampleType",  &tree_SampleType,  "tree_SampleType/I");
  TheTree->Branch("tree_Channel",     &tree_Channel,     "tree_Channel/I"   );

  theTree_list.push_back(TheTree);
  theTree_map[sample] = theTree_list.back();
  
      
  tree_cosThetaStar   = -10000;
  tree_EvtWeight      = -10000;    
  tree_topMass        = -10000;
  tree_totMass        = -10000; 
  tree_deltaPhilb     = -10000; 
  tree_deltaRlb       = -10000; 
  tree_deltaRTopZ     = -10000; 
  tree_asym	      = -10000; 
  tree_Zpt	      = -10000; 
  tree_ZEta	      = -10000; 
  tree_topPt	      = -10000; 
  tree_topEta	      = -10000; 
  tree_totPt	      = -10000;
  tree_totEta	      = -10000; 
  tree_deltaRZl       = -10000; 
  tree_deltaPhiZmet   = -10000; 
  tree_btagDiscri     = -10000; 
  tree_NJets	      = -10000; 
  tree_NBJets	      = -10000;
  tree_leptWPt        = -10000; 
  tree_leptWEta       = -10000;
  tree_leadJetPt      = -10000; 
  tree_leadJetEta     = -10000; 
  tree_deltaRZleptW   = -10000; 
  tree_deltaPhiZleptW = -10000; 
  tree_met            = -10000;
  tree_mTW            = -10000;
  tree_Channel        = -10000; 
  tree_SampleType     = -10000; 
  

  //--------------------------------------//
  //   Output CSV TTree
  //--------------------------------------//
   
//  TString treenameCSV = "TtreeCSV_"+sample;
//  cout << treenameCSV << endl;
//  TTree * TheTreeCSV = new TTree(treenameCSV.Data(),treenameCSV.Data());
//  TheTreeCSV->Branch("treeCSV_jetN",        &treeCSV_jetN,        "treeCSV_jetN/I"      );
//  TheTreeCSV->Branch("treeCSV_jetPt",       &treeCSV_jetPt,       "treeCSV_jetPt/F"     );
//  TheTreeCSV->Branch("treeCSV_jetEta",      &treeCSV_jetEta,      "treeCSV_jetEta/F"    );
//  TheTreeCSV->Branch("treeCSV_jetCSV",      &treeCSV_jetCSV,      "treeCSV_jetCSV/F"    );
//  TheTreeCSV->Branch("treeCSV_jetFlav",     &treeCSV_jetFlav,     "treeCSV_jetFlav/I"   );
     
//  TheTreeCSV->Branch("treeCSV_EvtWeight",   &treeCSV_EvtWeight,   "treeCSV_EvtWeight/F" );
//  TheTreeCSV->Branch("treeCSV_SampleType",  &treeCSV_SampleType,  "treeCSV_SampleType/I");
//  TheTreeCSV->Branch("treeCSV_Channel",     &treeCSV_Channel,     "treeCSV_Channel/I"   );
  
//  theTreeCSV_list.push_back(TheTreeCSV);
//  theTreeCSV_map[sample] = theTreeCSV_list.back();
   
//  treeCSV_jetN        = -10000;
   
//  for(int i=0;i<1000;i++)
//     {	
//	treeCSV_jetPt[i]              = -10000;
//	treeCSV_jetEta[i]      = -10000;
//	treeCSV_jetCSV[i]      = -10000;
//	treeCSV_jetFlav[i]      = -10000;
//     }
           
//  treeCSV_EvtWeight      = -10000;
//  treeCSV_Channel        = -10000;
//  treeCSV_SampleType     = -10000;
}


//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1F binning
//creates one histograms per channel
//-------------------------------------------------------------

void TreeReader::addHisto(TString var, TString selstep, TString sample, int nbins, float min, float max){
  
  
  TString name_mmm =  var+"_mumumu_"+selstep+"__"+sample;
  TH1F * thehisto_mmm = new TH1F(name_mmm,name_mmm,nbins,min,max);
  thehisto_mmm->Sumw2();
  histo_list_mmm.push_back(thehisto_mmm);
  histo_map_mmm[name_mmm.Data()] = numb_histo;
  
  
  TString name_mme =  var+"_mumue_"+selstep+"__"+sample;
  TH1F * thehisto_mme = new TH1F(name_mme,name_mme,nbins,min,max);
  thehisto_mme->Sumw2();
  histo_list_mme.push_back(thehisto_mme);
  histo_map_mme[name_mme.Data()] = numb_histo;
  
  TString name_eem =  var+"_eemu_"+selstep+"__"+sample;
  TH1F * thehisto_eem = new TH1F(name_eem,name_eem,nbins,min,max);
  thehisto_eem->Sumw2();
  histo_list_eem.push_back(thehisto_eem);
  histo_map_eem[name_eem.Data()] = numb_histo;
  
  TString name_eee =  var+"_eee_"+selstep+"__"+sample;
  TH1F * thehisto_eee = new TH1F(name_eee,name_eee,nbins,min,max);
  thehisto_eee->Sumw2();
  histo_list_eee.push_back(thehisto_eee);
  histo_map_eee[name_eee.Data()] = numb_histo;
  
  //cout << "adding an histo with name " << name_mmm  << " and map integer " << numb_histo << endl;
  //cout << "adding an histo with name " << name_mme  << " and map integer " << numb_histo << endl;
  //cout << "adding an histo with name " << name_eem  << " and map integer " << numb_histo << endl;
  //cout << "adding an histo with name " << name_eee  << " and map integer " << numb_histo << endl;
  
  numb_histo++;
  
  
}

//-------------------------------------------------------------
//fill histograms
//first parameter is the channel,
//second parameter is the variable name,
//third parameter is the selection step (like "afterleptsel")
//forths parameter is the sample name (like "Z)
//others are value and weight
//-------------------------------------------------------------
void TreeReader::fillHisto(TString channel, TString var, TString selstep, TString sample, float val, float weight){
  TString name = var+"_"+channel+"_"+selstep+"__"+sample;

  
  if(channel == "mumumu" && histo_map_mmm[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }
  if(channel == "mumue" && histo_map_mme[name.Data()] == 0) {
    cout << "  WARNING trying to fill a non existing histograms " << endl;
    cout << "  please check the naming conventions " << endl;
    cout << "  histo name "  << name << endl;
  }
  if(channel == "eemu" && histo_map_eem[name.Data()] == 0) {
    cout << "  WARNING trying to fill a non existing histograms " << endl;
    cout << "  please check the naming conventions " << endl;
    cout << "  histo name "  << name << endl;
  }
  if(channel == "eee" && histo_map_eee[name.Data()] == 0) {
    cout << "  WARNING trying to fill a non existing histograms " << endl;
    cout << "  please check the naming conventions " << endl;
    cout << "  histo name "  << name << endl;
  }
  
  
  /*if(sample == "tZq_leptup") {
    cout << "val " << val << " weight  " << weight << endl; 
    cout << "name " << name << endl;  
    cout << "channel " << channel << endl;
  }*/
  if(channel == "mumumu")     histo_list_mmm[histo_map_mmm[name.Data()]]->Fill(val, weight);
  else if(channel == "mumue") histo_list_mme[histo_map_mme[name.Data()]]->Fill(val, weight);
  else if(channel == "eemu")  histo_list_eem[histo_map_eem[name.Data()]]->Fill(val, weight);
  else if(channel == "eee")   histo_list_eee[histo_map_eee[name.Data()]]->Fill(val, weight);  
  
}

//-------------------------------------------------------------
//determine the decay channel
//-------------------------------------------------------------
TString TreeReader::determineChannel(int leptflav1, int leptflav2, int leptflav3){
  
  TString thechannel = "";
  
  if     (abs(leptflav1) == 13 && abs(leptflav2) == 13 && abs(leptflav3) == 13 ) thechannel = "mumumu";
  else if(abs(leptflav1) == 13 && abs(leptflav2) == 13 && abs(leptflav3) == 11 ) thechannel = "mumue";
  else if(abs(leptflav1) == 11 && abs(leptflav2) == 11 && abs(leptflav3) == 13 ) thechannel = "eemu";
  else if(abs(leptflav1) == 11 && abs(leptflav2) == 11 && abs(leptflav3) == 11 ) thechannel = "eee";
  
  //cout << "leptflav1 " << leptflav1 << "  leptflav2 " << leptflav2 << "  leptflav3 " << leptflav3 <<  "  the channel " << thechannel << endl;
  if(thechannel == "") cout << "WARNING no channel found, please check the lepton flavor " << endl;
  return thechannel;
  
  
}

void TreeReader::deleteHisto(){
   cout << __LINE__ << endl;
   /*for(unsigned int i=0; i<histo_list_mmm.size(); i++){
     
     delete  histo_list_mmm[i];
     delete  histo_list_mme[i];
     delete  histo_list_eem[i];
     delete  histo_list_eee[i];
     
     
   }*/
   cout << __LINE__ << endl;
  //delete TheTree;
   cout << __LINE__ << endl;
  
  
}

//SetUp CSV reweighting
void TreeReader::SetUpCSVreweighting(){

  // Do not set it up if we're running on collision data
  if(isData){ return; }
  f_CSVwgt_HF = new TFile ("../BTagCSV/CSVRW/csv_rwt_hf.root");
  f_CSVwgt_LF = new TFile ("../BTagCSV/CSVRW/csv_rwt_lf.root");


  // CSV reweighting
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<6; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<6; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)f_CSVwgt_LF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

}

// Get event weight
double TreeReader::GetCSVweight(const int iSys, int jet_n,
				float *jet_pt,float *jet_eta,float *jet_btagdiscri,int *jet_flav){
  if (isData) return 1.0;

  int iSysHF = 0;
  int iSysC = 0;
  int iSysLF = 0;
   
  // systematic variations
  // 1 HF/LF JESup iSysHF=1 iSysC=0 iSysLF=1
  // 2 HF/LF JESdown iSysHF=2 iSysC=0 iSysLF=2
  // 3 HF CSVLFup iSysHF=3 iSysC=0 iSysLF=0
  // 4 HF CSVLFdown iSysHF=4 iSysC=0 iSysLF=0
  // 5 HF CSVHFStats1up iSysHF=5 iSysC=0 iSysLF=0
  // 6 HF CSVHFStats1down iSysHF=6 iSysC=0 iSysLF=0
  // 7 HF CSVHFStats2up iSysHF=7 iSysC=0 iSysLF=0
  // 8 HF CSVHFStats2down iSysHF=8 iSysC=0 iSysLF=0
  // 9 C CSVCErr1up iSysHF=0 iSysC=1 iSysLF=0
  // 10 C CSVCErr1down iSysHF=0 iSysC=2 iSysLF=0
  // 11 C CSVCErr2up iSysHF=0 iSysC=3 iSysLF=0
  // 12 C CSVCErr2down iSysHF=0 iSysC=4 iSysLF=0
  // 13 LF CSVHFup iSysHF=0 iSysC=0 iSysLF=3
  // 14 LF CSVHFdown iSysHF=0 iSysC=0 iSysLF=4
  // 15 LF CSVLFStats1up iSysHF=0 iSysC=0 iSysLF=5
  // 16 LF CSVLFStats1down iSysHF=0 iSysC=0 iSysLF=6
  // 17 LF CSVLFStats2up iSysHF=0 iSysC=0 iSysLF=7
  // 18 LF CSVLFStats2down iSysHF=0 iSysC=0 iSysLF=8
  if( iSys == 1 ) {iSysHF=1; iSysC=0; iSysLF=1;}
  else if( iSys == 2 ) {iSysHF=2; iSysC=0; iSysLF=2;}
  else if( iSys == 3 ) {iSysHF=3; iSysC=0; iSysLF=0;}
  else if( iSys == 4 ) {iSysHF=4; iSysC=0; iSysLF=0;}
  else if( iSys == 5 ) {iSysHF=5; iSysC=0; iSysLF=0;}
  else if( iSys == 6 ) {iSysHF=6; iSysC=0; iSysLF=0;}
  else if( iSys == 7 ) {iSysHF=7; iSysC=0; iSysLF=0;}
  else if( iSys == 8 ) {iSysHF=8; iSysC=0; iSysLF=0;}
  else if( iSys == 9 ) {iSysHF=0; iSysC=1; iSysLF=0;}
  else if( iSys == 10 ) {iSysHF=0; iSysC=2; iSysLF=0;}
  else if( iSys == 11 ) {iSysHF=0; iSysC=3; iSysLF=0;}
  else if( iSys == 12 ) {iSysHF=0; iSysC=4; iSysLF=0;}
  else if( iSys == 13 ) {iSysHF=0; iSysC=0; iSysLF=3;}
  else if( iSys == 14 ) {iSysHF=0; iSysC=0; iSysLF=4;}
  else if( iSys == 15 ) {iSysHF=0; iSysC=0; iSysLF=5;}
  else if( iSys == 16 ) {iSysHF=0; iSysC=0; iSysLF=6;}
  else if( iSys == 17 ) {iSysHF=0; iSysC=0; iSysLF=7;}
  else if( iSys == 18 ) {iSysHF=0; iSysC=0; iSysLF=8;}

//  iSysHF=1; // JESup
//  iSysHF=2; // JESdown
//  iSysHF=3; // CSVLFup
//  iSysHF=4; // CSVLFdown
//  iSysHF=5; // CSVHFStats1up
//  iSysHF=6; // CSVHFStats1down
//  iSysHF=7; // CSVHFStats2up
//  iSysHF=8; // CSVHFStats2down

//  iSysC=1; // CSVCErr1up
//  iSysC=2; // CSVCErr1down
//  iSysC=3; // CSVCErr2up
//  iSysC=4; // CSVCErr2down

//  iSysLF=1; // JESup
//  iSysLF=2; // JESdown
//  iSysLF=3; // CSVHFup
//  iSysLF=4; // CSVHFdown
//  iSysLF=5; // CSVLFStats1up
//  iSysLF=6; // CSVLFStats1down
//  iSysLF=7; // CSVLFStats2up
//  iSysLF=8; // CSVLFStats2down


  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int ij=0;ij<jet_n;ij++){
     
    double csv = jet_btagdiscri[ij];
    double jetPt = jet_pt[ij];
    double jetAbsEta = fabs( jet_eta[ij] );
    int flavor = abs( jet_flav[ij] );

    if( jetPt < 20. || jetAbsEta > 2.5 ) continue;
     
    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100 && jetPt<160) iPt = 4;
    else if (jetPt >=160 && jetPt<10000) iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6) iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;
    // kskovpen hack - should change eta cut to 2.4 !
    else if ( jetAbsEta>=2.41 && jetAbsEta<=2.5) iEta = 2;
     
    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;

    }
    else {
      if (iPt >=3) iPt=3;       /// [>60] 
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;


    }
  }


  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  return csvWgtTotal;
}



