#include "TString.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>

bool doDD = true;

void setNjetSel(TH1F * thehisto, int njetsel){
  
  if(njetsel!=0) {
    for(int i=1; i<=njetsel; i++){
      thehisto->SetBinContent(i, 0.0000000001);
      thehisto->SetBinError(i, 0.0000000001);
    }
  }
  
}
 

bool ProdTemplate(TString inputdistrib, std::vector<TString> sampleList, std::vector<TString> stytList, TString intputfilename, TString outputfilename){
  
  
  TFile * inputfile	  = new TFile( intputfilename.Data() );

  TH1F * histBdt_DataMu   = (TH1F*)inputfile->Get(  (inputdistrib+"__DataMu").Data())->Clone();
  TH1F * histBdt_DataEG   = (TH1F*)inputfile->Get(  (inputdistrib+"__DataEG").Data())->Clone();
  TH1F * histBdt_DataMuEG = (TH1F*)inputfile->Get(  (inputdistrib+"__DataMuEG").Data())->Clone();
  
  histBdt_DataMu->Sumw2();
  histBdt_DataEG->Sumw2();
  histBdt_DataMuEG->Sumw2();
  histBdt_DataMu->Add(histBdt_DataEG);
  histBdt_DataMu->Add(histBdt_DataMuEG);
  
  
  std::vector<TH1F* > distrib_MC;
  std::vector<TH1F* > distrib_MC_sys;
  
  //deal with nominal templates
  for(unsigned int i=0; i<sampleList.size(); i++){ 
    cout << inputdistrib+"__"+sampleList[i]  << endl;
    if(doDD && sampleList[i] == "Zjets") continue;
    TH1F * tmp = (TH1F*)inputfile->Get( (inputdistrib+"__"+sampleList[i]).Data() )->Clone();
    if(tmp == 0){ cout << "non existing histogram " <<  inputdistrib+"__"+sampleList[i] << endl; return false;}
    tmp->Sumw2();
//    if(sampleList[i] == "tZq") tmp->Scale();
    distrib_MC.push_back( tmp ); 
  } 
  
  
  //deal with systematic templates
  for(unsigned int i=0; i<sampleList.size(); i++){
    if(sampleList[i] == "DataMu" || sampleList[i] == "DataEG" || sampleList[i] == "DataMuEG" ) continue;
    if(doDD && sampleList[i] == "Zjets") continue;
    for(unsigned int j=0; j<stytList.size(); j++){
        cout <<  inputdistrib+"__"+sampleList[i]+""+stytList[j]   << endl;
	TH1F * tmp = (TH1F*)inputfile->Get( (inputdistrib+"__"+sampleList[i]+""+stytList[j]).Data() )->Clone() ;
        tmp->Sumw2();
//	if(sampleList[i] == "tZq") tmp->Scale(3.);
	distrib_MC_sys.push_back( tmp);  
    }
  }  
  
  
  //deal with dd templates
  TH1F * dataMu_Zenriched   = (TH1F*)inputfile->Get( (inputdistrib+"__"+"DataMuZenriched").Data() );
  TH1F * dataMuEG_Zenriched = (TH1F*)inputfile->Get( (inputdistrib+"__"+"DataMuEGZenriched").Data() );
  TH1F * dataEG_Zenriched   = (TH1F*)inputfile->Get( (inputdistrib+"__"+"DataEGZenriched").Data() );
  dataMu_Zenriched->Sumw2();   
  dataMuEG_Zenriched->Sumw2(); 
  dataEG_Zenriched->Sumw2();   
  TH1F * zJets_forNorm      = (TH1F*)inputfile->Get( (inputdistrib+"__"+"Zjets").Data() );
   
  dataMu_Zenriched->Add(dataMuEG_Zenriched);
  dataMu_Zenriched->Add(dataEG_Zenriched);
  
  if(dataMu_Zenriched->Integral() > 0.0001) dataMu_Zenriched->Scale(zJets_forNorm->Integral()/dataMu_Zenriched->Integral());
   
  //for(unsigned int i=0; i<stytList.size(); i++) outputfilename+="_"+stytList[i];
  outputfilename = "TemplateRootFiles/"+inputdistrib+"_theta.root"; 
  
  
  
  TFile * outputfile = new TFile(outputfilename.Data(), "recreate" );
  outputfile->cd();

  //distrib__DATA->Write(); 
  histBdt_DataMu->Write((inputdistrib+"__DATA").Data() );
  
  
  for(unsigned int i=0; i<distrib_MC.size(); i++)               distrib_MC[i]->Write();
  if(doDD) dataMu_Zenriched->Write((inputdistrib+"__"+"Zjets").Data());
  for(unsigned int i=0; i<distrib_MC_sys.size(); i++)           distrib_MC_sys[i]->Write();
  
  return true;
}

 
void ProdTemplate(){
  
  std::vector<TString> sampleList;

  sampleList.push_back("tZq");
  sampleList.push_back("TTZ");
  sampleList.push_back("TTW");
  //  sampleList.push_back("TT");
  //  sampleList.push_back("DYToLL_M10-50");
  sampleList.push_back("Zjets");
  //  sampleList.push_back("Wjets");
  sampleList.push_back("WZ");
  sampleList.push_back("ZZ");
  //  sampleList.push_back("TbarsChan");
  //  sampleList.push_back("TsChan");
  //  sampleList.push_back("TtChan");
  //  sampleList.push_back("TbartChan");
  //  sampleList.push_back("TtW");
  //  sampleList.push_back("TbartW");
  //  samplelist.push_back("WW");
  
  std::vector<TString> systList;

  systList.push_back("");
  systList.push_back("__pileup__plus");
  systList.push_back("__pileup__minus");
  systList.push_back("__trig__plus");
  systList.push_back("__trig__minus");
  systList.push_back("__bTag__plus");
  systList.push_back("__bTag__minus");
  systList.push_back("__pdf__plus");
  systList.push_back("__pdf__minus");
  systList.push_back("__jes__plus");
  systList.push_back("__jes__minus");
  systList.push_back("__jer__plus"); 
  systList.push_back("__jer__minus");
  systList.push_back("__met__plus");
  systList.push_back("__met__minus"); 


  //Original Systematics List

  //  systList.push_back("");  
//  systList.push_back("__lept__plus");
//  systList.push_back("__lept__minus");
  //systList.push_back("__trig__plus");
  //systList.push_back("__trig__minus");
  //systList.push_back("__PDF__plus");
  //systList.push_back("__PDF__minus");
//  systList.push_back("__jes__plus");
//  systList.push_back("__jes__minus");
//  systList.push_back("__jer__plus"); 
//  systList.push_back("__jer__minus");
//  systList.push_back("__metuncls__plus");
//  systList.push_back("__metuncls__minus"); 
  

  ProdTemplate("MVA_all",    sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("MVA_mumumu", sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta");
  ProdTemplate("MVA_emumu",  sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta");
  ProdTemplate("MVA_eemu",   sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta");
  ProdTemplate("MVA_eee",    sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta");


  ProdTemplate("tree_topMass_all",          sampleList,  systList,  "outputroot/output_merged.root", "merged_theta"); 
  //ProdTemplate("tree_Zpt_all",              sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  ProdTemplate("tree_Zeta_all",             sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  //ProdTemplate("tree_topPt_all",            sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");    
  //  ProdTemplate("tree_topEta_all",           sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");   
  //ProdTemplate("tree_NJets_all",            sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");    
  //ProdTemplate("tree_NBJets_all",           sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");   
  ProdTemplate("tree_btagDiscri_all",       sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  ProdTemplate("tree_totPt_all", sampleList, systList, "outputroot/output_merged.root", "merged_theta");
  //ProdTemplate("tree_totEta_all", sampleList, systList, "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_leptWPt_all",          sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");  
  ProdTemplate("tree_leptWEta_all",         sampleList,  systList,  "outputroot/output_merged.root", "merged_theta"); 
  //  ProdTemplate("tree_leadJetPt_all",        sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  ProdTemplate("tree_leadJetEta_all",       sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  //  ProdTemplate("tree_met_all",        sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  //  ProdTemplate("tree_mTW_all",        sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  

  //User Added
  //  ProdTemplate("tree_totPtVec_all",        sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  //  ProdTemplate("tree_totVecM_all",        sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  //  ProdTemplate("tree_lepPt_all",        sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	  
  //  ProdTemplate("tree_lepMetPt_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_totPt2Jet_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_leadJetbTag_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_secJetbTag_all",   sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");   
  ProdTemplate("tree_secJetPt_all",              sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	
  //  ProdTemplate("tree_secJetEta_all",              sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");	
  ProdTemplate("tree_wzdelR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_jjdelR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_zjminR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_zlepWdelPhi_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_ZmetdelPhi_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_ZlepWdelR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_lbDelR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_lbDelPhi_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_zlb1DelR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  ProdTemplate("tree_zlb1DelPhi_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");  
  //  ProdTemplate("tree_zlb2DelR_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
    ProdTemplate("tree_zlb2DelPhi_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");  
  ProdTemplate("tree_totHt_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");  
  //  ProdTemplate("tree_lepHt_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_jetHt_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_lepMetHt_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_totHtOverPt_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  //  ProdTemplate("tree_zMass_all",     sampleList,  systList,  "outputroot/output_merged.root", "merged_theta");
  
  // mumumu

  ProdTemplate("tree_topMass_mumumu",          sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu"); 
  ProdTemplate("tree_Zeta_mumumu",             sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");	  
  ProdTemplate("tree_btagDiscri_mumumu",       sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");	  
  ProdTemplate("tree_totPt_mumumu", sampleList, systList, "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");
  ProdTemplate("tree_leptWPt_mumumu",          sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");  
  ProdTemplate("tree_leptWEta_mumumu",         sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu"); 
  ProdTemplate("tree_leadJetEta_mumumu",       sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");	  
  ProdTemplate("tree_secJetbTag_mumumu",   sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");   
  ProdTemplate("tree_secJetPt_mumumu",              sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");	
  ProdTemplate("tree_wzdelR_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");
  ProdTemplate("tree_jjdelR_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");
  ProdTemplate("tree_lbDelR_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");
  ProdTemplate("tree_lbDelPhi_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");
  ProdTemplate("tree_zlb1DelR_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");
  ProdTemplate("tree_zlb1DelPhi_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");  
  ProdTemplate("tree_zlb2DelPhi_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");  
  ProdTemplate("tree_totHt_mumumu",     sampleList,  systList,  "outputroot/output_merged_mumumu.root", "merged_theta_mumumu");  
 

  // emumu

  ProdTemplate("tree_topMass_emumu",          sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu"); 
  ProdTemplate("tree_Zeta_emumu",             sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");	  
  ProdTemplate("tree_btagDiscri_emumu",       sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");	  
  ProdTemplate("tree_totPt_emumu", sampleList, systList, "outputroot/output_merged_emumu.root", "merged_theta_emumu");
  ProdTemplate("tree_leptWPt_emumu",          sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");  
  ProdTemplate("tree_leptWEta_emumu",         sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu"); 
  ProdTemplate("tree_leadJetEta_emumu",       sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");	  
  ProdTemplate("tree_secJetbTag_emumu",   sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");   
  ProdTemplate("tree_secJetPt_emumu",              sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");	
  ProdTemplate("tree_wzdelR_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");
  ProdTemplate("tree_jjdelR_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");
  ProdTemplate("tree_lbDelR_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");
  ProdTemplate("tree_lbDelPhi_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");
  ProdTemplate("tree_zlb1DelR_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");
  ProdTemplate("tree_zlb1DelPhi_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");  
  ProdTemplate("tree_zlb2DelPhi_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");  
  ProdTemplate("tree_totHt_emumu",     sampleList,  systList,  "outputroot/output_merged_emumu.root", "merged_theta_emumu");  
 
  // eemu

  ProdTemplate("tree_topMass_eemu",          sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu"); 
  ProdTemplate("tree_Zeta_eemu",             sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");	  
  ProdTemplate("tree_btagDiscri_eemu",       sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");	  
  ProdTemplate("tree_totPt_eemu", sampleList, systList, "outputroot/output_merged_eemu.root", "merged_theta_eemu");
  ProdTemplate("tree_leptWPt_eemu",          sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");  
  ProdTemplate("tree_leptWEta_eemu",         sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu"); 
  ProdTemplate("tree_leadJetEta_eemu",       sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");	  
  ProdTemplate("tree_secJetbTag_eemu",   sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");   
  ProdTemplate("tree_secJetPt_eemu",              sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");	
  ProdTemplate("tree_wzdelR_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");
  ProdTemplate("tree_jjdelR_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");
  ProdTemplate("tree_lbDelR_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");
  ProdTemplate("tree_lbDelPhi_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");
  ProdTemplate("tree_zlb1DelR_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");
  ProdTemplate("tree_zlb1DelPhi_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");  
  ProdTemplate("tree_zlb2DelPhi_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");  
  ProdTemplate("tree_totHt_eemu",     sampleList,  systList,  "outputroot/output_merged_eemu.root", "merged_theta_eemu");  
 
  //eee

  ProdTemplate("tree_topMass_eee",          sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee"); 
  ProdTemplate("tree_Zeta_eee",             sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");	  
  ProdTemplate("tree_btagDiscri_eee",       sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");	  
  ProdTemplate("tree_totPt_eee", sampleList, systList, "outputroot/output_merged_eee.root", "merged_theta_eee");
  ProdTemplate("tree_leptWPt_eee",          sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");  
  ProdTemplate("tree_leptWEta_eee",         sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee"); 
  ProdTemplate("tree_leadJetEta_eee",       sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");	  
  ProdTemplate("tree_secJetbTag_eee",   sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");   
  ProdTemplate("tree_secJetPt_eee",              sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");	
  ProdTemplate("tree_wzdelR_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");
  ProdTemplate("tree_jjdelR_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");
  ProdTemplate("tree_lbDelR_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");
  ProdTemplate("tree_lbDelPhi_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");
  ProdTemplate("tree_zlb1DelR_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");
  ProdTemplate("tree_zlb1DelPhi_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");  
  ProdTemplate("tree_zlb2DelPhi_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");  
  ProdTemplate("tree_totHt_eee",     sampleList,  systList,  "outputroot/output_merged_eee.root", "merged_theta_eee");  

}
