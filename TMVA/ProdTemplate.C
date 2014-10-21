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



void setNjetSel(TH1F * thehisto, int njetsel){
  
  if(njetsel!=0) {
    for(int i=1; i<=njetsel; i++){
      thehisto->SetBinContent(i, 0.0000000001);
      thehisto->SetBinError(i, 0.0000000001);
    }
  }
  
}
 

void ProdTemplate(TString inputdistrib, std::vector<TString> sampleList, std::vector<TString> stytList, TString intputfilename, TString outputfilename){
  
  
  
  
  
  
  
  
  
  TFile * inputfile	  = new TFile( intputfilename.Data() );


  TH1F * histBdt_DataMu   = (TH1F*)inputfile->Get(  (inputdistrib+"__DataMu").Data())->Clone();
  TH1F * histBdt_DataEG   = (TH1F*)inputfile->Get(  (inputdistrib+"__DataEG").Data())->Clone();
  TH1F * histBdt_DataMuEG = (TH1F*)inputfile->Get(  (inputdistrib+"__DataMuEG").Data())->Clone();
  
  
  histBdt_DataMu->Add(histBdt_DataEG);
  histBdt_DataMu->Add(histBdt_DataMuEG);
  
  
  std::vector<TH1F* > distrib_MC;
  std::vector<TH1F* > distrib_MC_sys;
  
  
  for(unsigned int i=0; i<sampleList.size(); i++){ 
    cout << inputdistrib+"__"+sampleList[i]  << endl;
    distrib_MC.push_back( (TH1F*)inputfile->Get( (inputdistrib+"__"+sampleList[i]).Data() )->Clone() ); 
  } 
  
  
  
  for(unsigned int i=0; i<sampleList.size(); i++){
    if(sampleList[i] == "DataMu" || sampleList[i] == "DataEG" || sampleList[i] == "DataMuEG" ) continue;
    for(unsigned int j=0; j<stytList.size(); j++){
        cout <<  inputdistrib+"__"+sampleList[i]+""+stytList[j]   << endl;
	distrib_MC_sys.push_back( (TH1F*)inputfile->Get( (inputdistrib+"__"+sampleList[i]+""+stytList[j]).Data() )->Clone() );  
    }
  }  
  
  
  

   
  //for(unsigned int i=0; i<stytList.size(); i++) outputfilename+="_"+stytList[i];
  outputfilename ="template_theta.root"; 
  
  
  
  TFile * outputfile = new TFile(outputfilename, "recreate" );
  outputfile->cd();

  //distrib__DATA->Write(); 
  histBdt_DataMu->Write((inputdistrib+"__DATA").Data() );
  
  
  for(unsigned int i=0; i<distrib_MC.size(); i++)               distrib_MC[i]->Write();
  for(unsigned int i=0; i<distrib_MC_sys.size(); i++)           distrib_MC_sys[i]->Write();
  
  
}

 
void ProdTemplate(){
  
  std::vector<TString> sampleList;
  
  sampleList.push_back("tZq");
  sampleList.push_back("TTZ");
  sampleList.push_back("TTW");
  sampleList.push_back("TT");
  //sampleList.push_back("DYToLL_M10-50");
  sampleList.push_back("Zjets");
  //sampleList.push_back("Wjets");
  sampleList.push_back("WZ");
  sampleList.push_back("ZZ");
  /*sampleList.push_back("TbarsChan");
  sampleList.push_back("TsChan");
  sampleList.push_back("TtChan");
  sampleList.push_back("TbartChan");
  sampleList.push_back("TtW");
  sampleList.push_back("TbartW");*/
  
  
  std::vector<TString> systList;
  //systList.push_back("");  
  systList.push_back("__lept__plus");
  systList.push_back("__lept__minus");
  //systList.push_back("__trig__plus");
  //systList.push_back("__trig__minus");
  //systList.push_back("__PDF__plus");
  //systList.push_back("__PDF__minus");
  systList.push_back("__jes__plus");
  systList.push_back("__jes__minus");
  systList.push_back("__jer__plus"); 
  systList.push_back("__jer__minus");
  systList.push_back("__metuncls__plus");
  systList.push_back("__metuncls__minus"); 
  
  
  
  
  ProdTemplate("MVA", sampleList,  systList,  "outputroot/merged_all.root", "merged_theta");
  
  
  
  
  
}
