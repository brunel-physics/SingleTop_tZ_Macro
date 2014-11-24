

{

 gROOT->ProcessLine(".L PlotStack.C+");
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);
 
 
  
  //---------------------------
  //define list of data samples 
  std::vector<TString> dataSample_list;
  dataSample_list.push_back("DATA");
  
  //------------------------------------
  //define list of data-driven templates
  std::vector< TString > dataDrivenTemplates_list; 
  //dataDrivenTemplates_list.push_back( "DataEGZenriched"  );
  //dataDrivenTemplates_list.push_back( "DataMuEGZenriched");
  //dataDrivenTemplates_list.push_back( "DataMuZenriched"  );
  
  
 //-------------------------
 //define list of MC samples 
  std::vector<TString> mcSample_list;
  std::vector<int> colorVector;
  //mcSample_list.push_back("TbarsChan");     colorVector.push_back(kMagenta);
  //mcSample_list.push_back("TsChan");        colorVector.push_back(kMagenta);
  //mcSample_list.push_back("TtChan");        colorVector.push_back(kMagenta);
  //mcSample_list.push_back("TbartChan");     colorVector.push_back(kMagenta);
  //mcSample_list.push_back("TtW");           colorVector.push_back(kMagenta);
  //mcSample_list.push_back("TbartW");        colorVector.push_back(kMagenta);
  mcSample_list.push_back("TT");            colorVector.push_back(kRed-7);
  mcSample_list.push_back("TTZ");           colorVector.push_back(kRed+1);
  mcSample_list.push_back("TTW");           colorVector.push_back(kRed+1);
  mcSample_list.push_back("WZ");            colorVector.push_back(13);
  mcSample_list.push_back("ZZ");            colorVector.push_back(13);
  mcSample_list.push_back("WZHF");            colorVector.push_back(12);
  //mcSample_list.push_back("DYToLL_M10-50"); colorVector.push_back(kAzure-2);
  mcSample_list.push_back("Zjets");         colorVector.push_back(kAzure-2);
  mcSample_list.push_back("tZq");           colorVector.push_back(kGreen+2);
  //mcSample_list.push_back("WW");            colorVector.push_back(13);
  
  
  
 //-----------------------------
 //define list of signal samples 
  std::vector<TString> signalSample_list;
  
  
  
  
 //--------------------------
 //define list of systematics
 std::vector<TString> syst_list;
 syst_list.push_back("");
  
  
 
 //--------------------------
 //define list of systematics
 std::vector<TString> selectionStep_list;
 selectionStep_list.push_back("");
 
 
 //------------------------
 //define list of variables
 std::vector<TString> variables_list;
 variables_list.push_back("MVA_all");
 variables_list.push_back("MVA_mumumu");
 variables_list.push_back("MVA_mumue");
 variables_list.push_back("MVA_eemu");
 variables_list.push_back("MVA_eee");
 //variables_list.push_back("tree_cosThetaStar");  
/* variables_list.push_back("tree_topMass_all");       
 //variables_list.push_back("tree_totMass");       
 //variables_list.push_back("tree_deltaPhilb");    
 //variables_list.push_back("tree_deltaRlb");      
 //variables_list.push_back("tree_deltaRTopZ");    
 variables_list.push_back("tree_asym_all");	      
 //variables_list.push_back("tree_Zpt");	      
 //variables_list.push_back("tree_ZEta");	      
 //variables_list.push_back("tree_topPt");	      
 //variables_list.push_back("tree_topEta");        
 //variables_list.push_back("tree_NJets");	      
 variables_list.push_back("tree_NBJets_all");        
 //variables_list.push_back("tree_deltaRZl");      
 //variables_list.push_back("tree_deltaPhiZmet");  
 variables_list.push_back("tree_btagDiscri_all");    
 variables_list.push_back("tree_totPt_all");	      
 //variables_list.push_back("tree_totEta");        
 variables_list.push_back("tree_leptWPt_all");       
 //variables_list.push_back("tree_leptWEta");      
 variables_list.push_back("tree_leadJetPt_all");     
 //variables_list.push_back("tree_leadJetEta");    
 //variables_list.push_back("tree_deltaRZleptW");  
 //variables_list.push_back("tree_deltaPhiZleptW");
 //variables_list.push_back("tree_met");	      
 //variables_list.push_back("tree_mTW");	      

*/
   
   for(int ivar=0; ivar < variables_list.size(); ivar++){
        
       PlotStack(variables_list[ivar], 0, dataSample_list,  mcSample_list, signalSample_list, colorVector, dataDrivenTemplates_list);
       
    }
    
   
   
   
   
   
 }
