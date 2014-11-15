{
  gROOT->ProcessLine(".L TreeReader.C+");
  
  
  std::vector<TString > systlist;
  systlist.push_back("");		
  /*systlist.push_back("lept__plus");
  systlist.push_back("lept__minus" );    
  //systlist.push_back("trig__plus");
  //systlist.push_back("trig__minus" );  
  //systlist.push_back("PDF__plus");     
  //systlist.push_back("PDF__minus");    
  systlist.push_back("jes__plus");
  systlist.push_back("jes__minus");      
  systlist.push_back("jer__plus");
  systlist.push_back("jer__minus");   
  systlist.push_back("metuncls__plus");
  systlist.push_back("metuncls__minus");
  */
  
  
  std::vector<TString > datalist;
  datalist.push_back("");	
  
  TTree* tree=0;
  
  TreeReader * tree_DataMu = new TreeReader(tree, "DataMu", datalist);
  tree_DataMu.Loop("DataMu");
  delete tree_DataMu;
  
  
  
  TreeReader * tree_DataMuEG = new TreeReader(tree, "DataMuEG", datalist);
  tree_DataMuEG.Loop("DataMuEG");
  delete tree_DataMuEG;


  TreeReader * tree_DataEG = new TreeReader(tree, "DataEG", datalist);
  tree_DataEG.Loop("DataEG");
  delete tree_DataEG;
  
  
  
  TreeReader * tree_DataEGZenriched = new TreeReader(tree, "DataEGZenriched", datalist);
  tree_DataEGZenriched.Loop("DataEGZenriched");
  delete tree_DataEGZenriched;
  
  
  TreeReader * tree_DataMuEGZenriched = new TreeReader(tree, "DataMuEGZenriched", datalist);
  tree_DataMuEGZenriched.Loop("DataMuEGZenriched");
  delete tree_DataMuEGZenriched;
  
  
  TreeReader * tree_DataMuZenriched = new TreeReader(tree, "DataMuZenriched", datalist);
  tree_DataMuZenriched.Loop("DataMuZenriched");
  delete tree_DataMuZenriched;
  
  
  TreeReader * tree_tZq = new TreeReader(tree, "tZq", systlist);
  tree_tZq.Loop("tZq");
  delete tree_tZq;
  
  
  
  TreeReader * tree_TTZ = new TreeReader(tree, "TTZ", systlist);
  tree_TTZ.Loop("TTZ");
  delete tree_TTZ;
  
  
  TreeReader * tree_TTW = new TreeReader(tree, "TTW", systlist);
  tree_TTW.Loop("TTW");
  delete tree_TTW;
  
  
  TreeReader * tree_TT = new TreeReader(tree, "TT", systlist);
  tree_TT.Loop("TT");
  delete tree_TT;
  
  TreeReader * tree_DYToLL_M10_50 = new TreeReader(tree, "DYToLL_M10-50", systlist);
  tree_DYToLL_M10_50.Loop("DYToLL_M10-50");
  delete tree_DYToLL_M10_50;
  
  
  
  TreeReader * tree_Zjets = new TreeReader(tree, "Zjets", systlist);
  tree_Zjets.Loop("Zjets");
  delete tree_Zjets;
  
  TreeReader * tree_Wjets = new TreeReader(tree, "Wjets", systlist);
  tree_Wjets.Loop("Wjets");
  delete tree_Wjets;
  
  TreeReader * tree_WW = new TreeReader(tree, "WW", systlist);
  tree_WW.Loop("WW");
  delete tree_WW;
  
  
  //one for WZ +light
  TreeReader * tree_WZ = new TreeReader(tree, "WZ", systlist);
  tree_WZ.Loop("WZ");
  delete tree_WZ;
  
  //one for WZ + HF
  TreeReader * tree_WZ = new TreeReader(tree, "WZHF", systlist);
  tree_WZ.Loop("WZHF");
  delete tree_WZ;

  TreeReader * tree_ZZ = new TreeReader(tree, "ZZ", systlist);
  tree_ZZ.Loop("ZZ");
  delete tree_ZZ;

  
  TreeReader * tree_TbarsChan = new TreeReader(tree, "TbarsChan", systlist);
  tree_TbarsChan.Loop("TbarsChan");
  delete tree_TbarsChan;

  
  TreeReader * tree_TsChan = new TreeReader(tree, "TsChan", systlist);
  tree_TsChan.Loop("TsChan");
  delete tree_TsChan;

  
  TreeReader * tree_TtChan = new TreeReader(tree, "TtChan", systlist);
  tree_TtChan.Loop("TtChan");
  delete tree_TtChan;


  TreeReader * tree_TbartChan = new TreeReader(tree, "TbartChan", systlist);
  tree_TbartChan.Loop("TbartChan");
  delete tree_TbartChan;

  
  TreeReader * tree_TtW = new TreeReader(tree, "TtW", systlist);
  tree_TtW.Loop("TtW");
  delete tree_TtW;

  
  TreeReader * tree_TbartW = new TreeReader(tree, "TbartW", systlist);
  tree_TbartW.Loop("TbartW");
  delete tree_TbartW;

  
  TreeReader * tree_TbartW = new TreeReader(tree, "TbartW", systlist);
  tree_TbartW.Loop("TbartW");
  delete tree_TbartW;

  
  
  
}
