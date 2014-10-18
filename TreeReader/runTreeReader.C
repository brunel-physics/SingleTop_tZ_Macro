{
  gROOT->ProcessLine(".L TreeReader.C+");
  
  TTree* tree=0;
  
  TreeReader * tree_DataMu = new TreeReader(tree, "DataMu");
  tree_DataMu.Loop("DataMu");
  delete tree_DataMu;
  
  
  
  TreeReader * tree_DataMuEG = new TreeReader(tree, "DataMuEG");
  tree_DataMuEG.Loop("DataMuEG");
  delete tree_DataMuEG;


  TreeReader * tree_DataEG = new TreeReader(tree, "DataEG");
  tree_DataEG.Loop("DataEG");
  delete tree_DataEG;
  
  
  
  TreeReader * tree_DataEGZenriched = new TreeReader(tree, "DataEGZenriched");
  tree_DataEGZenriched.Loop("DataEGZenriched");
  delete tree_DataEGZenriched;
  
  
  TreeReader * tree_DataMuEGZenriched = new TreeReader(tree, "DataMuEGZenriched");
  tree_DataMuEGZenriched.Loop("DataMuEGZenriched");
  delete tree_DataMuEGZenriched;
  
  
  TreeReader * tree_DataMuZenriched = new TreeReader(tree, "DataMuZenriched");
  tree_DataMuZenriched.Loop("DataMuZenriched");
  delete tree_DataMuZenriched;
  
  
  TreeReader * tree_tZq = new TreeReader(tree, "tZq");
  tree_tZq.Loop("tZq");
  delete tree_tZq;
  
  
  
  TreeReader * tree_TTZ = new TreeReader(tree, "TTZ");
  tree_TTZ.Loop("TTZ");
  delete tree_TTZ;
  
  
  TreeReader * tree_TTW = new TreeReader(tree, "TTW");
  tree_TTW.Loop("TTW");
  delete tree_TTW;
  
  
  TreeReader * tree_TT = new TreeReader(tree, "TT");
  tree_TT.Loop("TT");
  delete tree_TT;
  
  TreeReader * tree_DYToLL_M10_50 = new TreeReader(tree, "DYToLL_M10-50");
  tree_DYToLL_M10_50.Loop("DYToLL_M10-50");
  delete tree_DYToLL_M10_50;
  
  
  
  TreeReader * tree_Zjets = new TreeReader(tree, "Zjets");
  tree_Zjets.Loop("Zjets");
  delete tree_Zjets;
  
  TreeReader * tree_Wjets = new TreeReader(tree, "Wjets");
  tree_Wjets.Loop("Wjets");
  delete tree_Wjets;
  
  TreeReader * tree_WW = new TreeReader(tree, "WW");
  tree_WW.Loop("WW");
  delete tree_WW;
  
  TreeReader * tree_WZ = new TreeReader(tree, "WZ");
  tree_WZ.Loop("WZ");
  delete tree_WZ;

  TreeReader * tree_ZZ = new TreeReader(tree, "ZZ");
  tree_ZZ.Loop("ZZ");
  delete tree_ZZ;

  
  TreeReader * tree_TbarsChan = new TreeReader(tree, "TbarsChan");
  tree_TbarsChan.Loop("TbarsChan");
  delete tree_TbarsChan;

  
  TreeReader * tree_TsChan = new TreeReader(tree, "TsChan");
  tree_TsChan.Loop("TsChan");
  delete tree_TsChan;

  
  TreeReader * tree_TtChan = new TreeReader(tree, "TtChan");
  tree_TtChan.Loop("TtChan");
  delete tree_TtChan;


  TreeReader * tree_TbartChan = new TreeReader(tree, "TbartChan");
  tree_TbartChan.Loop("TbartChan");
  delete tree_TbartChan;

  
  TreeReader * tree_TtW = new TreeReader(tree, "TtW");
  tree_TtW.Loop("TtW");
  delete tree_TtW;

  
  TreeReader * tree_TbartW = new TreeReader(tree, "TbartW");
  tree_TbartW.Loop("TbartW");
  delete tree_TbartW;

  
  TreeReader * tree_TbartW = new TreeReader(tree, "TbartW");
  tree_TbartW.Loop("TbartW");
  delete tree_TbartW;

  
  
  
}
