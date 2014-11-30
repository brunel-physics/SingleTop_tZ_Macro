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
  
  TreeReader * tree_FCNCzut = new TreeReader(tree, "FCNCzut", datalist);
  tree_FCNCzut.Loop("FCNCzut");
  delete tree_FCNCzut;   
   
}
