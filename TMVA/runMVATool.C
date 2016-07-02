{
  gROOT->ProcessLine(".L theMVAtool.C+");
  
   theMVAtool t;
   
  t.doTraining();
  t.doTraining("mumu");
  t.doTraining("ee");
  
  
  t.doReading();
  t.doReading(100, "mumu");
  t.doReading(100,"ee");
  
  
  
}
