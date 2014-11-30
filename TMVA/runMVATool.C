{
  gROOT->ProcessLine(".L theMVAtool.C+");
  
   theMVAtool t;
   
  t.doTraining();
  t.doTraining("mumumu");
  t.doTraining("mumue");
  t.doTraining("eemu");
  t.doTraining("eee");
  
  
  t.doReading();
  t.doReading(100, "mumumu");
  t.doReading(100,"mumue");
  t.doReading(100,"eemu");
  t.doReading(100,"eee");
  
  
  
}
