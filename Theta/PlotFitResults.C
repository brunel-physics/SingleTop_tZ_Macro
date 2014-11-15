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

bool postfit = true;

void PlotFitResults(bool usePostFit, TString distrib, TString inputfilename, TString outputplotname){
  
  
      
  Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;
  
  bool setlogy = 0;
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);
  
  
  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // ROOT . gStyle . SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(0); // kWhite
  gStyle->SetPadGridX(0); //false
  gStyle->SetPadGridY(0); //false
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  
  // For the histo:rebin
  // ROOT . gStyle . SetHistFillColor(1);
  // ROOT . gStyle . SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // ROOT . gStyle . SetLegoInnerR(Float_t rad = 0.5);
  // ROOT . gStyle . SetNumberContours(Int_t number = 20);
  
  gStyle->SetEndErrorSize(2);
  //ROOT . gStyle . SetErrorMarker(20);   /// I COMMENTED THIS OUT
  //ROOT . gStyle . SetErrorX(0.);
  
  //ROOT . gStyle . SetMarkerStyle(20);
  
  
  //For the fit/function:
  gStyle->SetOptFit(1011);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);
  
  //For the date:
  gStyle->SetOptDate(0);
  // ROOT . gStyle . SetDateX(Float_t x = 0.01);
  // ROOT . gStyle . SetDateY(Float_t y = 0.01);
  
  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(0); // kWhite
  gStyle->SetStatFont(42);
  //ROOT . gStyle . SetStatFontSize(0.025);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // ROOT . gStyle . SetStatStyle(Style_t style = 1001);
  // ROOT . gStyle . SetStatX(Float_t x = 0);
  // ROOT . gStyle . SetStatY(Float_t y = 0);
  
  // Margins:
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  //ROOT . gStyle . SetPadRightMargin(0.12);
  gStyle->SetPadRightMargin(0.03);
  
  // For the Global title:
  
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // ROOT . gStyle . SetTitleH(0); // Set the height of the title box
  // ROOT . gStyle . SetTitleW(0); // Set the width of the title box
  // ROOT . gStyle . SetTitleX(0); // Set the position of the title box
  // ROOT . gStyle . SetTitleY(0.985); // Set the position of the title box
  // ROOT . gStyle . SetTitleStyle(Style_t style = 1001);
  // ROOT . gStyle . SetTitleBorderSize(2);
  
  // For the axis titles:
  
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // ROOT . gStyle . SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // ROOT . gStyle . SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // ROOT . gStyle . SetTitleOffset(1.1, "Y"); // Another way to set the Offset
  
  // For the axis labels:
  
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  
  // For the axis:
  
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(1); // kTRUE
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  
  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);
  
  // Postscript options:
  gStyle->SetPaperSize(20.,20.);
  // ROOT . gStyle . SetLineScalePS(Float_t scale = 3);
  // ROOT . gStyle . SetLineStyleString(Int_t i, const char* text);
  // ROOT . gStyle . SetHeaderPS(const char* header);
  // ROOT . gStyle . SetTitlePS(const char* pstitle);
  
  // ROOT . gStyle . SetBarOffset(Float_t baroff = 0.5);
  // ROOT . gStyle . SetBarWidth(Float_t barwidth = 0.5);
  // ROOT . gStyle . SetPaintTextFormat(const char* format = "g");
  // ROOT . gStyle . SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // ROOT . gStyle . SetTimeOffset(Double_t toffset);
  // ROOT . gStyle . SetHistMinimumZero(kTRUE);
  
  
  //TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
  TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
  c1->SetBottomMargin(0.3);
  c1->SetLogy(setlogy);
  c1->cd();
   
  TFile * inputfile_fit ;
  if(usePostFit) inputfile_fit = new TFile("histos-mle_tZq.root");
  else  inputfile_fit = new TFile(inputfilename);
  
  TH1F * hist_NJetsNBjets__TT  = (TH1F*)inputfile_fit->Get( (distrib+"__TT").Data() )->Clone();
  TH1F * hist_NJetsNBjets__TTZ  = (TH1F*)inputfile_fit->Get( (distrib+"__TTZ").Data() )->Clone();
  TH1F * hist_NJetsNBjets__TTW  = (TH1F*)inputfile_fit->Get( (distrib+"__TTW").Data() )->Clone();
  TH1F * hist_NJetsNBjets__Zjets    = (TH1F*)inputfile_fit->Get( (distrib+"__Zjets").Data() )->Clone();
  TH1F * hist_NJetsNBjets__tZq = (TH1F*)inputfile_fit->Get( (distrib+"__tZq").Data() )->Clone();
  TH1F * hist_NJetsNBjets__WZ  = (TH1F*)inputfile_fit->Get( (distrib+"__WZ").Data() )->Clone();
  TH1F * hist_NJetsNBjets__WZHF    = (TH1F*)inputfile_fit->Get( (distrib+"__WZHF").Data() )->Clone();
  
  

  THStack* hs= new THStack();
 
  hist_NJetsNBjets__TTZ->SetFillColor(kRed+1);
  hist_NJetsNBjets__TT->SetFillColor(kRed-7);
  hist_NJetsNBjets__Zjets->SetFillColor(kAzure-2);
  hist_NJetsNBjets__TTW->SetFillColor(kRed+1);
  hist_NJetsNBjets__tZq->SetFillColor(kGreen+2);
  hist_NJetsNBjets__WZ->SetFillColor(13); 
  hist_NJetsNBjets__WZHF->SetFillColor(12); 
  
  hs->Add(hist_NJetsNBjets__tZq);
  hs->Add(hist_NJetsNBjets__Zjets);
  hs->Add(hist_NJetsNBjets__WZ);
  hs->Add(hist_NJetsNBjets__TT);
  hs->Add(hist_NJetsNBjets__TTZ);
  hs->Add(hist_NJetsNBjets__TTW);
  hs->Add(hist_NJetsNBjets__WZHF);
  
  hs->Draw("histo");
  hs->GetXaxis()->SetLabelSize(0.);
  hs->GetYaxis()->SetLabelSize(0.03);
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetTitleSize(0.04);
  if(distrib=="NJetsNBjets") hs->SetMaximum(30000);
  if(distrib=="InvMass") hs->SetMaximum(8000);
  if(distrib=="MVA_all") hs->SetMaximum(150);
  
  TFile * inputfile_data = new TFile(inputfilename);
  TH1F * hist_data	   = (TH1F*)inputfile_data->Get( (distrib+"__DATA").Data() )->Clone();
  
  
  /*cout <<  "hist_data               " << hist_data->GetXaxis()->GetNbins()              << " " << hist_data->GetXaxis()->GetXmin()               << " " << hist_data->GetXaxis()->GetXmax() << endl;
  cout <<  "hist_NJetsNBjets__TTZ  " << hist_NJetsNBjets__TTZ->GetXaxis()->GetNbins() << " " << hist_NJetsNBjets__TTZ->GetXaxis()->GetXmin()  << " " << hist_NJetsNBjets__TTZ->GetXaxis()->GetXmax() << endl;
  cout <<  "hist_NJetsNBjets__Zjets    " << hist_NJetsNBjets__Zjets->GetXaxis()->GetNbins()   << " " << hist_NJetsNBjets__Zjets->GetXaxis()->GetXmin()    << " " << hist_NJetsNBjets__Zjets->GetXaxis()->GetXmax() << endl;
  cout <<  "hist_NJetsNBjets__TTW  " << hist_NJetsNBjets__TTW->GetXaxis()->GetNbins() << " " << hist_NJetsNBjets__TTW->GetXaxis()->GetXmin()  << " " << hist_NJetsNBjets__TTW->GetXaxis()->GetXmax() << endl;
  cout <<  "hist_NJetsNBjets__tZq " << hist_NJetsNBjets__tZq->GetXaxis()->GetNbins()<< " " << hist_NJetsNBjets__tZq->GetXaxis()->GetXmin() << " " << hist_NJetsNBjets__tZq->GetXaxis()->GetXmax() << endl;
  cout <<  "hist_NJetsNBjets__WZ  " << hist_NJetsNBjets__WZ->GetXaxis()->GetNbins() << " " << hist_NJetsNBjets__WZ->GetXaxis()->GetXmin()  << " " << hist_NJetsNBjets__WZ->GetXaxis()->GetXmax() << endl;
  cout <<  "hist_NJetsNBjets__WZHF    " << hist_NJetsNBjets__WZHF->GetXaxis()->GetNbins()   << " " << hist_NJetsNBjets__WZHF->GetXaxis()->GetXmin()    << " " << hist_NJetsNBjets__WZHF->GetXaxis()->GetXmax() << endl;
  */
  hist_data->SetMarkerStyle(20);
  
  hist_data->Draw("epsame");
  
  
  //determine the error band
  //fixme need to be updated
  
  TH1* herrorband = (TH1F*) hist_NJetsNBjets__TTZ->Clone();
  herrorband->Add(hist_NJetsNBjets__TTW);
  herrorband->Add(hist_NJetsNBjets__TT);
  herrorband->Add(hist_NJetsNBjets__Zjets);
  herrorband->Add(hist_NJetsNBjets__tZq);
  herrorband->Add(hist_NJetsNBjets__WZ);
  herrorband->Add(hist_NJetsNBjets__WZHF);
  
  herrorband->SetMarkerStyle(21) ;
  herrorband->SetMarkerSize(1.2) ;  
  TGraphAsymmErrors *thegraph = new TGraphAsymmErrors(herrorband);
  thegraph->SetFillStyle(3005);
  thegraph->SetFillColor(1);
  thegraph->Draw("e2same");
  
  //ratio plot
  
  TH1F * histo_ratio = (TH1F*) hist_data->Clone();
  TH1F * histo_mc    = (TH1F*) hist_NJetsNBjets__TTZ->Clone();
  histo_mc->Add(hist_NJetsNBjets__TTW);
  histo_mc->Add(hist_NJetsNBjets__TT);
  histo_mc->Add(hist_NJetsNBjets__Zjets);
  histo_mc->Add(hist_NJetsNBjets__tZq);
  histo_mc->Add(hist_NJetsNBjets__WZ);
  histo_mc->Add(hist_NJetsNBjets__WZHF);
  histo_ratio->Divide(histo_ratio, histo_mc, 1, 1);
  histo_ratio->GetXaxis()->SetLabelSize(0.07);
  //if(distrib=="NJetsNBjets") histo_ratio->GetXaxis()->SetRange(1, 15);
  histo_ratio->GetXaxis()->SetTitle("BDT output");
  
     
  TGraphAsymmErrors *theRatio = new TGraphAsymmErrors(histo_ratio);
  for (int ierr=1; ierr<=histo_ratio->GetNbinsX(); ierr++) {
    
    double num     = hist_data->GetBinContent(ierr);
    double num_err = hist_data->GetBinError(ierr);
    
    double denom         = herrorband->GetBinContent(ierr);
    double denom_err     = 0;
    double error = pow(
       pow(num_err/denom, 2)+
       pow(num*denom_err/(denom*denom), 2)
    , 0.5);
    
    cout << "***************"  << ierr << endl;
    cout << "num " << num << " pm " << num_err << endl;
    cout << "denom " << denom << " pm " << denom_err << endl;
    theRatio->SetPointEYhigh( ierr-1, error);
    theRatio->SetPointEYlow(  ierr-1, error);
    
    cout << "num/denom " << num/denom <<" PM  " << error<< endl;
    cout << "num_err/denom " << num_err/denom << endl;
    cout << "num*denom_err/(denom*denom) " << num*denom_err/(denom*denom) << endl;
    //cout << "get error check " << theRatio->GetErrorYhigh(ierr-1) << endl;
    //cout << "get value check " << histo_ratio->GetBinContent(ierr) << endl;
  
  }
  
  
  
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2", 0.0, 0.0, 1.0, 1.0);
  canvas_2->SetTopMargin(0.7);
  canvas_2->SetFillColor(0);
  canvas_2->SetFillStyle(0);
  canvas_2->SetGridy(1);
  canvas_2->Draw();
  canvas_2->cd(0);
  histo_ratio->SetTitle("");
  
  histo_ratio->SetMarkerStyle(20);
  histo_ratio->SetMarkerSize(1.2);
  histo_ratio->SetMaximum( 1.5 );
  histo_ratio->SetMinimum(0.5);
  histo_ratio->GetYaxis()->SetTitle("");
  histo_ratio->GetXaxis()->SetLabelSize(0.04);
  histo_ratio->GetYaxis()->SetLabelSize(0.04);
  histo_ratio->GetYaxis()->SetNdivisions(6);

  
  histo_ratio->GetYaxis()->SetTitleSize(0.03);
  histo_ratio->SetMarkerSize(1.2);
  //histo_ratio->GetYaxis()->SetNdivisions(5);
  //ratio.Draw("e")
  
  histo_ratio->SetMinimum(0.49);
  histo_ratio->SetMaximum(1.51);
  histo_ratio->SetLineColor(0);
  histo_ratio->SetMarkerColor(0);
  
  histo_ratio->Draw("ep");
  
  
  
  
  theRatio->SetMarkerSize(1.2);
  theRatio->Draw("ep");
   
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(31); 
  latex->DrawLatex(0.45, 0.95, "CMS Preliminary");
  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.04);
  latex2->SetTextAlign(31); 
  latex2->DrawLatex(0.87, 0.95, "20 fb^{-1} at #sqrt{s} = 8 TeV");

  TString info_data = "e#mu channel";; 
  
  
  //text2 = new TLatex(0.15,0.8, info_data);
  TLatex * text2 = new TLatex(0.80,0.98, info_data);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.60);
  text2->SetY(0.92);
  //text2->SetLineWidth(2);
  text2->SetTextFont(42);
  text2->SetTextSize(0.0610687);
  //    text2->SetTextSizePixels(24);// dflt=28
  text2->Draw();
  
  
  
  
  
  TLegend* qw = new TLegend(.70,.60,.95,.80);
  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  qw->AddEntry(hist_data,                 "Data","ep");
  qw->AddEntry(hist_NJetsNBjets__TT,    "TT"   ,"f");
  qw->AddEntry(hist_NJetsNBjets__TTZ,    "TTZ"   ,"f");
  qw->AddEntry(hist_NJetsNBjets__TTW,    "TTW"   ,"f");
  qw->AddEntry(hist_NJetsNBjets__Zjets,      "Z/#gamma^{*}"   ,"f");
  qw->AddEntry(hist_NJetsNBjets__tZq,   "tZq"   ,"f");
  qw->AddEntry(hist_NJetsNBjets__WZ,    "VV"   ,"f");
  qw->AddEntry(hist_NJetsNBjets__WZHF,      "WZHF"   ,"f");
  
  qw->Draw();
  

  c1->SaveAs( (outputplotname+".gif").Data());
 
}


void PlotFitResults(){

   
  //PlotFitResults(postfit, "NJetsNBjets", "RootFiles/Jun22_Jet30_CSVT_TopPtLepSysMll_DD.root", "PreFit_NJetsNBjets_CSVT");
  PlotFitResults(postfit, "MVA_all", "../TMVA/TemplateRootFiles/MVA_all_theta.root", "PostFitMVA");
  
  //PlotFitResults(true, "InvMass");

}
