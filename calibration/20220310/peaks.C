/***********************************************
* SIDDHARTA-2 Experiment
* Aleksander K.                 2022-07
* Licensed under the Apache License, Version 2.0
***********************************************/

// Macro for SDDs calibration

#include  <TROOT.h>
#include  <Riostream.h>
#include  <stdlib.h>
#include  <TTree.h>
#include  <TBranch.h>
#include  <TLeaf.h>
#include  <TMinuit.h>
#include  <TFitResult.h>
#include  <TFitResultPtr.h>

#include  "fitting_func_new.h"
#include  "calib_sdd.h"

#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

Int_t npeaks;
Double_t findPeaks(Double_t *x, Double_t *par) {
  Double_t result = par[0]+par[1]*x[0] + TMath::Exp(par[2]+par[3]*x[0]); //background
  for (Int_t p=0; p<npeaks; p++) { //gauss
    Double_t norm  = par[3*p+4];
    Double_t mean  = par[3*p+5];
    Double_t sigma = par[3*p+6];
    result += norm*TMath::Gaus(x[0],mean,sigma); //background + gauss
  }
  return result;
}

void peaks(Int_t np=6) {

  npeaks = np;

  busNumber = 1;
  sddNumber = 25;

  // Read ROOT file
  f[0] = new TFile("../../rootfiles/SIDDHARTA2_xray/output/20220310/hist_20220310_1900_0310_2000_xray_25kv_70ua_tube1_p1.root","READ");

  // Take SDDs histograms
  hSDD = (TH1F*)f[0]->Get(Form("bus%d_sdd%d",busNumber,sddNumber));
  hCrosstalk = (TH1F*)f[0]->Get(Form("bus%d_sdd%d_crosstalk",busNumber,sddNumber));

  hSDD_copy = (TH1F*)hSDD->Clone(Form("hADC_bus%d_sdd%d_copy",busNumber,sddNumber));
  hSDD_copy->Rebin(4);
  hSDD_fit = (TH1F*)hSDD_copy->Clone(Form("hADC_bus%d_sdd%d_fit",busNumber,sddNumber));

  hEnergySDD = (TH1F*)hSDD->Clone(Form("hEnergy_bus%d_sdd%d",busNumber,sddNumber));

  //// Search for peaks ////

  Double_t par[100];

  hSDD_copy->SetAxisRange(1700,3700,"X");
  hSDD_fit->SetAxisRange(1700,3700,"X");

  // Use TSpectrum to find the peak candidates
  TSpectrum *s = new TSpectrum(2*npeaks);
  Int_t nfound = s->Search(hSDD_fit,12,"nobackground",0.0065);
  cout<<Form("Found %d candidate peaks to fitn",nfound)<<endl;

  // Estimate background
  TF1 *fBkgnd = new TF1("fBkgnd",fitBkgndFunc,1500,4500,4);
  hSDD_fit->Fit("fBkgnd","qn");

  // Loop on all found peaks. Eliminate peaks at the background level
  par[0] = fBkgnd->GetParameter(0);
  par[1] = fBkgnd->GetParameter(1);
  par[2] = fBkgnd->GetParameter(2);
  par[3] = fBkgnd->GetParameter(3);
  npeaks = 0;
  Double_t *xpeaks = s->GetPositionX();

  // Evaluation of the centroids for each single peak
  for (Int_t p=0; p<nfound; p++) {
    Float_t xp = xpeaks[p];
    Int_t bin = hSDD_copy->GetXaxis()->FindBin(xp);
    Float_t yp = hSDD_copy->GetBinContent(bin);
    if (yp-TMath::Sqrt(yp) < fBkgnd->Eval(xp)) continue;
    par[3*npeaks+4] = yp;
    par[3*npeaks+5] = xp;
    par[3*npeaks+6] = 15;
    npeaks++;
  }

  cout<<Form("Found %d useful peaks to fitn",npeaks)<<endl;

  TF1 *fit = new TF1("fit",findPeaks,1500,4500,4+3*npeaks);
  fit->SetParameters(par);
  fit->SetNpx(1000);
  hSDD_fit->Fit("fit");

  cout<<"chi2/NDf: "<<(fit->GetChisquare())/(fit->GetNDF())<<endl;

//////////////////////////////////////////////////////

  fitFuncTotal = new TF1(Form("fitFuncTotal_%d_%d",busNumber,sddNumber),fitTotalFunc,1500,4500,30);

  fitFuncTotal->SetParameter(0,fit->GetParameter(4));
  fitFuncTotal->SetParameter(1,fit->GetParameter(5));
  fitFuncTotal->SetParameter(4,fit->GetParameter(7));
  fitFuncTotal->SetParameter(5,fit->GetParameter(8));
  fitFuncTotal->SetParameter(8,fit->GetParameter(10));
  fitFuncTotal->SetParameter(9,fit->GetParameter(11));
  fitFuncTotal->SetParameter(12,fit->GetParameter(13));
  fitFuncTotal->SetParameter(13,fit->GetParameter(14));
  fitFuncTotal->SetParameter(16,fit->GetParameter(16));
  fitFuncTotal->SetParameter(17,fit->GetParameter(17));
  fitFuncTotal->SetParameter(20,fit->GetParameter(19));
  fitFuncTotal->SetParameter(21,fit->GetParameter(20));

  fitFuncTotal->SetParLimits(2,0.00325,0.5);
  fitFuncTotal->SetParLimits(6,0.025,0.5);
  fitFuncTotal->SetParLimits(10,0.004,0.5);
  fitFuncTotal->SetParLimits(14,0.021,0.5);
  fitFuncTotal->SetParLimits(18,0.004,0.5);
  fitFuncTotal->SetParLimits(22,0.021,0.5);
  fitFuncTotal->SetParLimits(3,0.78,100);
  fitFuncTotal->SetParLimits(7,0.42,100);
  fitFuncTotal->SetParLimits(11,0.38,100);
  fitFuncTotal->SetParLimits(15,0.21,100);
  fitFuncTotal->SetParLimits(19,0.38,100);
  fitFuncTotal->SetParLimits(23,0.21,100);

  fitFuncTotal->SetParLimits(24,0.01,0.15);
  fitFuncTotal->SetParLimits(25,20,200);

  for (Int_t k = 0; k < 4; k++) {
    fitFuncTotal->SetParameter(k+26,fit->GetParameter(k));
  }

  fitFuncTotal->SetNpx(1000);

  hSDD_fit->Fit(fitFuncTotal,"","",1500,4500);

  //////////////////////////////

  hResidual = new TH1F("hResidual","",2500,0,10000);
  hSDD_copy->SetAxisRange(0,10000,"X");

  for (Int_t i=1; i<2501; i++) {

    Double_t sg = TMath::Sqrt( fitFuncTotal->GetParameter(24)*3.71*i + TMath::Power(fitFuncTotal->GetParameter(25)/2.35,2) );
    Double_t diff = hSDD_copy->GetBinContent(i) - fitFuncTotal->Eval(hSDD_copy->GetBinCenter(i));
    hResidual->SetBinContent(i,diff/sg);
  }

  //////////////////////////////

  // Parameters from the global function

  fitFuncGaussCuKalpha = new TF1(Form("fitFuncGaussCuKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussCuKalpha->SetParameter(0,fitFuncTotal->GetParameter(0));
  fitFuncGaussCuKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncGaussCuKalpha->SetParameter(2,fitFuncTotal->GetParameter(1));
  fitFuncGaussCuKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncGaussCuKalpha->SetNpx(1000);

  fitFuncTailCuKalpha = new TF1(Form("fitFuncTailCuKalpha_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailCuKalpha->SetParameter(0,(fitFuncTotal->GetParameter(0)*fitFuncTotal->GetParameter(2)));
  fitFuncTailCuKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncTailCuKalpha->SetParameter(2,fitFuncTotal->GetParameter(1));
  fitFuncTailCuKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncTailCuKalpha->SetParameter(4,fitFuncTotal->GetParameter(3));
  fitFuncTailCuKalpha->SetNpx(1000);

  fitFuncGaussCuKbeta = new TF1(Form("fitFuncGaussCuKbeta_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussCuKbeta->SetParameter(0,fitFuncTotal->GetParameter(4));
  fitFuncGaussCuKbeta->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncGaussCuKbeta->SetParameter(2,fitFuncTotal->GetParameter(5));
  fitFuncGaussCuKbeta->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncGaussCuKbeta->SetNpx(1000);

  fitFuncTailCuKbeta = new TF1(Form("fitFuncTailCuKbeta_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailCuKbeta->SetParameter(0,(fitFuncTotal->GetParameter(4)*fitFuncTotal->GetParameter(6)));
  fitFuncTailCuKbeta->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncTailCuKbeta->SetParameter(2,fitFuncTotal->GetParameter(5));
  fitFuncTailCuKbeta->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncTailCuKbeta->SetParameter(4,fitFuncTotal->GetParameter(7));
  fitFuncTailCuKbeta->SetNpx(1000);

  /////

  fitFuncGaussTiKalpha = new TF1(Form("fitFuncGaussTiKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussTiKalpha->SetParameter(0,fitFuncTotal->GetParameter(8));
  fitFuncGaussTiKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncGaussTiKalpha->SetParameter(2,fitFuncTotal->GetParameter(9));
  fitFuncGaussTiKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncGaussTiKalpha->SetNpx(1000);

  fitFuncTailTiKalpha = new TF1(Form("fitFuncTailTiKalpha_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailTiKalpha->SetParameter(0,(fitFuncTotal->GetParameter(8)*fitFuncTotal->GetParameter(10)));
  fitFuncTailTiKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncTailTiKalpha->SetParameter(2,fitFuncTotal->GetParameter(9));
  fitFuncTailTiKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncTailTiKalpha->SetParameter(4,fitFuncTotal->GetParameter(11));
  fitFuncTailTiKalpha->SetNpx(1000);

  fitFuncGaussTiKbeta = new TF1(Form("fitFuncGaussTiKbeta_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussTiKbeta->SetParameter(0,fitFuncTotal->GetParameter(12));
  fitFuncGaussTiKbeta->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncGaussTiKbeta->SetParameter(2,fitFuncTotal->GetParameter(13));
  fitFuncGaussTiKbeta->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncGaussTiKbeta->SetNpx(1000);

  fitFuncTailTiKbeta = new TF1(Form("fitFuncTailTiKbeta_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailTiKbeta->SetParameter(0,(fitFuncTotal->GetParameter(12)*fitFuncTotal->GetParameter(14)));
  fitFuncTailTiKbeta->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncTailTiKbeta->SetParameter(2,fitFuncTotal->GetParameter(13));
  fitFuncTailTiKbeta->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncTailTiKbeta->SetParameter(4,fitFuncTotal->GetParameter(15));
  fitFuncTailTiKbeta->SetNpx(1000);

  fitFuncGaussMnKalpha = new TF1(Form("fitFuncGaussMnKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussMnKalpha->SetParameter(0,fitFuncTotal->GetParameter(16));
  fitFuncGaussMnKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncGaussMnKalpha->SetParameter(2,fitFuncTotal->GetParameter(17));
  fitFuncGaussMnKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncGaussMnKalpha->SetNpx(1000);

  fitFuncGaussFeKalpha = new TF1(Form("fitFuncGaussFeKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussFeKalpha->SetParameter(0,fitFuncTotal->GetParameter(20));
  fitFuncGaussFeKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncGaussFeKalpha->SetParameter(2,fitFuncTotal->GetParameter(21));
  fitFuncGaussFeKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncGaussFeKalpha->SetNpx(1000);

  fitFuncTailMnKalpha = new TF1(Form("fitFuncTailMnKalpha_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailMnKalpha->SetParameter(0,fitFuncTotal->GetParameter(16)*fitFuncTotal->GetParameter(18));
  fitFuncTailMnKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncTailMnKalpha->SetParameter(2,fitFuncTotal->GetParameter(17));
  fitFuncTailMnKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncTailMnKalpha->SetParameter(4,fitFuncTotal->GetParameter(19));
  fitFuncTailMnKalpha->SetNpx(1000);

  fitFuncTailFeKalpha = new TF1(Form("fitFuncTailFeKalpha_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailFeKalpha->SetParameter(0,fitFuncTotal->GetParameter(20)*fitFuncTotal->GetParameter(22));
  fitFuncTailFeKalpha->SetParameter(1,fitFuncTotal->GetParameter(24));
  fitFuncTailFeKalpha->SetParameter(2,fitFuncTotal->GetParameter(21));
  fitFuncTailFeKalpha->SetParameter(3,fitFuncTotal->GetParameter(25));
  fitFuncTailFeKalpha->SetParameter(4,fitFuncTotal->GetParameter(23));
  fitFuncTailFeKalpha->SetNpx(1000);

  /////

  fitFuncBkgrnd = new TF1(Form("fitFuncBkgrnd_%d_%d",busNumber,sddNumber),fitBkgndFunc,1500,4500,4);
  for(int l=0; l<4; l++) {
    fitFuncBkgrnd->SetParameter(l,fitFuncTotal->GetParameter(l+26));
  }

/////////////////////////////////////////HISTOGRAMS/////////////////////////////////////////

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1,0);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextFont(42);

  //
  myCanvas[3] = new TCanvas;
  myCanvas[3]->SetLogy();
  myCanvas[3]->SetGrid();

  hSDD_copy->SetTitle(Form("BUS: %d, SDD: %d",busNumber,sddNumber));
  hSDD_copy->GetXaxis()->SetTitle("ADC [channel]");
  hSDD_copy->GetXaxis()->SetTitleSize(0.05);
  hSDD_copy->GetXaxis()->SetTitleOffset(1.);
  hSDD_copy->GetXaxis()->SetLabelSize(0.04);
  hSDD_copy->SetAxisRange(1700,3700,"X");
  hSDD_copy->GetYaxis()->SetTitle("counts / 4 channel");
  hSDD_copy->GetYaxis()->SetTitleSize(0.05);
  hSDD_copy->GetYaxis()->SetTitleOffset(1.);
  hSDD_copy->GetYaxis()->SetLabelSize(0.04);

  hSDD_copy->SetLineWidth(1);
  hSDD_copy->SetLineColor(1);
  hSDD_copy->Draw();

  fitFuncTotal->SetLineColor(2);
  fitFuncTotal->SetLineWidth(2);
  fitFuncTotal->Draw("same");

  fit->SetLineColor(4);
  fit->SetLineWidth(2);
  //fit->Draw("same");

  fitFuncBkgrnd->SetLineColor(7);
  fitFuncBkgrnd->SetLineWidth(1);
  fitFuncBkgrnd->Draw("same");

  fitFuncGaussTiKalpha->SetLineColor(51);
  fitFuncGaussTiKalpha->SetLineWidth(1);
  fitFuncGaussTiKalpha->Draw("same");

  fitFuncGaussTiKbeta->SetLineColor(51);
  fitFuncGaussTiKbeta->SetLineWidth(1);
  fitFuncGaussTiKbeta->Draw("same");

  fitFuncTailTiKalpha->SetLineColor(8);
  fitFuncTailTiKalpha->SetLineWidth(1);
  fitFuncTailTiKalpha->Draw("same");

  fitFuncTailTiKbeta->SetLineColor(8);
  fitFuncTailTiKbeta->SetLineWidth(1);
  fitFuncTailTiKbeta->Draw("same");

  fitFuncGaussCuKalpha->SetLineColor(51);
  fitFuncGaussCuKalpha->SetLineWidth(1);
  fitFuncGaussCuKalpha->Draw("same");

  fitFuncGaussCuKbeta->SetLineColor(51);
  fitFuncGaussCuKbeta->SetLineWidth(1);
  fitFuncGaussCuKbeta->Draw("same");

  fitFuncTailCuKalpha->SetLineColor(8);
  fitFuncTailCuKalpha->SetLineWidth(1);
  fitFuncTailCuKalpha->Draw("same");

  fitFuncTailCuKbeta->SetLineColor(8);
  fitFuncTailCuKbeta->SetLineWidth(1);
  fitFuncTailCuKbeta->Draw("same");

  fitFuncGaussMnKalpha->SetLineColor(51);
  fitFuncGaussMnKalpha->SetLineWidth(1);
  fitFuncGaussMnKalpha->Draw("same");

  fitFuncGaussFeKalpha->SetLineColor(51);
  fitFuncGaussFeKalpha->SetLineWidth(1);
  fitFuncGaussFeKalpha->Draw("same");

  fitFuncTailMnKalpha->SetLineColor(8);
  fitFuncTailMnKalpha->SetLineWidth(1);
  fitFuncTailMnKalpha->Draw("same");

  fitFuncTailFeKalpha->SetLineColor(8);
  fitFuncTailFeKalpha->SetLineWidth(1);
  fitFuncTailFeKalpha->Draw("same");

  myLegend[3] = new TLegend(0.285, 0.695, 0.565, 0.900);
  myLegend[3]->SetFillStyle(1001); myLegend[3]->SetLineColor(1); myLegend[3]->SetFillColor(0); myLegend[3]->SetTextSize(0.037);
  myLegend[3]->AddEntry(hSDD, "Data", "l");
  myLegend[3]->AddEntry(fitFuncTotal, "Global fit function", "l");
  myLegend[3]->AddEntry(fitFuncGaussTiKalpha, "Gaussian function", "l");
  myLegend[3]->AddEntry(fitFuncTailTiKalpha, "Tail function", "l");
  myLegend[3]->AddEntry(fitFuncBkgrnd, "Background", "l");
  myLegend[3]->Draw("same");

  myCanvas[3]->Print(Form("plots/hADC_Ylog_bus%d_sdd%d.png",busNumber,sddNumber), "png");

  hResidual->SetMarkerStyle(20);
  hResidual->SetMarkerSize(0.5);
  //hResidual->GetXaxis()->SetRangeUser(1700,3700);

  myCanvas[7] = new TCanvas();

  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  pad1->SetBottomMargin(0.008);
  pad1->SetBorderMode(0);
  pad1->SetLogy();
  pad1->SetGrid();
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.33);
  //pad2->SetBorderMode(0);
  pad2->SetGrid();
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  hSDD_copy->GetYaxis()->SetTitleOffset(0.8);
  hSDD_copy->GetYaxis()->SetTitleSize(0.07);
  hSDD_copy->GetYaxis()->SetLabelSize(0.06);
  hSDD_copy->GetYaxis()->SetTitleOffset(0.67);

  hSDD_copy->Draw();
  fitFuncTotal->Draw("same");
  fitFuncBkgrnd->Draw("same");
  fitFuncGaussTiKalpha->Draw("same");
  fitFuncGaussTiKbeta->Draw("same");
  fitFuncTailTiKalpha->Draw("same");
  fitFuncTailTiKbeta->Draw("same");
  fitFuncGaussCuKalpha->Draw("same");
  fitFuncGaussCuKbeta->Draw("same");
  fitFuncTailCuKalpha->Draw("same");
  fitFuncTailCuKbeta->Draw("same");
  fitFuncGaussMnKalpha->Draw("same");
  fitFuncTailMnKalpha->Draw("same");
  fitFuncGaussFeKalpha->Draw("same");
  fitFuncTailFeKalpha->Draw("same");
  myLegend[3]->Draw("same");

  pad2->cd();
  hResidual->GetXaxis()->SetTitle("ADC [channel]");
  hResidual->GetXaxis()->SetTitleSize(0.15);
  hResidual->GetXaxis()->SetTitleOffset(1.);
  hResidual->GetXaxis()->SetLabelSize(0.12);
  hResidual->SetAxisRange(1700,3700,"X");
  hResidual->GetYaxis()->SetTitle("residual / #sigma");
  hResidual->GetYaxis()->SetTitleSize(0.15);
  hResidual->GetYaxis()->SetTitleOffset(0.33);
  hResidual->GetYaxis()->SetLabelSize(0.12);
  hResidual->GetYaxis()->SetNdivisions(5,5,0, kTRUE);
  hResidual->Draw("same p");

  myCanvas[7]->Print(Form("plots/hADC_Ylog_res_bus%d_sdd%d.png",busNumber,sddNumber), "png");

}
