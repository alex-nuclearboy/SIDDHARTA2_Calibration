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

void calib_sdd_new() {

  busNumber = 1;
  sddNumber = 25;

  beginWindow[0] = 1925;          //Ti K_alpha
  endWindow[0] = 2060;
  beginWindow[1] = endWindow[0];  //Ti K_beta
  endWindow[1] = 2200;

  beginWindow[2] = 2400;          //Mn K_alpha
  endWindow[2] = 2525;

  beginWindow[3] = endWindow[2];  //Fe K_alpha
  endWindow[3] = 2675;

  beginWindow[4] = 3070;          //Cu K_alpha
  endWindow[4] = 3170;

  beginWindow[5] = 3345;          //Cu K_beta
  endWindow[5] = 3445;

  //Read ROOT file
  f[0] = new TFile("../../rootfiles/SIDDHARTA2_xray/output/20220310/hist_20220310_1900_0310_2000_xray_25kv_70ua_tube1_p1.root","READ");

  hSDD = (TH1F*)f[0]->Get(Form("bus%d_sdd%d",busNumber,sddNumber));  //take SDD histograms
  hCrosstalk = (TH1F*)f[0]->Get(Form("bus%d_sdd%d_crosstalk",busNumber,sddNumber));  //take crosstalk histograms

  //hSDD->Rebin(4);
  hSDD_copy = (TH1F*)hSDD->Clone(Form("hADC_bus%d_sdd%d_copy",busNumber,sddNumber));
  hSDD_copy->Rebin(4);
  hSDD_fit = (TH1F*)hSDD_copy->Clone(Form("hADC_bus%d_sdd%d_fit",busNumber,sddNumber));

  hEnergySDD = (TH1F*)hSDD->Clone(Form("hEnergy_bus%d_sdd%d",busNumber,sddNumber));

  // Find peaks
  for(Int_t i=0; i<6; i++){
    //hSDD_copy->GetXaxis()->SetRange(beginWindow[0],endWindow[0]);
    hSDD_copy->SetAxisRange(beginWindow[i],endWindow[i],"X");
    binPeak[i]=hSDD_copy->GetXaxis()->GetBinCenter(hSDD_copy->GetMaximumBin());
    peakValue[i]=hSDD_copy->GetMaximum();
  }

  //// Fitting ////

  // Background
  fitFuncBkgrndStart = new TF1(Form("fitFuncBkgrndStart_%d_%d",busNumber,sddNumber), fitBkgndFunc, 1000,5000,4);
  hSDD_fit->Fit(Form("fitFuncBkgrndStart_%d_%d",busNumber,sddNumber),"","",1500,4500);
/*
  // Ti //

  // K_alpha
  fitFuncGaussTailTiKalphaStart = new TF1(Form("fitFuncGaussTailTiKalphaStart_%d_%d",busNumber,sddNumber), fitGaussTailFunc, 1900,2100, 6);
  fitFuncGaussTailTiKalphaStart->SetParName(0,"Gain");
  fitFuncGaussTailTiKalphaStart->SetParName(1,"Fano");
  fitFuncGaussTailTiKalphaStart->SetParName(2,"E");
  fitFuncGaussTailTiKalphaStart->SetParName(3,"Noise");
  fitFuncGaussTailTiKalphaStart->SetParName(4,"Fraction");
  fitFuncGaussTailTiKalphaStart->SetParName(5,"Tail width");
  fitFuncGaussTailTiKalphaStart->SetParameters(2.6e5,0.114,binPeak[0],111,9e-4,1.52);
  fitFuncGaussTailTiKalphaStart->SetParLimits(1,0.10,0.17);
  fitFuncGaussTailTiKalphaStart->SetParLimits(3,25,250);
  fitFuncGaussTailTiKalphaStart->SetParLimits(4,-7.6e-3,0.25);
  fitFuncGaussTailTiKalphaStart->SetParLimits(5,0.5,100);

  hSDD_fit->Fit(Form("fitFuncGaussTailTiKalphaStart_%d_%d",busNumber,sddNumber),"","",binPeak[0]-35,binPeak[0]+35);

  // K_beta
  fitFuncGaussTailTiKbetaStart = new TF1(Form("fitFuncGaussTailTiKbetaStart_%d_%d",busNumber,sddNumber), fitGaussTailFunc, 2100,2200, 6);
  fitFuncGaussTailTiKbetaStart->SetParName(0,"Gain");
  fitFuncGaussTailTiKbetaStart->SetParName(1,"Fano");
  fitFuncGaussTailTiKbetaStart->SetParName(2,"E");
  fitFuncGaussTailTiKbetaStart->SetParName(3,"Noise");
  fitFuncGaussTailTiKbetaStart->SetParName(4,"Fraction");
  fitFuncGaussTailTiKbetaStart->SetParName(5,"Tail width");
  fitFuncGaussTailTiKbetaStart->SetParameters(4e4,0.114,binPeak[1],210,7e-4,1.95);
  fitFuncGaussTailTiKbetaStart->SetParLimits(1,0.10,0.15);
  fitFuncGaussTailTiKbetaStart->SetParLimits(3,100,250);
  fitFuncGaussTailTiKbetaStart->SetParLimits(4,-5.26e-3,0.25);
  fitFuncGaussTailTiKbetaStart->SetParLimits(5,1.5,100);

  hSDD_fit->Fit(Form("fitFuncGaussTailTiKbetaStart_%d_%d",busNumber,sddNumber),"","",binPeak[1]-25,binPeak[1]+25);

  // Mn & Fe //
  fitFuncGaussMnKalphaStart = new TF1(Form("fitFuncGaussMnKalphaStart_%d_%d",busNumber,sddNumber), fitGaussFunc, 2400,2525, 4);
  fitFuncGaussMnKalphaStart->SetParameters(1000,0.114,binPeak[2],111);
  fitFuncGaussMnKalphaStart->SetParLimits(1,0.1,0.15);
  fitFuncGaussMnKalphaStart->SetParLimits(3,100,250);

  hSDD_fit->Fit(Form("fitFuncGaussMnKalphaStart_%d_%d",busNumber,sddNumber),"","",binPeak[2]-40,binPeak[2]+40);

  fitFuncGaussFeKalphaStart = new TF1(Form("fitFuncGaussFeKalphaStart_%d_%d",busNumber,sddNumber), fitGaussFunc, 2525,2675, 4);
  fitFuncGaussFeKalphaStart->SetParameters(100,0.114,binPeak[3],111);
  fitFuncGaussFeKalphaStart->SetParLimits(1,0.1,0.15);
  fitFuncGaussFeKalphaStart->SetParLimits(3,100,250);


  hSDD_fit->Fit(Form("fitFuncGaussFeKalphaStart_%d_%d",busNumber,sddNumber),"","",binPeak[3]-50,binPeak[3]+50);
*/
  // Cu //

  // K_alpha
  fitFuncGaussTailCuKalphaStart = new TF1(Form("fitFuncGaussTailCuKalphaStart_%d_%d",busNumber,sddNumber), fitGaussTailFunc, 3025,3225, 6);
  fitFuncGaussTailCuKalphaStart->SetParName(0,"Gain");
  fitFuncGaussTailCuKalphaStart->SetParName(1,"Fano");
  fitFuncGaussTailCuKalphaStart->SetParName(2,"E");
  fitFuncGaussTailCuKalphaStart->SetParName(3,"Noise");
  fitFuncGaussTailCuKalphaStart->SetParName(4,"Fraction");
  fitFuncGaussTailCuKalphaStart->SetParName(5,"Tail width");
  fitFuncGaussTailCuKalphaStart->SetParameters(2.6e5,0.114,binPeak[4],111,8e-4,1.52);
  fitFuncGaussTailCuKalphaStart->SetParLimits(1,0.10,0.15);
  fitFuncGaussTailCuKalphaStart->SetParLimits(3,100,250);
  fitFuncGaussTailCuKalphaStart->SetParLimits(4,7.6e-4,0.25);
  fitFuncGaussTailCuKalphaStart->SetParLimits(5,1.5,100);

  hSDD_fit->Fit(Form("fitFuncGaussTailCuKalphaStart_%d_%d",busNumber,sddNumber),"","",binPeak[4]-40,binPeak[4]+40);

  // K_beta
  fitFuncGaussTailCuKbetaStart = new TF1(Form("fitFuncGaussTailCuKbetaStart_%d_%d",busNumber,sddNumber), fitGaussTailFunc, 3300,3500, 6);
  fitFuncGaussTailCuKbetaStart->SetParName(0,"Gain");
  fitFuncGaussTailCuKbetaStart->SetParName(1,"Fano");
  fitFuncGaussTailCuKbetaStart->SetParName(2,"E");
  fitFuncGaussTailCuKbetaStart->SetParName(3,"Noise");
  fitFuncGaussTailCuKbetaStart->SetParName(4,"Fraction");
  fitFuncGaussTailCuKbetaStart->SetParName(5,"Tail width");
  fitFuncGaussTailCuKbetaStart->SetParameters(4e4,0.114,binPeak[5],210,7e-4,1.95);
  fitFuncGaussTailCuKbetaStart->SetParLimits(1,0.10,0.15);
  fitFuncGaussTailCuKbetaStart->SetParLimits(3,100,250);
  fitFuncGaussTailCuKbetaStart->SetParLimits(4,5.26e-4,0.25);
  fitFuncGaussTailCuKbetaStart->SetParLimits(5,1.5,100);

  hSDD_fit->Fit(Form("fitFuncGaussTailCuKbetaStart_%d_%d",busNumber,sddNumber),"","",binPeak[5]-40,binPeak[5]+40);

  //// Global fit function ////

  fitFuncTotal = new TF1(Form("fitFuncTotal_%d_%d",busNumber,sddNumber),fitTotalFunc,1500,4500,14);
/*
  fitFuncTotal->SetParName(0,"Ti_a: Gain");
  fitFuncTotal->SetParName(1,"Fano");
  fitFuncTotal->SetParName(2,"E");
  fitFuncTotal->SetParName(3,"Noise");
  fitFuncTotal->SetParName(4,"Fraction");
  fitFuncTotal->SetParName(5,"Tail width");
  fitFuncTotal->SetParName(6,"Ti_b: Gain");
  fitFuncTotal->SetParName(7,"Fano");
  fitFuncTotal->SetParName(8,"E");
  fitFuncTotal->SetParName(9,"Noise");
  fitFuncTotal->SetParName(10,"Fraction");
  fitFuncTotal->SetParName(11,"Tail width");
  fitFuncTotal->SetParName(12,"Mn: Gain");
  fitFuncTotal->SetParName(13,"Fano");
  fitFuncTotal->SetParName(14,"E");
  fitFuncTotal->SetParName(15,"Noise");
  fitFuncTotal->SetParName(16,"Fe: Gain");
  fitFuncTotal->SetParName(17,"Fano");
  fitFuncTotal->SetParName(18,"E");
  fitFuncTotal->SetParName(19,"Noise");
  fitFuncTotal->SetParName(20,"Cu_a: Gain");
  fitFuncTotal->SetParName(21,"Fano");
  fitFuncTotal->SetParName(22,"E");
  fitFuncTotal->SetParName(23,"Noise");
  fitFuncTotal->SetParName(24,"Fraction");
  fitFuncTotal->SetParName(25,"Tail width");
  fitFuncTotal->SetParName(26,"Cu_b: Gain");
  fitFuncTotal->SetParName(27,"Fano");
  fitFuncTotal->SetParName(28,"E");
  fitFuncTotal->SetParName(29,"Noise");
  fitFuncTotal->SetParName(30,"Fraction");
  fitFuncTotal->SetParName(31,"Tail width");
  fitFuncTotal->SetParName(32,"a");
  fitFuncTotal->SetParName(33,"b");
  fitFuncTotal->SetParName(34,"c");
  fitFuncTotal->SetParName(35,"d");
*/
  //take parameters from previous fits

  fitFuncTotal->SetParameter(0,4000);
  fitFuncTotal->SetParameter(1,binPeak[4]);
  fitFuncTotal->SetParLimits(2,-5.e-3,0.3);
  fitFuncTotal->SetParLimits(3,0.5,10);
  fitFuncTotal->SetParameter(4,1700);
  fitFuncTotal->SetParameter(5,binPeak[5]);
  fitFuncTotal->SetParLimits(6,5.e-3,2.5);
  fitFuncTotal->SetParLimits(7,0.5,100);
  fitFuncTotal->SetParLimits(8,0.01,0.15);
  fitFuncTotal->SetParLimits(9,50,250);
  fitFuncTotal->SetParameter(10,-1.5);
  fitFuncTotal->SetParameter(11,1.0E-4);
  fitFuncTotal->SetParameter(12,1.73);
  fitFuncTotal->SetParameter(13,-4.8e-8);

  hSDD_fit->Fit(fitFuncTotal,"","",1500,4500);

  Double_t chi2 = (fitFuncTotal->GetChisquare())/(fitFuncTotal->GetNDF());

  ////////////////////////////////////////
/*
  // Parameters from the global function

  fitFuncGaussTiKalpha = new TF1(Form("fitFuncGaussTiKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  for(int l=0; l<4; l++) {
    fitFuncGaussTiKalpha->SetParameter(l,fitFuncTotal->GetParameter(l));
  }

  fitFuncTailTiKalpha = new TF1(Form("fitFuncTailTiKalpha_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailTiKalpha->SetParameter(0,(fitFuncTotal->GetParameter(0)*fitFuncTotal->GetParameter(4)));
  fitFuncTailTiKalpha->SetParameter(1,fitFuncTotal->GetParameter(1));
  fitFuncTailTiKalpha->SetParameter(2,fitFuncTotal->GetParameter(2));
  fitFuncTailTiKalpha->SetParameter(3,fitFuncTotal->GetParameter(3));
  fitFuncTailTiKalpha->SetParameter(4,fitFuncTotal->GetParameter(5));

  fitFuncGaussTiKbeta = new TF1(Form("fitFuncGaussTiKbeta_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  for(int l=0; l<4; l++) {
    fitFuncGaussTiKbeta->SetParameter(l,fitFuncTotal->GetParameter(l+6));
  }

  fitFuncTailTiKbeta = new TF1(Form("fitFuncTailTiKbeta_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailTiKbeta->SetParameter(0,(fitFuncTotal->GetParameter(6)*fitFuncTotal->GetParameter(10)));
  fitFuncTailTiKbeta->SetParameter(1,fitFuncTotal->GetParameter(7));
  fitFuncTailTiKbeta->SetParameter(2,fitFuncTotal->GetParameter(8));
  fitFuncTailTiKbeta->SetParameter(3,fitFuncTotal->GetParameter(9));
  fitFuncTailTiKbeta->SetParameter(4,fitFuncTotal->GetParameter(11));

  fitFuncGaussMnKalpha = new TF1(Form("fitFuncGaussMnKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  for(int l=0; l<4; l++) {
    fitFuncGaussMnKalpha->SetParameter(l,fitFuncTotal->GetParameter(l+12));
  }

  fitFuncGaussFeKalpha = new TF1(Form("fitFuncGaussFeKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  for(int l=0; l<4; l++) {
    fitFuncGaussFeKalpha->SetParameter(l,fitFuncTotal->GetParameter(l+16));
  }
*/
  fitFuncGaussCuKalpha = new TF1(Form("fitFuncGaussCuKalpha_%d_%d",busNumber,sddNumber),fitGaussFunc,1500,4500,4);
  fitFuncGaussCuKalpha->SetParameter(0,fitFuncTotal->GetParameter(0));
  fitFuncGaussCuKalpha->SetParameter(1,fitFuncTotal->GetParameter(8));
  fitFuncGaussCuKalpha->SetParameter(2,fitFuncTotal->GetParameter(1));
  fitFuncGaussCuKalpha->SetParameter(3,fitFuncTotal->GetParameter(9));

  fitFuncTailCuKalpha = new TF1(Form("fitFuncTailCuKalpha_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailCuKalpha->SetParameter(0,(fitFuncTotal->GetParameter(0)*fitFuncTotal->GetParameter(4)));
  fitFuncTailCuKalpha->SetParameter(1,fitFuncTotal->GetParameter(8));
  fitFuncTailCuKalpha->SetParameter(2,fitFuncTotal->GetParameter(1));
  fitFuncTailCuKalpha->SetParameter(3,fitFuncTotal->GetParameter(9));
  fitFuncTailCuKalpha->SetParameter(4,fitFuncTotal->GetParameter(3));

  fitFuncGaussCuKbeta = new TF1(Form("fitFuncGaussCuKbeta_%d_%d",busNumber,sddNumber),fitGaussFunc,3300,3500,4);
  fitFuncGaussCuKbeta->SetParameter(0,fitFuncTotal->GetParameter(4));
  fitFuncGaussCuKbeta->SetParameter(1,fitFuncTotal->GetParameter(8));
  fitFuncGaussCuKbeta->SetParameter(2,fitFuncTotal->GetParameter(5));
  fitFuncGaussCuKbeta->SetParameter(3,fitFuncTotal->GetParameter(9));

  fitFuncTailCuKbeta = new TF1(Form("fitFuncTailCuKbeta_%d_%d",busNumber,sddNumber),fitTailFunc,1500,4500,5);
  fitFuncTailCuKbeta->SetParameter(0,(fitFuncTotal->GetParameter(4)*fitFuncTotal->GetParameter(6)));
  fitFuncTailCuKbeta->SetParameter(1,fitFuncTotal->GetParameter(8));
  fitFuncTailCuKbeta->SetParameter(2,fitFuncTotal->GetParameter(5));
  fitFuncTailCuKbeta->SetParameter(3,fitFuncTotal->GetParameter(9));
  fitFuncTailCuKbeta->SetParameter(4,fitFuncTotal->GetParameter(7));

  fitFuncBkgrnd = new TF1(Form("fitFuncBkgrnd_%d_%d",busNumber,sddNumber),fitBkgndFunc,1500,4500,4);
  for(int l=0; l<4; l++) {
    fitFuncBkgrnd->SetParameter(l,fitFuncTotal->GetParameter(l+10));
  }
/*
  //// Linearity ////
  Double_t x[2], y[2];
  x[0]=fitFuncGaussTiKalpha->GetMaximumX();
  y[0]=4509;
  x[1]=fitFuncGaussCuKalpha->GetMaximumX();
  y[1]=8041;

  gLinearity = new TGraph(2,x,y);
  gLinearity->Fit("pol1");

  Double_t slope = gLinearity->GetFunction("pol1")->GetParameter(1);
  Double_t slope_err = gLinearity->GetFunction("pol1")->GetParError(1);
  Double_t offset = gLinearity->GetFunction("pol1")->GetParameter(0);
  Double_t offset_err = gLinearity->GetFunction("pol1")->GetParError(0);

  //// Energy calibration ////

  Double_t binContent = 0.;
  Double_t E_x[10000], E_y[10000];
  for(Int_t l = 1; l < 10001; l++) {
    binContent = hSDD->GetBinContent(l);
    E_x[l-1]=offset+slope*(l-1);
    E_y[l-1]=hSDD->GetBinContent(l);
  }

  gSDDenergy = new TGraph(10000,E_x,E_y);

  TAxis *axis = hEnergySDD->GetXaxis();
  axis->SetLimits(axis->GetXmin()*slope+offset, axis->GetXmax()*slope+offset);

  hEnergySDD->SetAxisRange(0,25000,"X");

  outfile= new TFile(Form("/home/nuclearboy/SIDDHARTA2/SIDDHARTA2_Calibration/calibration/20220310/files/hEnergySDD_bus%d_sdd%d.root",busNumber,sddNumber),"RECREATE");

  outfile->cd();
  hSDD->Write("hADC");
  hEnergySDD->Write("hEnergy");
  outfile->Close();

  //new file with results
  ofstream calib_file;
  calib_file.open(Form("files/calib_bus%d_sdd%d.dat",busNumber,sddNumber), ios::trunc);
  calib_file<<Form("%i %i %g %g %g %g",busNumber,sddNumber,slope,slope_err,offset,offset_err)<<endl;
  calib_file.close();

  //new file with results
  ofstream parameters_file;
  parameters_file.open(Form("files/parameters_bus%d_sdd%d.dat",busNumber,sddNumber), ios::trunc);
  for (Int_t i=0; i<36; i++) {
    parameters_file<<Form("%i %g %g",i, fitFuncTotal->GetParameter(i), fitFuncTotal->GetParError(i))<<endl;
  }
  parameters_file.close();

  hEnergySDD_fit = (TH1F*)hEnergySDD->Clone(Form("hEnergy_bus%d_sdd%d_fit",busNumber,sddNumber));

  fitFuncEnergyTotal = new TF1(Form("fitFuncEnergyTotal_%d_%d",busNumber,sddNumber),"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)+expo(18)+pol1(20)",0,12000);
  fitFuncEnergyTotal->SetParameter(0,54.7);
  fitFuncEnergyTotal->SetParameter(1,4508.);
  fitFuncEnergyTotal->SetParameter(2,66.5);
  fitFuncEnergyTotal->SetParameter(3,11.9);
  fitFuncEnergyTotal->SetParameter(4,4925.);
  fitFuncEnergyTotal->SetParameter(5,77.4);
  fitFuncEnergyTotal->SetParameter(6,4.6);
  fitFuncEnergyTotal->SetParameter(7,5895.);
  fitFuncEnergyTotal->SetParameter(8,84.);
  fitFuncEnergyTotal->SetParameter(9,6.);
  fitFuncEnergyTotal->SetParameter(10,6395.);
  fitFuncEnergyTotal->SetParameter(11,124.);
  fitFuncEnergyTotal->SetParameter(12,364.6);
  fitFuncEnergyTotal->SetParameter(13,8041.);
  fitFuncEnergyTotal->SetParameter(14,78.5);
  fitFuncEnergyTotal->SetParameter(15,65.);
  fitFuncEnergyTotal->SetParameter(16,8905.);
  fitFuncEnergyTotal->SetParameter(17,86.);
  fitFuncEnergyTotal->SetParameter(18,-35.2);
  fitFuncEnergyTotal->SetParameter(19,0.0030);
  fitFuncEnergyTotal->SetParameter(20,1.46);
  fitFuncEnergyTotal->SetParameter(21,-6.62e-05);

  hEnergySDD_fit->Fit(Form("fitFuncEnergyTotal_%d_%d",busNumber,sddNumber),"","",2000,12000);

  Double_t chi2_E = (fitFuncEnergyTotal->GetChisquare())/(fitFuncEnergyTotal->GetNDF());
*/
/////////////////////////////////////////HISTOGRAMS/////////////////////////////////////////

  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1,0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTextFont(42);

  //
  myCanvas[0] = new TCanvas;
  //myCanvas[0]->SetLogy();

  hSDD->SetTitle(Form("BUS: %d, SDD: %d",busNumber,sddNumber));
  hSDD->GetXaxis()->SetTitle("ADC [channel]");
  hSDD->GetXaxis()->SetTitleSize(0.05);
  hSDD->GetXaxis()->SetTitleOffset(1.);
  hSDD->GetXaxis()->SetLabelSize(0.04);
  //hSDD->SetAxisRange(1500,4500,"X");
  //hSDD->GetXaxis()->SetRangeUser(1500.,4500.);
  hSDD->GetYaxis()->SetTitle("counts");
  hSDD->GetYaxis()->SetTitleSize(0.05);
  hSDD->GetYaxis()->SetTitleOffset(1.);
  hSDD->GetYaxis()->SetLabelSize(0.04);

  hSDD->SetLineWidth(1);
  hSDD->SetLineColor(1);
  hSDD->Draw();

  hCrosstalk->SetLineWidth(1);
  hCrosstalk->SetLineColor(4);
  hCrosstalk->Draw("same");

  myLegend[0] = new TLegend(0.625, 0.790, 0.900, 0.900);
  myLegend[0]->SetFillStyle(1001); myLegend[0]->SetLineColor(1); myLegend[0]->SetFillColor(0); myLegend[0]->SetTextSize(0.037);
  myLegend[0]->AddEntry(hSDD, "Data", "l");
  myLegend[0]->AddEntry(hCrosstalk, "Crosstalk", "l");
  myLegend[0]->Draw("same");

  myCanvas[0]->Print(Form("plots/hADC_bus%d_sdd%d.png",busNumber,sddNumber), "png");
/*
  //
  myCanvas[1] = new TCanvas;
  myCanvas[1]->SetLogy();

  hSDD_copy->SetTitle(Form("BUS: %d, SDD: %d",busNumber,sddNumber));
  hSDD_copy->GetXaxis()->SetTitle("ADC [channel]");
  hSDD_copy->GetXaxis()->SetTitleSize(0.05);
  hSDD_copy->GetXaxis()->SetTitleOffset(1.);
  hSDD_copy->GetXaxis()->SetLabelSize(0.04);
  hSDD_copy->SetAxisRange(1700,2400,"X");
  //hSDD_copy->GetXaxis()->SetRangeUser(1700.,2400.);
  hSDD_copy->GetYaxis()->SetTitle("counts / 4 channel");
  hSDD_copy->GetYaxis()->SetTitleSize(0.05);
  hSDD_copy->GetYaxis()->SetTitleOffset(1.);
  hSDD_copy->GetYaxis()->SetLabelSize(0.04);

  hSDD_copy->SetLineWidth(1);
  hSDD_copy->SetLineColor(1);
  hSDD_copy->Draw();
*/
  fitFuncTotal->SetLineColor(2);
  fitFuncTotal->SetLineWidth(2);
  //fitFuncTotal->Draw("same");

  fitFuncBkgrnd->SetLineColor(7);
  fitFuncBkgrnd->SetLineWidth(1);
  //fitFuncBkgrnd->Draw("same");

  fitFuncBkgrndStart->SetLineColor(6);
  fitFuncBkgrndStart->SetLineWidth(1);
  //fitFuncBkgrndStart->Draw("same");
/*
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
/*
  TPaveText *kTi_a = new TPaveText(binPeak[0]*1.02,peakValue[0], binPeak[0]*1.02,peakValue[0],"TiK_a");
  kTi_a->SetTextSize(0.03);
  kTi_a->SetTextColor(1);
  kTi_a->SetTextAlign(22);
  kTi_a->AddText("TiK_{#alpha}");
  kTi_a->Draw("same");

  TPaveText *kTi_b = new TPaveText(binPeak[1]*1.02,peakValue[1], binPeak[1]*1.02,peakValue[1],"TiK_b");
  kTi_b->SetTextSize(0.03);
  kTi_b->SetTextColor(1);
  kTi_b->SetTextAlign(22);
  kTi_b->AddText("TiK_{#beta}");
  kTi_b->Draw("same");

  myLegend[1] = new TLegend(0.600, 0.695, 0.900, 0.900);
  myLegend[1]->SetFillStyle(1001); myLegend[1]->SetLineColor(1); myLegend[1]->SetFillColor(0); myLegend[1]->SetTextSize(0.037);
  myLegend[1]->AddEntry(hSDD, "Data", "l");
  myLegend[1]->AddEntry(fitFuncTotal, "Global fit function", "l");
  myLegend[1]->AddEntry(fitFuncGaussTiKalpha, "Gaussian function", "l");
  myLegend[1]->AddEntry(fitFuncTailTiKalpha, "Tail function", "l");
  myLegend[1]->AddEntry(fitFuncBkgrnd, "Background", "l");
  myLegend[1]->Draw("same");

  myCanvas[1]->Print(Form("plots/hADC_Ylog_bus%d_sdd%d_Ti.png",busNumber,sddNumber), "png");
*/
  //
  myCanvas[2] = new TCanvas;
  myCanvas[2]->SetLogy();

  hSDD_copy->SetAxisRange(2700,3700,"X");
  //hSDD_copy->GetXaxis()->SetRangeUser(2700.,3700.);
  hSDD_copy->Draw("same");

  fitFuncTotal->Draw("same");

  fitFuncBkgrnd->Draw("same");
  fitFuncBkgrndStart->Draw("same");

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

  TPaveText *kCu_a = new TPaveText(binPeak[4]*1.02,peakValue[4], binPeak[4]*1.02,peakValue[4],"CuK_a");
  kCu_a->SetTextSize(0.03);
  kCu_a->SetTextColor(1);
  kCu_a->SetTextAlign(22);
  kCu_a->AddText("CuK_{#alpha}");
  kCu_a->Draw("same");

  TPaveText *kCu_b = new TPaveText(binPeak[5]*1.02,peakValue[5]*0.75, binPeak[5]*1.02,peakValue[5]*0.75,"TiK_b");
  kCu_b->SetTextSize(0.03);
  kCu_b->SetTextColor(1);
  kCu_b->SetTextAlign(22);
  kCu_b->AddText("CuK_{#beta}");
  kCu_b->Draw("same");

  myLegend[2] = new TLegend(0.600, 0.695, 0.900, 0.900);
  myLegend[2]->SetFillStyle(1001); myLegend[2]->SetLineColor(1); myLegend[2]->SetFillColor(0); myLegend[2]->SetTextSize(0.037);
  myLegend[2]->AddEntry(hSDD, "Data", "l");
  myLegend[2]->AddEntry(fitFuncTotal, "Global fit function", "l");
  myLegend[2]->AddEntry(fitFuncGaussCuKalpha, "Gaussian function", "l");
  myLegend[2]->AddEntry(fitFuncTailCuKalpha, "Tail function", "l");
  myLegend[2]->AddEntry(fitFuncBkgrnd, "Background", "l");
  myLegend[2]->Draw("same");

  myCanvas[2]->Print(Form("plots/hADC_Ylog_bus%d_sdd%d_Cu.png",busNumber,sddNumber), "png");
/*
  //
  myCanvas[3] = new TCanvas;
  //myCanvas[3]->SetLogy();

  hSDD_copy->SetTitle(Form("BUS: %d, SDD: %d, #chi^{2}/NDf: %f",busNumber,sddNumber,chi2));
  hSDD_copy->SetAxisRange(1700,3700,"X");
  hSDD_copy->Draw("same");

  fitFuncTotal->Draw("same");

  fitFuncBkgrnd->Draw("same");
  fitFuncBkgrndStart->Draw("same");

  fitFuncGaussTiKalpha->Draw("same");
  fitFuncGaussTiKbeta->Draw("same");

  fitFuncGaussMnKalpha->SetLineColor(51);
  fitFuncGaussMnKalpha->SetLineWidth(1);
  fitFuncGaussMnKalpha->Draw("same");

  fitFuncGaussFeKalpha->SetLineColor(51);
  fitFuncGaussFeKalpha->SetLineWidth(1);
  fitFuncGaussFeKalpha->Draw("same");

  fitFuncGaussCuKalpha->Draw("same");
  fitFuncGaussCuKbeta->Draw("same");

  fitFuncTailTiKalpha->Draw("same");
  fitFuncTailTiKbeta->Draw("same");
  fitFuncTailCuKalpha->Draw("same");
  fitFuncTailCuKbeta->Draw("same");

  kTi_a = new TPaveText(binPeak[0],peakValue[0]*1.3, binPeak[0],peakValue[0]*1.3,"TiK_a");
  kTi_a->SetTextSize(0.03);
  kTi_a->SetTextColor(1);
  kTi_a->SetTextAlign(22);
  kTi_a->AddText("TiK_{#alpha}");
  kTi_a->Draw("same");

  kTi_b = new TPaveText(binPeak[1],peakValue[1]*1.3, binPeak[1],peakValue[1]*1.3,"TiK_b");
  kTi_b->SetTextSize(0.03);
  kTi_b->SetTextColor(1);
  kTi_b->SetTextAlign(22);
  kTi_b->AddText("TiK_{#beta}");
  kTi_b->Draw("same");

  TPaveText *kMn_a = new TPaveText(binPeak[2],peakValue[2]*1.3, binPeak[2],peakValue[2]*1.3,"MnK_a");
  kMn_a->SetTextSize(0.03);
  kMn_a->SetTextColor(1);
  kMn_a->SetTextAlign(22);
  kMn_a->AddText("MnK_{#alpha}");
  kMn_a->Draw("same");

  TPaveText *kFe_a = new TPaveText(binPeak[3],peakValue[3]*1.3, binPeak[3],peakValue[3]*1.3,"FeK_a");
  kFe_a->SetTextSize(0.03);
  kFe_a->SetTextColor(1);
  kFe_a->SetTextAlign(22);
  kFe_a->AddText("FeK_{#alpha}");
  kFe_a->Draw("same");

  kCu_a = new TPaveText(binPeak[4]*1.03,peakValue[4], binPeak[4]*1.03,peakValue[4],"CuK_a");
  kCu_a->SetTextSize(0.03);
  kCu_a->SetTextColor(1);
  kCu_a->SetTextAlign(22);
  kCu_a->AddText("CuK_{#alpha}");
  kCu_a->Draw("same");

  kCu_b = new TPaveText(binPeak[5],peakValue[5]*1.3, binPeak[5],peakValue[5]*1.3,"CuK_b");
  kCu_b->SetTextSize(0.03);
  kCu_b->SetTextColor(1);
  kCu_b->SetTextAlign(22);
  kCu_b->AddText("CuK_{#beta}");
  kCu_b->Draw("same");

  myLegend[3] = new TLegend(0.285, 0.695, 0.565, 0.900);
  myLegend[3]->SetFillStyle(1001); myLegend[3]->SetLineColor(1); myLegend[3]->SetFillColor(0); myLegend[3]->SetTextSize(0.037);
  myLegend[3]->AddEntry(hSDD, "Data", "l");
  myLegend[3]->AddEntry(fitFuncTotal, "Global fit function", "l");
  myLegend[3]->AddEntry(fitFuncGaussTiKalpha, "Gaussian function", "l");
  myLegend[3]->AddEntry(fitFuncTailTiKalpha, "Tail function", "l");
  myLegend[3]->AddEntry(fitFuncBkgrnd, "Background", "l");
  myLegend[3]->Draw("same");

  myCanvas[3]->Print(Form("plots/hADC_Ylog_bus%d_sdd%d.png",busNumber,sddNumber), "png");

  //

  myCanvas[4] = new TCanvas;

  gLinearity->SetTitle(Form("BUS: %d, SDD: %d",busNumber,sddNumber));
  gLinearity->GetXaxis()->SetTitle("ADC [channel]");
  gLinearity->GetXaxis()->SetTitleSize(0.05);
  gLinearity->GetXaxis()->SetTitleOffset(1.);
  gLinearity->GetXaxis()->SetLabelSize(0.04);
  //gLinearity->GetXaxis()->SetRangeUser(1000.,5000.);
  gLinearity->GetYaxis()->SetTitle("energy [eV]");
  gLinearity->GetYaxis()->SetTitleSize(0.05);
  gLinearity->GetYaxis()->SetTitleOffset(1.1);
  gLinearity->GetYaxis()->SetLabelSize(0.04);
  gLinearity->GetYaxis()->SetRangeUser(4150,8500);

  //gLinearity->SetLineWidth(2);
  //gLinearity->SetLineColor(5);
  gLinearity->SetMarkerStyle(23);
  gLinearity->SetMarkerSize(1);
  gLinearity->SetMarkerColor(1);
  gLinearity->Draw();

  TPaveText *Lin_Ti = new TPaveText(binPeak[0],4750, binPeak[0],4750,"Ti_a");
  Lin_Ti->SetTextSize(0.04);
  Lin_Ti->SetTextColor(4);
  Lin_Ti->SetTextAlign(22);
  Lin_Ti->AddText("TiK_{#alpha}");
  Lin_Ti->Draw("same");

  TPaveText *Lin_Cu = new TPaveText(binPeak[4],8300, binPeak[4],8300,"Cu_a");
  Lin_Cu->SetTextSize(0.04);
  Lin_Cu->SetTextColor(4);
  Lin_Cu->SetTextAlign(22);
  Lin_Cu->AddText("CuK_{#alpha}");
  Lin_Cu->Draw("same");

  TPaveText *descr = new TPaveText(2000.,7300,2500,8000);
  TText *tl02 = descr->AddText(Form("offset = %.4f #pm %.4f",offset,offset_err));
  TText *tl01 = descr->AddText(Form("slope = %.4f #pm %.4f",slope,slope_err));
  descr->SetFillColor(0);
  tl01->SetTextFont(42);
  tl02->SetTextFont(42);
  tl01->SetTextSize(0.04);
  tl02->SetTextSize(0.04);
  tl01->SetTextAlign(12);
  tl02->SetTextAlign(12);
  descr->Draw("same");

  myCanvas[4]->Print(Form("plots/gLinearity_bus%d_sdd%d.png",busNumber,sddNumber), "png");

  //

  myCanvas[5] = new TCanvas;

  gSDDenergy->SetTitle(Form("BUS: %d, SDD: %d",busNumber,sddNumber));
  gSDDenergy->GetXaxis()->SetTitle("energy [eV]");
  gSDDenergy->GetXaxis()->SetTitleSize(0.05);
  gSDDenergy->GetXaxis()->SetTitleOffset(1.);
  gSDDenergy->GetXaxis()->SetLabelSize(0.04);
  gSDDenergy->GetXaxis()->SetRangeUser(2000.,12000.);
  gSDDenergy->GetYaxis()->SetTitle("counts");
  gSDDenergy->GetYaxis()->SetTitleSize(0.05);
  gSDDenergy->GetYaxis()->SetTitleOffset(1.);
  gSDDenergy->GetYaxis()->SetLabelSize(0.04);
  gSDDenergy->GetYaxis()->SetRangeUser(0.,1.05*hSDD->GetMaximum());

  gSDDenergy->SetLineWidth(1);
  gSDDenergy->SetLineColor(1);
  gSDDenergy->Draw();

  kTi_a = new TPaveText(4508,(peakValue[0]/4)*1.5, 4508,(peakValue[0]/4)*1.5,"TiK_a");
  kTi_a->SetTextSize(0.03);
  kTi_a->SetTextColor(1);
  kTi_a->SetTextAlign(22);
  kTi_a->AddText("TiK_{#alpha}");
  kTi_a->Draw("same");

  kTi_b = new TPaveText(5000,peakValue[1]/2, 5000,peakValue[1]/2,"TiK_b");
  kTi_b->SetTextSize(0.03);
  kTi_b->SetTextColor(1);
  kTi_b->SetTextAlign(22);
  kTi_b->AddText("TiK_{#beta}");
  kTi_b->Draw("same");

  kMn_a = new TPaveText(5900,peakValue[2]*1.3, 5900,peakValue[2]*1.3,"MnK_a");
  kMn_a->SetTextSize(0.03);
  kMn_a->SetTextColor(1);
  kMn_a->SetTextAlign(22);
  kMn_a->AddText("MnK_{#alpha}");
  kMn_a->Draw("same");

  kFe_a = new TPaveText(6395,peakValue[3]*1.3, 6395,peakValue[3]*1.3,"FeK_a");
  kFe_a->SetTextSize(0.03);
  kFe_a->SetTextColor(1);
  kFe_a->SetTextAlign(22);
  kFe_a->AddText("FeK_{#alpha}");
  kFe_a->Draw("same");

  kCu_a = new TPaveText(8490,peakValue[4]/3.9, 8490,peakValue[4]/3.9,"CuK_a");
  kCu_a->SetTextSize(0.03);
  kCu_a->SetTextColor(1);
  kCu_a->SetTextAlign(22);
  kCu_a->AddText("CuK_{#alpha}");
  kCu_a->Draw("same");

  kCu_b = new TPaveText(8910,(peakValue[5]/4)*1.3, 8910,(peakValue[5]/4)*1.3,"CuK_b");
  kCu_b->SetTextSize(0.03);
  kCu_b->SetTextColor(1);
  kCu_b->SetTextAlign(22);
  kCu_b->AddText("CuK_{#beta}");
  kCu_b->Draw("same");

  myLegend[5] = new TLegend(0.775, 0.845, 0.900, 0.900);
  myLegend[5]->SetFillStyle(1001); myLegend[5]->SetLineColor(1); myLegend[5]->SetFillColor(0); myLegend[5]->SetTextSize(0.04);
  myLegend[5]->AddEntry(gSDDenergy, "data", "l");
  myLegend[5]->Draw("same");

  myCanvas[5]->Print(Form("plots/hEnergy_bus%d_sdd%d.png",busNumber,sddNumber), "png");

  //

  myCanvas[6] = new TCanvas;

  hEnergySDD->SetTitle(Form("BUS: %d, SDD: %d, #chi^{2}/NDf: %f",busNumber,sddNumber,chi2_E));
  hEnergySDD->GetXaxis()->SetTitle("energy [eV]");
  hEnergySDD->GetXaxis()->SetTitleSize(0.05);
  hEnergySDD->GetXaxis()->SetTitleOffset(1.);
  hEnergySDD->GetXaxis()->SetLabelSize(0.04);
  hEnergySDD->SetAxisRange(2000,12000,"X");
  //hEnergySDD->GetXaxis()->SetRangeUser(2800.,3700.);
  hEnergySDD->GetYaxis()->SetTitle("counts");
  hEnergySDD->GetYaxis()->SetTitleSize(0.05);
  hEnergySDD->GetYaxis()->SetTitleOffset(1.);
  hEnergySDD->GetYaxis()->SetLabelSize(0.04);
  //hEnergySDD->GetYaxis()->SetRangeUser(0.,1.05*hEnergySDD->GetMaximum());

  hEnergySDD->SetLineWidth(1);
  hEnergySDD->SetLineColor(1);
  hEnergySDD->Draw();

  fitFuncEnergyTotal->SetLineWidth(1);
  fitFuncEnergyTotal->SetLineColor(2);
  fitFuncEnergyTotal->Draw("same");

  kTi_a->Draw("same");
  kTi_b->Draw("same");
  kMn_a->Draw("same");
  kFe_a->Draw("same");
  kCu_a->Draw("same");
  kCu_b->Draw("same");

  myLegend[6] = new TLegend(0.745, 0.805, 0.900, 0.900);
  myLegend[6]->SetFillStyle(1001); myLegend[6]->SetLineColor(1); myLegend[6]->SetFillColor(0); myLegend[6]->SetTextSize(0.04);
  myLegend[6]->AddEntry(gSDDenergy, "data", "l");
  myLegend[6]->AddEntry(fitFuncEnergyTotal, "fit", "l");
  myLegend[6]->Draw("same");

  myCanvas[6]->Print(Form("plots/hEnergySDD_bus%d_sdd%d_test.png",busNumber,sddNumber), "png");
*/
}
