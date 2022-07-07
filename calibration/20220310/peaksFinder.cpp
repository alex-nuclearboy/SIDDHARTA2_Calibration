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

void peaksFinder() {

  int nBUS = 7;
  int nSDD = 65;

  ofstream new_file;

  // Read ROOT file
  f[0] = new TFile("/home/nuclearboy/SIDDHARTA2/SIDDHARTA2_Calibration/rootfiles/SIDDHARTA2_xray/output/20220310/hist_20220310_1900_0310_2000_xray_25kv_70ua_tube1_p1.root","READ");

  //ADC histo booking
  int nbinsadc = 10000;
  int minadc = 0.;
  int maxadc = 10000;
  //rebin:
  int rebinFactor = 8;
  nbinsadc = nbinsadc/rebinFactor;

  TH1F* hADC[nBUS][nSDD];
  TH1F* hADCfit[nBUS][nSDD];
  TH1F* hEnergy[nBUS][nSDD];

  for(int iBUS=1;iBUS<nBUS;iBUS++) {
    for(int iSDD=1;iSDD<nSDD;iSDD++) {
      // Take SDDs histograms
      hADC[iBUS][iSDD] = (TH1F*)f[0]->Get(Form("bus%d_sdd%d",iBUS,iSDD));
      if(!hADC[iBUS][iSDD]) continue;
      hADC[iBUS][iSDD]->Rebin(rebinFactor);
      hADCfit[iBUS][iSDD] = (TH1F*)hADC[iBUS][iSDD]->Clone(Form("hADCfit_bus%d_sdd%d",iBUS,iSDD));
      hEnergy[iBUS][iSDD] = (TH1F*)hADC[iBUS][iSDD]->Clone(Form("hEnergy_bus%d_sdd%d",iBUS,iSDD));
      hADC[iBUS][iSDD]->SetName(Form("ADC_bus%i_sdd%i",iBUS,iSDD));
      hADC[iBUS][iSDD]->GetXaxis()->SetTitle("ADC [channel]");
      hADC[iBUS][iSDD]->SetLineColor(1);
      hADC[iBUS][iSDD]->SetLineWidth(1);
    }
  }

  double eTi_KA = 4508.83;
  double eTi_KB = 4931.81;
  double eCu_KA = 8041.05;
  double eCu_KB = 8905.29;
  double eMn_KA = 5895.23;
  double eFe_KA = 6399.47;
  double eFe_KB = 7058.0;

  const int nPFPeaksMAX = 13; //number of MAX peaks for peak finder
  const int nPFPeaks = 4;     //number DESIRED of peaks for peak finder
  // ==> SELECT HERE THE RANGE OF THE PeakFinder
  float xminPeakFinder = 1600;  // Region of peak finding
  float xmaxPeakFinder = 3900;

  int sigmaPeakFinder = 20;   //sigma for peak Finder, in adc channels
  sigmaPeakFinder = sigmaPeakFinder/rebinFactor;

  float InitThresholdPeakFinder = 0.01; //initial thresholdPar for peak finder, std in TSpectrum is 0.05
  float InitTolerance = 0.05; // 5% -> Tolerance to check that peak assumption is correct

  Double_t *xpeaks;
  Double_t *ypeaks;

  float xadc[nPFPeaksMAX] = {};
  float yadc[nPFPeaksMAX] = {};
  float PFPeakE[nPFPeaks] = {};

  TF1 *fBkgnd[nBUS][nSDD];
  TF1 *fitFuncTotalN[nBUS][nSDD];

  for(int iBUS=1;iBUS<nBUS;iBUS++) {

    for(int iSDD=1;iSDD<nSDD;iSDD++) {

      TH1F* thehisto;

      if(!hADC[iBUS][iSDD]) continue;
      thehisto = (TH1F*) hADC[iBUS][iSDD]->Clone("thehisto");  thehisto->SetName("thehisto");
      int minstatsForCalib = 10000;
      if(thehisto->GetEntries()<minstatsForCalib) continue;

      thehisto->GetXaxis()->SetRangeUser(xminPeakFinder,xmaxPeakFinder);

      TSpectrum *spectrum = new TSpectrum(nPFPeaksMAX);

      PFPeakE[0] = eTi_KA;
      PFPeakE[1] = eCu_KA;
      PFPeakE[2] = eCu_KB;

      Int_t nfound = 0;
      float thresholdPeakFinder = InitThresholdPeakFinder;
      int nPFtries = 0;
      while(nfound<nPFPeaks&&nPFtries<50) {
        nfound = spectrum->Search(thehisto,sigmaPeakFinder,"",thresholdPeakFinder);
        printf("B: %d S: %d. Found %d candidate peaks to fit\n",iBUS,iSDD,nfound);
        thresholdPeakFinder = thresholdPeakFinder*.1; //TSpectrum std=0.05. Change it til it finds the peaks
        nPFtries++;
      }

      if(nPFtries>=50) {
        cout<<" PEAK FINDER DIDNT WORK, ==> CONTINUE"<<endl;
        continue;
      }
      xpeaks = spectrum->GetPositionX();  //array with X-positions of the centroids found by TSpectrum
      ypeaks = spectrum->GetPositionY();  //array with X-positions of the centroids found by TSpectrum

      //reorder in adc counts and check if compatible with maximum peaks assumption
      Float_t themin = 999999.;
      Int_t imin = 0;
      for(int i =0;i<nfound;i++) { //find the smallest, write it in xadc and remove it:
        themin = 999999.;
        for(int j =0;j<nfound;j++) {
          if(xpeaks[j]<themin) {
            themin = xpeaks[j];
            imin = j;
          }
        }
        xadc[i] = xpeaks[imin];
        yadc[i] = ypeaks[imin];
        xpeaks[imin]=999999.;
      }
      //print them, now ordered:
      for(int i =0;i<nfound;i++){cout<<" found peak ADC="<<xadc[i]<<" height "<<yadc[i]<<endl;}
      cout<<endl;

      //check if the peaks found are compatible with the assumption:
      float eDist10 = PFPeakE[1]-PFPeakE[0];
      float eDist21 = PFPeakE[2]-PFPeakE[1];
      float Erelation = eDist21/eDist10;
      //make trios out of all found peaks:
      float GPF =0.;
      float G0PF =0;
      bool TestPassed = false;
      int ipeak0 = -1;
      int ipeak1 = -1;
      int ipeak2 = -1;
      for(int i0=0; i0<nfound; i0++) {
        for(int i1=i0+1; i1<nfound; i1++) {
          for(int i2=i1+1; i2<nfound; i2++) {
            cout<<endl<<" -> trying trio "<<i0<<" "<<i1<<" "<<i2<<endl;
            float xDist10 = xadc[i1]-xadc[i0];
            float xDist21 = xadc[i2]-xadc[i1];
            float ADCrelation = xDist21/xDist10;
            cout<<" Check assumption: Energy relation "<<Erelation<<" vs ADC relation "<<ADCrelation<<endl;

            //Define tolerance parameter
            float tol = InitTolerance; // 5%
            bool TolerancePass = true;
            if(fabs(1.-(Erelation/ADCrelation))>tol) TolerancePass = false;
            if(!TolerancePass) cout<<" Tolerance not passed: tolERROR!!!  "<<endl;

            // Get Peak Finder calibration Offset G0PF and Slope GPF
            float Dadc = xadc[i0]-xadc[i1];
            float De = PFPeakE[0]-PFPeakE[1];
            GPF = De/Dadc;
            G0PF = -1.*xadc[i0]*GPF+PFPeakE[0];
            cout<<"PeakFinder cal, offset G0PF "<<G0PF<<" slope GPF "<<GPF<<endl;

            //define an acceptable gain and offset and check if tolerance, G, and G0 are ok:
            float mingoodG = 2.9;float maxgoodG = 3.9;
            float mingoodG0 = -3000;float maxgoodG0 = -1000;
            if(GPF<maxgoodG&&GPF>mingoodG&&G0PF<maxgoodG0&&G0PF>mingoodG0&&TolerancePass) {
              cout<<" -- TEST PASSED! --"<<endl<<endl;
              TestPassed = true;

              ipeak0 = i0;
              ipeak1 = i1;
              ipeak2 = i2;

            }
            if(TestPassed) break;
          }
          if(TestPassed) break;
        }
        if(TestPassed) break;
      }
      cout<<"peak#"<<ipeak0<<": "<<xadc[ipeak0]<<" => "<< xadc[ipeak0]*GPF+G0PF<<endl;
      cout<<"peak#"<<ipeak1<<": "<<xadc[ipeak1]<<" => "<< xadc[ipeak1]*GPF+G0PF<<endl;
      cout<<"peak#"<<ipeak2<<": "<<xadc[ipeak2]<<" => "<< xadc[ipeak2]*GPF+G0PF<<endl;

      float diff;
      float dMin=0.97;
      float dMax=1.03;
      float peakTiKa, peakTiKb, peakCuKa, peakCuKb, peakMnKa, peakFeKa;
      for(int k=0;k<nfound;k++){
        diff = xadc[k]*GPF+G0PF;
        if (diff/eTi_KA > dMin && diff/eTi_KA < dMax ) {peakTiKa = xadc[k];} else
        if (diff/eTi_KB > dMin && diff/eTi_KB < dMax ) {peakTiKb = xadc[k];} else
        if (diff/eCu_KA > dMin && diff/eCu_KA < dMax ) {peakCuKa = xadc[k];} else
        if (diff/eCu_KB > dMin && diff/eCu_KB < dMax ) {peakCuKb = xadc[k];}
      }

      new_file.open("check.dat", ios::app);
      new_file<<Form("B#%d S#%d peakTiKa (%g): %g",iBUS,iSDD,peakTiKa*GPF+G0PF,peakTiKa)<<endl;
      new_file<<Form("B#%d S#%d peakTiKb (%g): %g",iBUS,iSDD,peakTiKb*GPF+G0PF,peakTiKb)<<endl;
      new_file<<Form("B#%d S#%d peakCuKa (%g): %g",iBUS,iSDD,peakCuKa*GPF+G0PF,peakCuKa)<<endl;
      new_file<<Form("B#%d S#%d peakCuKb (%g): %g\n",iBUS,iSDD,peakCuKb*GPF+G0PF,peakCuKb)<<endl;
      new_file.close();

      thehisto->Delete();

      // Estimate background

      fBkgnd[iBUS][iSDD] = new TF1(Form("fBkgnd_bus%d_sdd%d",iBUS,iSDD),"pol0(0)+expo(1)",1500,4500);
      hADCfit[iBUS][iSDD]->Fit(Form("fBkgnd_bus%d_sdd%d",iBUS,iSDD),"","",1500,4500);

      fitFuncTotalN[iBUS][iSDD] = new TF1(Form("fitFuncTotal_%d_%d",iBUS,iSDD),"gaus(0)+gaus(3)+gaus(6)+gaus(9)+pol0(12)+expo(13)",1500,4500);

      fitFuncTotalN[iBUS][iSDD]->SetParameter(1,peakTiKa);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(2,20);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(4,peakTiKb);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(5,20);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(7,peakCuKa);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(8,20);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(10,peakCuKb);
      fitFuncTotalN[iBUS][iSDD]->SetParameter(11,20);

      for (Int_t k = 0; k < 3; k++) {
        fitFuncTotalN[iBUS][iSDD]->SetParameter(k+12,fBkgnd[iBUS][iSDD]->GetParameter(k));
      }

      fitFuncTotalN[iBUS][iSDD]->SetNpx(1000);

      hADCfit[iBUS][iSDD]->Fit(fitFuncTotalN[iBUS][iSDD],"","",1500,4500);

    }
  }

  myCanvas[3] = new TCanvas;
  myCanvas[3]->SetLogy();
  myCanvas[3]->SetGrid();

  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      if(!hADC[iBUS][iSDD]) continue;
      if(hADCfit[iBUS][iSDD]->GetEntries()<10000) continue;
      if(!fBkgnd[iBUS][iSDD]) continue;
      hADC[iBUS][iSDD]->SetTitle(Form("BUS: %d, SDD: %d",iBUS,iSDD));
      hADC[iBUS][iSDD]->SetAxisRange(1500,4500,"X");
      hADC[iBUS][iSDD]->GetYaxis()->SetTitle("counts / 8 channel");

      hADC[iBUS][iSDD]->Draw();

      fBkgnd[iBUS][iSDD]->SetLineColor(7);
      fBkgnd[iBUS][iSDD]->SetLineWidth(1);
      fBkgnd[iBUS][iSDD]->Draw("same");

      fitFuncTotalN[iBUS][iSDD]->SetLineColor(2);
      fitFuncTotalN[iBUS][iSDD]->SetLineWidth(2);
      fitFuncTotalN[iBUS][iSDD]->Draw("same");

      //myCanvas[3]->Print(Form("plots2/hADC_Ylog_bus%d_sdd%d.png",iBUS,iSDD), "png");
    }
  }



}
