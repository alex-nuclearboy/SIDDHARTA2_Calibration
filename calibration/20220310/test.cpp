/***********************************************
* SIDDHARTA-2 Experiment
* Aleksander K.                 2022-07
* Licensed under the Apache License, Version 2.0
***********************************************/

// Peak Finder

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

//Energies of X-ray emmission lines
Double_t eTiKa = 4508.83;
Double_t eTiKb = 4931.81;
Double_t eCuKa = 8041.05;
Double_t eCuKb = 8905.29;
Double_t eMnKa = 5895.23;
Double_t eFeKa = 6399.47;
Double_t eFeKb = 7058.0;

Int_t npeaks;
Double_t findPeaks(Double_t *x, Double_t *par) {
  Double_t result = par[0] + TMath::Exp(par[1]+par[2]*x[0]); //background
  for (Int_t p=0; p<npeaks; p++) { //gauss
    Double_t norm  = par[3*p+3];
    Double_t mean  = par[3*p+4];
    Double_t sigma = par[3*p+5];
    result += norm*TMath::Gaus(x[0],mean,sigma); //background + gauss
  }
  return result;
}

void peakFinder(){

  const Int_t nBUS = 7;   //last BUS #6
  const Int_t nSDD = 65;  //last SDD #64

  ofstream new_file,new_file2;

  // Read ROOT file
  f[0] = new TFile("/home/nuclearboy/SIDDHARTA2/SIDDHARTA2_Calibration/rootfiles/SIDDHARTA2_xray/output/20220310/hist_20220310_1900_0310_2000_xray_25kv_70ua_tube1_p1.root","READ");

  //ADC histograms
  Int_t nbinsadc  = 10000;
  Int_t minadc    = 0.;
  Int_t maxadc    = 10000.;

  Int_t rebinFactor = 4; //rebining
  nbinsadc = nbinsadc/rebinFactor;

  TH1F* hADC[nBUS][nSDD];
  TH1F* hADCfit[nBUS][nSDD];

  TF1* fitFuncTotalN[nBUS][nSDD];
  TF1* fBkgnd[nBUS][nSDD];
  TF1 *fPreFit[nBUS][nSDD];

  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      hADC[iBUS][iSDD] = (TH1F*)f[0]->Get(Form("bus%d_sdd%d",iBUS,iSDD));
      if(!hADC[iBUS][iSDD]) continue;
      hADC[iBUS][iSDD]->Rebin(rebinFactor);
      hADCfit[iBUS][iSDD] = (TH1F*)hADC[iBUS][iSDD]->Clone(Form("hADCfit_bus%d_sdd%d",iBUS,iSDD));
      hADC[iBUS][iSDD]->GetXaxis()->SetTitle("ADC [channel]");
      hADC[iBUS][iSDD]->SetLineColor(1);
      hADC[iBUS][iSDD]->SetLineWidth(2);
    }
  }

  const Int_t nPFPeaksMAX = 13; //MAX number of peaks for the Peak Finder
  const Int_t nPFPeaks    = 4;  //DESIRED number of peaks

  // RANGE of ADC for the Peak Finder
  Float_t xminPeakFinder  = 1700;
  Float_t xmaxPeakFinder  = 3700;

  Int_t sigmaPeakFinder   = 20; //sigma for tge Peak Finder (in ADC channels)
  sigmaPeakFinder = sigmaPeakFinder/rebinFactor;

  Float_t initThresholdPeakFinder = 0.01; //initial threshold parameter for the Peak Finder (std in TSpectrum is 0.05)
  Float_t initTolerance = 0.05;           //tolerance to check that the peak assumption is correct (5%)

  Double_t *xpeaks;
  Double_t *ypeaks;

  Float_t xadc[nPFPeaksMAX] ={};
  Float_t yadc[nPFPeaksMAX] ={};
  Float_t PFPeakE[nPFPeaks] ={};

  TGraph* gPeakFinderG[nBUS];
  TGraph* gPeakFinderG0[nBUS];

  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    gPeakFinderG[iBUS]= new TGraph();
    gPeakFinderG[iBUS]->SetName(Form("gPeakFinderG_%i",iBUS));
    gPeakFinderG0[iBUS]= new TGraph();
    gPeakFinderG0[iBUS]->SetName(Form("gPeakFinderG0_%i",iBUS));
  }

  // Assumption for the Peak Finder peaks:  TiKa - CuKa - CuKb
  PFPeakE[0] = eTiKa;
  PFPeakE[1] = eCuKa;
  PFPeakE[2] = eCuKb;

  Int_t ipoint = 0;

  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      if(!hADC[iBUS][iSDD]) continue;
      TH1F* histoPF;
      histoPF = (TH1F*) hADC[iBUS][iSDD]->Clone("histoPF"); histoPF->SetName("histoPF");
      Int_t minStatsForCalib = 10000;   //min statistics for calibration
      if(histoPF->GetEntries()<minStatsForCalib) continue;
      histoPF->SetAxisRange(xminPeakFinder,xmaxPeakFinder,"X");

      // Use TSpectrum to find the peak candidates
      Int_t nfound = 0;
      Float_t thresholdPeakFinder = initThresholdPeakFinder;
      Int_t nPFtries = 0;
      TSpectrum *spectrum = new TSpectrum(nPFPeaksMAX);
      while(nfound<nPFPeaks && nPFtries<15) {
        nfound = spectrum->Search(histoPF,sigmaPeakFinder,"",thresholdPeakFinder);
        printf("BUS#%d SDD#%d. Found %d candidate peaks to fit\n",iBUS,iSDD,nfound);
        thresholdPeakFinder = thresholdPeakFinder*.1; //initial = 0.01. It changes until it finds peaks
        nPFtries++;
      }
      if(nPFtries>=15) {
        cout<<"PEAK FINDER DOES NOT WORK, ==> CONTINUE"<<endl;
        continue;
      }

      xpeaks = spectrum->GetPositionX();  //array with X-positions of the centroids found by TSpectrum
      ypeaks = spectrum->GetPositionY();  //array with Y-positions of the centroids found by TSpectrum

      Double_t par[100];

      histoPF->Draw();

      // Estimate background
      fBkgnd[iBUS][iSDD] = new TF1(Form("fBkgnd_bus%d_sdd%d",iBUS,iSDD),"pol0(0)+expo(1)",1500,4500);
      hADCfit[iBUS][iSDD]->Fit(Form("fBkgnd_bus%d_sdd%d",iBUS,iSDD),"","",1500,4500);

      // Loop on all found peaks. Eliminate peaks at the background level
      par[0] = fBkgnd[iBUS][iSDD]->GetParameter(0);
      par[1] = fBkgnd[iBUS][iSDD]->GetParameter(1);
      par[2] = fBkgnd[iBUS][iSDD]->GetParameter(2);
      npeaks = 0;

      // Evaluation of the centroids for each single peak
      for (Int_t p=0; p<nfound; p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = hADCfit[iBUS][iSDD]->GetXaxis()->FindBin(xp);
        Double_t yp = hADCfit[iBUS][iSDD]->GetBinContent(bin);
        if (yp-TMath::Sqrt(yp) < fBkgnd[iBUS][iSDD]->Eval(xp)) continue;
        par[3*npeaks+3] = yp;
        par[3*npeaks+4] = xp;
        par[3*npeaks+5] = 15;
        npeaks++;
      }

      fPreFit[iBUS][iSDD] = new TF1(Form("fit_bus%d_sdd%d",iBUS,iSDD),findPeaks,1500,4500,3+3*npeaks);
      fPreFit[iBUS][iSDD]->SetParameters(par);
      fPreFit[iBUS][iSDD]->SetNpx(1000);
      //hADCfit[iBUS][iSDD]->Fit(Form("fit_bus%d_sdd%d",iBUS,iSDD));

      float themin = 999999.;
      int imin = 0;
      for(int i =0;i<nfound;i++){ //find the smallest, write it in xadc and remove it:
        themin = 999999.;
        for(int j =0;j<nfound;j++){
          if(xpeaks[j]<themin){
            themin = xpeaks[j];
            imin = j;
          }
        }
        xadc[i]=xpeaks[imin];
        yadc[i]=ypeaks[imin];
        xpeaks[imin]=999999.;
      }
      //print them, now ordered:
      for(int i =0;i<nfound;i++){cout<<"Found peak at ADC = "<<xadc[i]<<". Height = "<<yadc[i]<<endl;}
      //cout<<endl;








      //check if the peaks found are compatible with the assumption:
   //------------------------------------------------------------
   Float_t eDist10 = PFPeakE[1]-PFPeakE[0];
   Float_t eDist21 = PFPeakE[2]-PFPeakE[1];
   Float_t Erelation = eDist21/eDist10;
   //make a triad from all the peaks found:
   Float_t GPF  = 0.; //gain
   Float_t G0PF = 0.; //offset
   Bool_t TestPassed = false;
   Int_t ipeak0 = -1;
   Int_t ipeak1 = -1;
   Int_t ipeak2 = -1;
   for(Int_t i0=0;i0<nfound;i0++) {
     for(Int_t i1=i0+1;i1<nfound;i1++) {
       for(Int_t i2=i1+1;i2<nfound;i2++) {
         cout<<endl<<"-> Trying the triad: "<<i0<<" "<<i1<<" "<<i2<<endl;
         Float_t xDist10 = xadc[i1]-xadc[i0];
         Float_t xDist21 = xadc[i2]-xadc[i1];
         Float_t ADCrelation = xDist21/xDist10;
         cout<<"Checking assumption: Energy relation "<<Erelation<<" vs ADC relation "<<ADCrelation<<endl;
         //cout<<endl;

         //Define tolerance parameter
         Float_t tol = initTolerance; // 5%
         Bool_t TolerancePass = true;
         if(fabs(1.-(Erelation/ADCrelation))>tol) TolerancePass = false;
         if(!TolerancePass) cout<<"The test was not passed: ERROR!!!"<<endl;
         cout<<endl;
         //get the Peak Finder calibration offset GPFo and slope GPFs
         Float_t Dadc = xadc[i0]-xadc[i1];
         Float_t De   = PFPeakE[0]-PFPeakE[1];
         GPF  = De/Dadc;
         G0PF = -1.*xadc[i0]*GPF + PFPeakE[0];
         cout<<"Peak Finder: offset G0PF = "<<G0PF<<" and gain GPF = "<<GPF<<endl;

         //define an acceptable gain and offset; check if conditions (tolerance, G and G0) are met
         Float_t mingoodG   = 2.9;    Float_t maxgoodG  = 3.9;
         Float_t mingoodG0  = -3000;  Float_t maxgoodG0 = -1000;
         if(GPF<maxgoodG && GPF>mingoodG && G0PF<maxgoodG0 && G0PF>mingoodG0 && TolerancePass) {
           cout<<" -- TEST PASSED! --"<<endl;
           cout<<endl;
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

   new_file2.open("peaks.dat", ios::app);
   new_file2<<Form("BUS#%d SDD#%d has %d peaks",iBUS,iSDD,nfound)<<endl;
   new_file2.close();

   Float_t highestPeak = 0.;
   if(ipeak0>-1) {
     if(yadc[ipeak0]>highestPeak) highestPeak = yadc[ipeak0];
     if(yadc[ipeak1]>highestPeak) highestPeak = yadc[ipeak1];
     if(yadc[ipeak2]>highestPeak) highestPeak = yadc[ipeak2];
   }

   float diff;
   float dMin=0.97;
   float dMax=1.03;
   float peakTiKa, peakTiKb, peakCuKa, peakCuKb, peakMnKa, peakFeKa;
   for(int k=0;k<nfound;k++){
     diff = xadc[k]*GPF+G0PF;
     if (diff/eTiKa > dMin && diff/eTiKa < dMax ) {peakTiKa = xadc[k];} else
     if (diff/eTiKb > dMin && diff/eTiKb < dMax ) {peakTiKb = xadc[k];} else
     if (diff/eCuKa > dMin && diff/eCuKa < dMax ) {peakCuKa = xadc[k];} else
     if (diff/eCuKb > dMin && diff/eCuKb < dMax ) {peakCuKb = xadc[k];}
   }

   histoPF->Delete();



   fitFuncTotalN[iBUS][iSDD] = new TF1(Form("fitFuncTotal_%d_%d",iBUS,iSDD),"pol0(0)+expo(1)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",1500,4500);

   for (Int_t k = 0; k < 3; k++) {
     fitFuncTotalN[iBUS][iSDD]->SetParameter(k,fBkgnd[iBUS][iSDD]->GetParameter(k));
   }

   fitFuncTotalN[iBUS][iSDD]->SetParameter(4,peakTiKa);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(5,20);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(7,peakTiKb);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(8,20);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(10,peakCuKa);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(11,20);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(13,peakCuKb);
   fitFuncTotalN[iBUS][iSDD]->SetParameter(14,20);



   fitFuncTotalN[iBUS][iSDD]->SetNpx(1000);

   hADCfit[iBUS][iSDD]->Fit(fitFuncTotalN[iBUS][iSDD],"","",1500,4500);

   new_file.open("check.dat", ios::app);
   new_file<<Form("B#%d S#%d peakTiKa (%g): %g | %g",iBUS,iSDD,peakTiKa*GPF+G0PF,fitFuncTotalN[iBUS][iSDD]->GetParameter(1),peakTiKa)<<endl;
   new_file<<Form("B#%d S#%d peakTiKb (%g): %g | %g",iBUS,iSDD,peakTiKb*GPF+G0PF,fitFuncTotalN[iBUS][iSDD]->GetParameter(4),peakTiKb)<<endl;
   new_file<<Form("B#%d S#%d peakCuKa (%g): %g | %g",iBUS,iSDD,peakCuKa*GPF+G0PF,fitFuncTotalN[iBUS][iSDD]->GetParameter(7),peakCuKa)<<endl;
   new_file<<Form("B#%d S#%d peakCuKb (%g): %g | %g\n",iBUS,iSDD,peakCuKb*GPF+G0PF,fitFuncTotalN[iBUS][iSDD]->GetParameter(10),peakCuKb)<<endl;
   new_file.close();





 }

}

  ////////////////////////
  //
  myCanvas[0] = new TCanvas;
  //myCanvas[0]->SetLogy();
  //myCanvas[0]->SetGrid();

  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      if(!hADC[iBUS][iSDD] ) continue;
      hADC[iBUS][iSDD]->SetTitle(Form("BUS: %d, SDD: %d",iBUS,iSDD));
      hADC[iBUS][iSDD]->Draw();
      //myCanvas[0]->Print(Form("plots2/hADC_bus%d_sdd%d.png",iBUS,iSDD), "png");
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

      //fPreFit[iBUS][iSDD]->SetLineColor(4);
      //fPreFit[iBUS][iSDD]->SetLineWidth(2);
      //fPreFit[iBUS][iSDD]->Draw("same");

      fitFuncTotalN[iBUS][iSDD]->SetLineColor(2);
      fitFuncTotalN[iBUS][iSDD]->SetLineWidth(2);
      fitFuncTotalN[iBUS][iSDD]->Draw("same");

      //myCanvas[3]->Print(Form("plots2/hADC_Ylog_bus%d_sdd%d.png",iBUS,iSDD), "png");
    }
  }
/*
  myCanvas[1] = new TCanvas;
  //myCanvas[0]->SetLogy();
  //myCanvas[0]->SetGrid();

  for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
    for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
      if(!hADC[iBUS][iSDD]) continue;
      if(hADC[iBUS][iSDD]->GetEntries()<10000) continue;
      hADCpeaks[iBUS][iSDD]->SetTitle(Form("BUS: %d, SDD: %d",iBUS,iSDD));
      hADCpeaks[iBUS][iSDD]->Draw();

      fBkgnd->Draw("same");

      myCanvas[1]->Print(Form("plots2/hADCpeaks_bus%d_sdd%d.png",iBUS,iSDD), "png");
    }
  }
*/
  
  TCanvas* myCanvas[nBUS][4];

  for(Int_t k=1;k<nBUS;k++) {
    for(Int_t j=0;j<4;j++) {
      myCanvas[k][j] = new TCanvas;
      myCanvas[k][j]->Divide(4,4);  //column, line
      for(Int_t i=1;i<17;i++) {
        myCanvas[k][j]->cd(i);
        if(!hADC[k][i+16*j]) continue;
        hADC[k][i+16*j]->SetTitle(Form("BUS: %d, SDD: %d",k,i+16*j));
        hADC[k][i+16*j]->Draw();
      }
      myCanvas[k][j]->Print(Form("plots/hADC_Ylog_bus%d_sdd%d_%d.png",k,j*16,16*(j+1)), "png");
      myCanvas[k][j]->Print(Form("plots/hADC_Ylog_bus%d_sdd%d_%d.pdf",k,j*16,16*(j+1)), "pdf");
    }
  }


}
