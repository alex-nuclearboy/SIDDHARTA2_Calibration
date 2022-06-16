#ifndef _CALIB_SDD_H
#define _CALIB_SDD_H

#include  <TFile.h>
#include  <TGraph.h>
#include  <TGraphErrors.h>
#include  <TClonesArray.h>
#include  <TF1.h>
#include  <TH1.h>
#include  <TH2.h>
#include  <TH1F.h>
#include  <TH2F.h>
#include  <TH1D.h>
#include  <TH2D.h>
#include  <TCanvas.h>
#include  <TLegend.h>
#include  <TPaveText.h>
#include  <TStyle.h>

TFile* f[2];

Int_t busNumber,sddNumber;

Int_t beginWindow[10],endWindow[10];
Int_t binPeak[10];
Double_t peakValue[10];

TH1F  *hSDD,*hSDD_copy,*hSDD_fit;
TH1F  *hCrosstalk;
TH1F  *hEnergySDD,*hEnergySDD_fit;

TGraph  *gLinearity;
TGraph  *gSDDenergy;

TF1   *fitFuncBkgrndStart,*fitFuncBkgrnd;

TF1   *fitFuncGaussTiKalphaStart,*fitFuncGaussTiKalpha,*fitFuncGaussTiKbetaStart,*fitFuncGaussTiKbeta;
TF1   *fitFuncGaussMnKalphaStart,*fitFuncGaussMnKalpha;
TF1   *fitFuncGaussFeKalphaStart,*fitFuncGaussFeKalpha;
TF1   *fitFuncGaussCuKalphaStart,*fitFuncGaussCuKalpha,*fitFuncGaussCuKbetaStart,*fitFuncGaussCuKbeta;
TF1   *fitFuncTailTiKalphaStart,*fitFuncTailTiKalpha,*fitFuncTailTiKbetaStart,*fitFuncTailTiKbeta;
TF1   *fitFuncTailCuKalphaStart,*fitFuncTailCuKalpha,*fitFuncTailCuKbetaStart,*fitFuncTailCuKbeta;

TF1   *fitFuncGaussTailTiKalphaStart,*fitFuncGaussTailTiKalpha,*fitFuncGaussTailTiKbetaStart,*fitFuncGaussTailTiKbeta;
TF1   *fitFuncGaussTailCuKalphaStart,*fitFuncGaussTailCuKalpha,*fitFuncGaussTailCuKbetaStart,*fitFuncGaussTailCuKbeta;

TF1   *fitFuncTotal;
TF1   *fitFuncEnergyTotal;

TCanvas* myCanvas[10];
TLegend* myLegend[10];

#endif
