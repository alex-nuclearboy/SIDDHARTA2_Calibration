#ifndef _FITTING_FUNC_H
#define _FITTING_FUNC_H

#include  <TMath.h>

const Double_t PI       = TMath::Pi();
const Double_t EPSILON  = 3.71; //electron-hole pair creation energy

//Background: a constant function + an exponential function
Double_t  fitBkgndFunc(Double_t *x, Double_t *par) {
  return  par[0]+par[1]*x[0] + TMath::Exp(par[2]+par[3]*x[0]);
}

//Gaussian function
Double_t fitGaussFunc(Double_t *x, Double_t *par) {
  Double_t  sigma     = TMath::Sqrt( par[1]*EPSILON*par[2] + TMath::Power(par[3]/2.35,2) );
  Double_t  amplitude = par[0]/(TMath::Sqrt(2*PI)*sigma);
  Double_t  arg       = -0.5*TMath::Power(x[0]-par[2],2)/TMath::Power(sigma,2);
  return    amplitude*TMath::Exp(arg);
}

//Tail function
Double_t fitTailFunc(Double_t *x, Double_t *par) {
  Double_t  sigma      = TMath::Sqrt( par[1]*EPSILON*par[2] + TMath::Power(par[3]/2.35,2) );
  Double_t  amplitude  = par[0]/(2*sigma*par[4]);
  Double_t  arg        = (x[0]-par[2])/(sigma*par[4]) + 1/(2*TMath::Power(par[4],2));
  Double_t  errorFunc  = (x[0]-par[2])/(TMath::Sqrt(2)*sigma) + 1/(TMath::Sqrt(2)*par[4]);
  return    amplitude*TMath::Exp(arg)*TMath::Erfc(errorFunc);
}

// Gaussian + Tail functions
Double_t fitGaussTailFunc(double *x, double *par) {
  Double_t  sigma      = TMath::Sqrt( par[1]*EPSILON*par[2] + TMath::Power(par[3]/2.35,2) );
  Double_t  amplitudeG = par[0]/(TMath::Sqrt(2*PI)*sigma);
  Double_t  amplitudeT = par[0]/(2*sigma*par[5]);
  Double_t  argG       = -0.5*TMath::Power(x[0]-par[2],2)/TMath::Power(sigma,2);
  Double_t  argT       = (x[0]-par[2])/(sigma*par[5]) + 1/(2*TMath::Power(par[5],2));
  Double_t  errorFunc  = (x[0]-par[2])/(TMath::Sqrt(2)*sigma) + 1/(TMath::Sqrt(2)*par[5]);
  return    amplitudeG*TMath::Exp(argG) + amplitudeT*par[4]*TMath::Exp(argT)*TMath::Erfc(errorFunc);
}

// Global function
Double_t fitTotalFunc(Double_t *x, Double_t *par) {
      return  fitGaussTailFunc(x,par) + fitGaussTailFunc(x,&par[6]) +         // Ti
              fitGaussFunc(x,&par[12]) + fitGaussFunc(x,&par[16]) +           // Mn & Fe
              fitGaussTailFunc(x,&par[20]) + fitGaussTailFunc(x,&par[26]) +   // Cu
              fitBkgndFunc(x,&par[32]);                                       // background
}

#endif
