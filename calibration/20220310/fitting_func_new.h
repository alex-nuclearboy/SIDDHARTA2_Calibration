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
    
    //Double_t fr      = par[8];
    //Double_t beta    = par[9];
    Double_t FF      = par[24];
    Double_t noise   = par[25];
    Double_t  sigma[6];
    sigma[0] = TMath::Sqrt( FF*EPSILON*par[1] + TMath::Power(noise/2.35,2) );
    sigma[1] = TMath::Sqrt( FF*EPSILON*par[5] + TMath::Power(noise/2.35,2) );
    sigma[2] = TMath::Sqrt( FF*EPSILON*par[9] + TMath::Power(noise/2.35,2) );
    sigma[3] = TMath::Sqrt( FF*EPSILON*par[13] + TMath::Power(noise/2.35,2) );
    sigma[4] = TMath::Sqrt( FF*EPSILON*par[17] + TMath::Power(noise/2.35,2) );
    sigma[5] = TMath::Sqrt( FF*EPSILON*par[21] + TMath::Power(noise/2.35,2) );
    Double_t argG[6];
    Double_t argT[6];
    argG[0] = -0.5*TMath::Power(x[0]-par[1],2)/TMath::Power(sigma[0],2);
    argT[0] = (x[0]-par[1])/(sigma[0]*par[3]) + 1/(2*TMath::Power(par[3],2));
    argG[1] = -0.5*TMath::Power(x[0]-par[5],2)/TMath::Power(sigma[1],2);
    argT[1] = (x[0]-par[5])/(sigma[1]*par[7]) + 1/(2*TMath::Power(par[7],2));
    argG[2] = -0.5*TMath::Power(x[0]-par[9],2)/TMath::Power(sigma[2],2);
    argT[2] = (x[0]-par[9])/(sigma[2]*par[11]) + 1/(2*TMath::Power(par[11],2));
    argG[3] = -0.5*TMath::Power(x[0]-par[13],2)/TMath::Power(sigma[3],2);
    argT[3] = (x[0]-par[13])/(sigma[3]*par[15]) + 1/(2*TMath::Power(par[15],2));
    argG[4] = -0.5*TMath::Power(x[0]-par[17],2)/TMath::Power(sigma[4],2);
    argT[4] = (x[0]-par[17])/(sigma[4]*par[19]) + 1/(2*TMath::Power(par[19],2));
    argG[5] = -0.5*TMath::Power(x[0]-par[21],2)/TMath::Power(sigma[5],2);
    argT[5] = (x[0]-par[21])/(sigma[5]*par[23]) + 1/(2*TMath::Power(par[23],2));
    Double_t errorFunc[6];
    errorFunc[0] = (x[0]-par[1])/(TMath::Sqrt(2)*sigma[0]) + 1/(TMath::Sqrt(2)*par[3]);
    errorFunc[1] = (x[0]-par[5])/(TMath::Sqrt(2)*sigma[1]) + 1/(TMath::Sqrt(2)*par[7]);
    errorFunc[2] = (x[0]-par[9])/(TMath::Sqrt(2)*sigma[2]) + 1/(TMath::Sqrt(2)*par[11]);
    errorFunc[3] = (x[0]-par[13])/(TMath::Sqrt(2)*sigma[3]) + 1/(TMath::Sqrt(2)*par[15]);
    errorFunc[4] = (x[0]-par[17])/(TMath::Sqrt(2)*sigma[4]) + 1/(TMath::Sqrt(2)*par[19]);
    errorFunc[5] = (x[0]-par[21])/(TMath::Sqrt(2)*sigma[5]) + 1/(TMath::Sqrt(2)*par[23]);
    
    return  (par[0]/(TMath::Sqrt(2*PI)*sigma[0]))*TMath::Exp(argG[0]) + (par[2]*par[0]/(2*sigma[0]*par[3]))*TMath::Exp(argT[0])*TMath::Erfc(errorFunc[0]) + 
            (par[4]/(TMath::Sqrt(2*PI)*sigma[1]))*TMath::Exp(argG[1]) + (par[6]*par[4]/(2*sigma[1]*par[7]))*TMath::Exp(argT[1])*TMath::Erfc(errorFunc[1]) + 
            (par[8]/(TMath::Sqrt(2*PI)*sigma[2]))*TMath::Exp(argG[2]) + (par[10]*par[8]/(2*sigma[2]*par[11]))*TMath::Exp(argT[2])*TMath::Erfc(errorFunc[2]) + 
            (par[12]/(TMath::Sqrt(2*PI)*sigma[3]))*TMath::Exp(argG[3]) + (par[14]*par[12]/(2*sigma[3]*par[15]))*TMath::Exp(argT[3])*TMath::Erfc(errorFunc[3]) +
            (par[16]/(TMath::Sqrt(2*PI)*sigma[4]))*TMath::Exp(argG[4]) + (par[18]*par[16]/(2*sigma[4]*par[19]))*TMath::Exp(argT[4])*TMath::Erfc(errorFunc[4]) + 
            (par[20]/(TMath::Sqrt(2*PI)*sigma[5]))*TMath::Exp(argG[5]) + (par[22]*par[20]/(2*sigma[5]*par[23]))*TMath::Exp(argT[5])*TMath::Erfc(errorFunc[5]) +
            //par[18]*TMath::Gaus(x[0],par[19],par[20]) +
            //par[21]*TMath::Gaus(x[0],par[22],par[23]) +
            fitBkgndFunc(x,&par[26]);
 
}

#endif
