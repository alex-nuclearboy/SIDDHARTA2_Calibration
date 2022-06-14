// Created by Francesco Sgaramella
// Modified 2022-06 by Aleksander Khreptak

#include  "TFile.h"
#include  "TTree.h"
#include  "TBranch.h"
#include  "TString.h"
#include  "TH2F.h"
#include  "TH2D.h"
#include  "TH1F.h"
#include  "TCanvas.h"
#include  "TSystemFile.h"
#include  "TSystemDirectory.h"
#include  "TLeaf.h"
#include  "stdio.h"
#include  <Riostream.h>

TTree *tr;
TH2F *h_chn_time[6][64];
TH1F *h_sdd[6][64], *h_crosstalk[6][64];
bool empty[6][64]; //, ISGOODtime;

void SDDhitmap(int sddnumber, int busnumber, int &column, int &row);
bool CrossTalkTiming(short drift, short drift_pre);
int SFERAnumber(int sdd);

void hist () {

    int nfile =0;
    int npath=1;
    char nbuff[100], nbuff1[100], nbuff2[100];
    bool ISGOODtime;
    int column, row;
    double runtime;

    TFile *f[nfile];
    TString path[npath];

    path[0] = "/home/nuclearboy/SIDDHARTA2/SIDDHARTA2_Calibration/rootfiles/SIDDHARTA2_xray/input/20220310/";

    for(int i=0; i<npath; i++){TSystemDirectory dir(path[i], path[i]);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file; TString fname; TIter next(files);
       while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          if (!file->IsDirectory() && fname.EndsWith(".root")) {
              f[nfile] = new TFile(path[i]+fname.Data(),"read");
              nfile++;
              cout << fname.Data() << endl;
            }
        }
    }}
    cout <<"number of root file: "<< nfile << endl;

    TFile *outfile= new TFile("/home/nuclearboy/SIDDHARTA2/SIDDHARTA2_Calibration/rootfiles/SIDDHARTA2_xray/output/20220310/hist_20220310_1900_0310_2000_xray_25kv_70ua_tube1_p1.root","RECREATE");

    const static Int_t MAXHITS = 40000;
    Int_t buf=0;
    Short_t kt[4] = {0};
    Int_t   nhits = 0 ;
    Int_t   bus[MAXHITS]   = {0};
    Short_t evnr[MAXHITS]  = {0};
    Short_t ht[MAXHITS]    = {0};
    Short_t trigg[MAXHITS] = {0};
    Short_t sdd[MAXHITS]   = {0};
    Short_t adc[MAXHITS]   = {0};
    Short_t drift[MAXHITS] = {0};
    Int_t   date = 0;
    Short_t ie   = 0;
    Short_t ip   = 0;
    Short_t dum  = 0;

    int ent = 0;

    for(int i=0; i<6; i++){
       for(int j=0; j<64; j++){
        sprintf(nbuff,"bus%d_sdd%d",i+1,j+1);
        sprintf(nbuff1,"bus%d_sdd%d_crosstalk",i+1,j+1);
        sprintf(nbuff2,"bus%d_sdd%d_time",i+1,j+1);
        h_sdd[i][j] = new TH1F(nbuff,nbuff,10000,0,10000);
        h_crosstalk[i][j] = new TH1F(nbuff1,nbuff,10000,0,10000);
        //h_chn_time[i][j] = new TH2F(nbuff2,nbuff2,10000,1634903017,1635162217,10000,0,10000);
        empty[i][j]=true;
       }
    }

    TH2F *sddmap = new TH2F ("sdd_map","sdd_position_map",48,1,49,8,1,9);//map of the sdd position around the target
    TH2D *sddrate = new TH2D ("sdd_rate","sdd_rate_map",48,1,49,8,1,9);//map of the sdd rate around the target
    TH2D *sddrate_sig = new TH2D ("sdd_rate_sig","sdd_rate_sig",48,1,49,8,1,9);//map of the sdd rate around the target
    TH2D *sddrate_noise = new TH2D ("sdd_rate_noise","sdd_rate_noise",48,1,49,8,1,9);//map of the sdd rate around the target

    for(Int_t j=0; j<nfile; j++){ //start the loop on the root input file
        f[j]->cd();
        tr = (TTree*)f[j]->Get("ft10");
        ent = (Int_t)tr->GetEntries();

        tr->GetEntry(0);
        int timestart=(int)tr->GetBranch("date")->GetLeaf("date")->GetValue(0);
        tr->GetEntry(ent-1);
        int timeend=(int)tr->GetBranch("date")->GetLeaf("date")->GetValue(0);
        int DeltaTime= timeend-timestart;
        runtime=runtime+DeltaTime;
        cout << runtime << endl;

        tr->SetBranchAddress("nhits",&nhits);
        tr->SetBranchAddress("evnr",&evnr);
        tr->SetBranchAddress("sdd",&sdd);
        tr->SetBranchAddress("adc",&adc);
        tr->SetBranchAddress("drift",&drift);
        tr->SetBranchAddress("date",&date);
        tr->SetBranchAddress("bus",&bus);
        tr->SetBranchAddress("ie",&ie);
        tr->SetBranchAddress("ip",&ip);
        tr->SetBranchAddress("dum",&dum);
        tr->SetBranchAddress("kt",&kt);
        tr->SetBranchAddress("trigg",&trigg);
        tr->SetBranchAddress("buf",&buf);
        tr->SetBranchAddress("ht",&ht);

        for( Int_t ient = 0; ient < ent; ient ++ ){
            if(ient%100==0) {cout<<"reading file "<< j+1 << " entry "<< ient <<"/"<< ent << endl;}
            //if(ient ==2513) continue;
            tr->GetEntry(ient);

            for( Int_t ihit = 0; ihit < nhits; ihit ++ ){

                ISGOODtime=true;
                if(ihit>0){ISGOODtime = CrossTalkTiming(drift[ihit],drift[ihit-1]);}
                SDDhitmap(sdd[ihit], bus[ihit], column, row);
                sddmap->SetBinContent(column,row,sdd[ihit]);
                sddrate->Fill(column,row);
                if(adc[ihit]>1200){sddrate_sig->Fill(column,row);}
                if(adc[ihit]<=1200){sddrate_noise->Fill(column,row);}

                if(ISGOODtime){h_sdd[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]); empty[bus[ihit]-1][sdd[ihit]-1]=false;}
                if(!ISGOODtime){h_crosstalk[bus[ihit]-1][sdd[ihit]-1]->Fill(adc[ihit]);}
                //h_chn_time[bus[ihit]-1][sdd[ihit]-1]->Fill(date,adc[ihit]);
            }
        }
    }

    sddrate->Scale(1/runtime); sddrate_sig->Scale(1/runtime); sddrate_noise->Scale(1/runtime);
        for(int i=1;i<=49;i++){
        for(int j=1;j<9;j++){
            if(sddrate->GetBinContent(i,j)<0.001){
                sddrate->SetBinContent(i,j,0);
                sddrate_sig->SetBinContent(i,j,0);
                sddrate_noise->SetBinContent(i,j,0);
            }
        }
    }


    outfile->cd();
    for(Int_t i=0; i<6; i++){
       for(Int_t j=0; j<64; j++){
            if(h_sdd[i][j]->GetEntries()>0){
                h_sdd[i][j]->Write();
                h_crosstalk[i][j]->Write();
                //h_chn_time[i][j]->Write();
            }
        }
    }
    sddmap->Write();
    sddrate->Write();
    sddrate_sig->Write();
    sddrate_noise->Write();

    outfile->Close();

}


void SDDhitmap(int sddnumber, int busnumber, int &column, int &row) {
    row=0; column=0;
    int SFERA=0;
    //if(busnumber==5){busnumber=1;}
    //else if(busnumber==1){busnumber=2;}
    //else if(busnumber==3){busnumber=4;}
    //else if(busnumber==2){busnumber=5;}
    //else{cout << "busnumber unknown: bus" << busnumber << endl; getchar();}

    //Back side view
    if(sddnumber<=16){sddnumber=sddnumber;SFERA=0;}
    if(sddnumber>=17 && sddnumber<=32){sddnumber=sddnumber-16;SFERA=2;}
    if(sddnumber>=33 && sddnumber<=48){sddnumber=sddnumber-(16*2);SFERA=4;}
    if(sddnumber>=49 && sddnumber<=64){sddnumber=sddnumber-(16*3);SFERA=6;}

    if(sddnumber==1){row=1; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber== 9){row=5; column=1+SFERA+((busnumber-1)*8);}
    if(sddnumber==2){row=4; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber==10){row=8; column=1+SFERA+((busnumber-1)*8);}
    if(sddnumber==3){row=2; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber==11){row=7; column=1+SFERA+((busnumber-1)*8);}
    if(sddnumber==4){row=3; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber==12){row=6; column=1+SFERA+((busnumber-1)*8);}
    if(sddnumber==5){row=3; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==13){row=6; column=2+SFERA+((busnumber-1)*8);}
    if(sddnumber==6){row=2; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==14){row=7; column=2+SFERA+((busnumber-1)*8);}
    if(sddnumber==7){row=1; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==15){row=5; column=2+SFERA+((busnumber-1)*8);}
    if(sddnumber==8){row=4; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==16){row=8; column=2+SFERA+((busnumber-1)*8);}
}

bool CrossTalkTiming(short drift, short drift_pre){
    bool ISGOODtime = true;
    double t1=0, t2=0, timediff=0;
    t1=drift+32768;
    t2=drift_pre+32768;
    if(t1>t2) {timediff = t1-t2;}
    if(t2>t1) {timediff = t1+(32768.*2)-t2;}
    if((timediff>0. && timediff<625)) {ISGOODtime = false;}//625 = 5 microseconds
    return ISGOODtime;
}

int SFERAnumber(int sdd){
    int SFERA=0;
    if(sdd<=16){SFERA=4;}
    if(sdd>=17 && sdd<=32){SFERA=3;}
    if(sdd>=33 && sdd<=48){SFERA=2;}
    if(sdd>=49 && sdd<=64){SFERA=1;}

    return SFERA;
}
