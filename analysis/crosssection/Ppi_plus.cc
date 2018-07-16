#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <TPad.h>
#include <cmath>
#include "Math/IFunction.h"
#include "TSystem.h"

void Ppi_plus()
{

  TFile *file = new TFile("HEPData-ins97248-v1-root.root","READ");
  TFile *tfile = new TFile("test.root","recreate");

  TDirectory *dir[16];
  TGraphAsymmErrors *gr[16];

  TCanvas *can1 = new TCanvas("can1","",2000,1200);
  can1 -> Divide(4,4);
  TCanvas *can2 = new TCanvas("can2","",2000,1200);
  TCanvas *can3 = new TCanvas("can3","",2000,1200);

  for(int i=0;i<16;i++){
    dir[i] = (TDirectory*)file->GetDirectory(Form("Table %d",i+1));
    dir[i] -> GetObject("Graph1D_y1",gr[i]);
    //dir[i] -> ls();
    can1 -> cd(i+1);
    gr[i] -> Draw();

  }

  TDirectory *dir_amp;
  TGraphAsymmErrors *gr_amp[16];
  TCanvas *can_amp = new TCanvas("can_amp","",2000,1200);
  can_amp -> Divide(4,4);
  dir_amp = (TDirectory*)file->GetDirectory("Table 17");
  dir_amp -> ls();

  double error_high[16][9]={0};
  double error_low[16][9]={0};
  double P[16][9]={0};
  double amplitude[16][9]={0};
  int start;

  for(int i=0;i<9;i++){
    dir_amp -> GetObject(Form("Graph1D_y%d",i+1),gr_amp[i]);
    can_amp -> cd(i+1);
    gr_amp[i] -> Draw();
    if(i==6||i==7) start=1;
    else if(i==8) start=2;
    else start=0;
    for(int j=start;j<16;j++){
      gr_amp[i] -> GetPoint(j-start,P[j][i],amplitude[j][i]);
      error_high[j][i] = gr_amp[i] -> GetErrorXhigh(j);
      error_low[j][i] = gr_amp[i] -> GetErrorXlow(j);
    }
  }


  TCanvas *can_legen = new TCanvas("can_legen","",2000,1200);
  can_legen -> Divide(4,4);

  TF1 *legen[10];
  legen[0] = new TF1("legen[0]","1",-1,1);
  legen[1] = new TF1("legen[1]","x",-1,1);
  legen[2] = new TF1("legen[2]","(3*TMath::Power(x,2)-1)/2 ",-1,1);
  legen[3] = new TF1("legen[3]","(5*TMath::Power(x,3)-3*x)/2 ",-1,1);
  legen[4] = new TF1("legen[4]","(35*TMath::Power(x,4)-30*TMath::Power(x,2)+3)/8 ",-1,1);
  legen[5] = new TF1("legen[5]","(63*TMath::Power(x,5)-70*TMath::Power(x,3)+15*x)/8 ",-1,1);
  legen[6] = new TF1("legen[6]","(231*TMath::Power(x,6)-315*TMath::Power(x,4)+105*TMath::Power(x,2)-5)/16 ",-1,1);
  legen[7] = new TF1("legen[7]","(429*TMath::Power(x,7)-693*TMath::Power(x,5)+315*TMath::Power(x,3)-35*x)/16 ",-1,1);
  legen[8] = new TF1("legen[8]","(6435*TMath::Power(x,8)-12012*TMath::Power(x,6)+6930*TMath::Power(x,4)-1260*TMath::Power(x,2)+35)/128 ",-1,1);
  legen[9] = new TF1("legen[9]","(12155*TMath::Power(x,9)-25740*TMath::Power(x,7)+18018*TMath::Power(x,5)-4620*TMath::Power(x,3)+315*x)/128 ",-1,1);


  TF1 *diff[16];
  for(int i=0;i<16;i++){
    diff[i]=  new TF1(Form("diff_%d",i),Form("legen[0]*%f+legen[1]*%f+legen[2]*%f+legen[3]*%f+legen[4]*%f+legen[5]*%f+legen[6]*%f+legen[7]*%f+legen[8]*%f+legen[9]*%f",1.0,amplitude[i][0],amplitude[i][1],amplitude[i][2],amplitude[i][3],amplitude[i][4],amplitude[i][5],amplitude[i][6],amplitude[i][7],amplitude[i][8]),-0.9,0.9);
  }

  const int ch=11;
  for(int i=0;i<16;i++){
    can_legen -> cd(i+1);
    diff[i] -> Draw();
    if(i==ch){
      can2 -> cd();
      diff[ch] -> Draw();
      diff[ch] -> GetHistogram() -> SetTitle("Differential cross section");
      diff[ch] -> GetHistogram() -> GetXaxis()-> SetTitle("Cos#theta_{C.M.}");
      diff[ch] -> GetHistogram() -> GetYaxis()-> SetTitle("d#sigma/d(Cos#theta)(arb.)");
      can3 -> cd();
      gr[ch] -> Draw();
      gr[ch] -> GetHistogram() -> SetTitle("Differential cross section");
      gr[ch] -> GetHistogram() -> GetXaxis()-> SetTitle("Cos#theta_{C.M.}");
      gr[ch] -> GetHistogram() -> GetYaxis()-> SetTitle("d#sigma/d(Cos#theta)(arb.)");
    }
  }

  std::cout<<"max"<<diff[ch]->GetMaximum()<<std::endl;
  for(int j=0;j<9;j++){
    //std::cout<<"amplitude : "<<" "<<amplitude[11][j]<<std::endl;
    //std::cout<<"P : "<<" "<<P[11][j]<<std::endl;
    std::cout<<", "<<amplitude[ch][j]<<std::endl;
  }

}//end
