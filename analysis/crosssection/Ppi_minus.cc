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

void Ppi_minus()
{

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

  TF1 *diff = new TF1("diff","legen[0]*[0]+legen[1]*[1]+legen[2]*[2]+legen[3]*[3]+legen[4]*[4]+legen[5]*[5]+legen[6]*[6]+legen[7]*[7]+legen[8]*[8]+legen[9]*[9]",-1.0,1.0);
  TF1 *diff2 = new TF1("diff2","legen[0]*[0]+legen[1]*[1]+legen[2]*[2]+legen[3]*[3]+legen[4]*[4]+legen[5]*[5]+legen[6]*[6]+legen[7]*[7]+legen[8]*[8]+legen[9]*[9]",-1.0,1.0);

  TFile *file = new TFile("HEPData-ins1104030-v1-root.root","READ");
  TFile *tfile = new TFile("test.root","recreate");

  TDirectory *dir[31];
  TGraphAsymmErrors *gr[31];

  TCanvas *can1 = new TCanvas("can1","",1200,1200);
  can1 -> Divide(6,6);

  TCanvas *can2 = new TCanvas("can2","",1200,1200);
  TCanvas *can3 = new TCanvas("can3","",1200,1200);

  double amplitude[31][10]={0};
  double max[31]={0};
  const int ch=18;
  for(int i=0;i<31;i++){
    dir[i] = (TDirectory*)file->GetDirectory(Form("Table %d",i+1));
    dir[i] -> GetObject("Graph1D_y1",gr[i]);
    //dir[i] -> ls();
    can1 -> cd(i+1);
    gr[i] -> Draw();
    gr[i] -> Fit(diff,"REQ");
    diff -> GetParameters(amplitude[i]);
    max[i] = diff -> GetMaximum();
    if(i==ch){
      can2 -> cd();
      diff2 -> SetParameters(amplitude[i]);
      diff2 -> Draw();
      diff2 -> GetHistogram() -> SetTitle("Differential cross section");
      diff2 -> GetHistogram() -> GetXaxis()-> SetTitle("Cos#theta_{C.M.}");
      diff2 -> GetHistogram() -> GetYaxis()-> SetTitle("d#sigma/d(Cos#theta)(arb.)");
      can3 -> cd();
      gr[ch] -> Draw();
      gr[ch] -> GetHistogram() -> SetTitle("Differential cross section");
      gr[ch] -> GetHistogram() -> GetXaxis()-> SetTitle("Cos#theta_{C.M.}");
      gr[ch] -> GetHistogram() -> GetYaxis()-> SetTitle("d#sigma/d#Omega");
    }
  }

  std::cout<<"max"<<max[ch]<<std::endl;
  for(int j=0;j<10;j++){
    std::cout<<", "<<amplitude[ch][j]<<std::endl;
  }

}//end
