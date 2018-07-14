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

void test()
{

  TFile *file = new TFile("HEPData-ins97248-v1-root.root","READ");
  TFile *tfile = new TFile("test.root","recreate");

  TDirectory *dir[16];
  TGraph *gr[16];
  TGraph *gr1[16];
  TGraph *gr2[16];
  /*
    TGraphAsymmErrors *gr[16];
    TGraphAsymmErrors *gr1[16];
    TGraphAsymmErrors *gr2[16];
  */
  TCanvas *can1 = new TCanvas("can1","",1200,1200);
  can1 -> Divide(4,4);
  TCanvas *can2 = new TCanvas("can2","",1200,1200);
  can2 -> Divide(4,4);
  TCanvas *can3 = new TCanvas("can3","",1200,1200);
  can3 -> Divide(4,4);

  TF1 *extra1 = new TF1("extra1","[0]*x*x+[1]*x+[2]",0.8,1.0);
  TF1 *extra2 = new TF1("extra2","[0]*x*x+[1]*x+[2]",-1.0,-0.8);

  //for(int i=0;i<16;i++){
  for(int i=11;i<12;i++){

    dir[i] = (TDirectory*)file->GetDirectory(Form("Table %d",i+1));
    dir[i] -> GetObject("Graph1D_y1",gr1[i]);
    dir[i] -> GetObject("Graph1D_y1",gr2[i]);
    //dir[i] -> ls();

    int number=gr1[i] -> GetN();
    //gr[i] = new TGraph(number+2);
    gr[i] = new TGraph();
    can1 -> cd(i+1);
    gr1[i] -> Draw();
    gr1[i] -> Fit(extra1,"REQ");

    can2 -> cd(i+1);
    gr2[i] -> Draw();
    gr2[i] -> Fit(extra2,"REQ");

    double dummy_x=2, dummy_y=2;
    double dummy_cos1, dummy_cos2;
    double max=extra1->Eval(1.0);
    gr[i]->SetPoint(0,1.0,1.0);
    for(int j=0;j<number;j++){
      gr1[i]->GetPoint(j,dummy_x,dummy_y);
      gr[i]->SetPoint(j+1,dummy_x,dummy_y/max);
    }
    gr[i]->SetPoint(number+1,-1.0,(extra2->Eval(-1.0))/max);
    can3 -> cd(i+1);
    gr[i] -> Draw();
    double x,y;
    std::cout<<"x"<<std::endl;
    for(int j=0;j<number+2;j++){
      gr[i]->GetPoint(j,x,y);
      std::cout<<" ,"<<x<<std::endl;
    }
    std::cout<<"y"<<std::endl;
    for(int j=0;j<number+2;j++){
      gr[i]->GetPoint(j,x,y);
      std::cout<<" ,"<<y<<std::endl;
    }

    gr[i] -> SetTitle(Form("gr_%d",i));
  }
  tfile -> cd();
  gr[11] -> Write(Form("gr_%d",11));
  /*
  TDirectory *dir_amp;
  TGraphAsymmErrors *gr_amp[16];
  //TCanvas *can_amp = new TCanvas("can_amp","",1200,1200);
  //can_amp -> Divide(4,4);
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


  TCanvas *can_legen = new TCanvas("can_legen","",1200,1200);
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
    //diff[i]=  new TF1(Form("diff_%d",i),Form("TMath::Power(legen[0]*%f+legen[1]*%f+legen[2]*%f+legen[3]*%f+legen[4]*%f+legen[5]*%f+legen[6]*%f+legen[7]*%f+legen[8]*%f+legen[9]*%f,2)",1.0,amplitude[0][i],amplitude[1][i],amplitude[2][i],amplitude[3][i],amplitude[4][i],amplitude[5][i],amplitude[6][i],amplitude[7][i],amplitude[8][i]),-0.9,0.9);

    diff[i] = new TF1(Form("diff_%d",i),Form("TMath::Power(legen[0]*%f+legen[1]*%f+legen[2]*%f+legen[3]*%f+legen[4]*%f+legen[5]*%f+legen[6]*%f+legen[7]*%f+legen[8]*%f+legen[9]*%f,2)",1.0,amplitude[0][i]*3,amplitude[1][i]*5,amplitude[2][i]*7,amplitude[3][i]*9,amplitude[4][i]*11,amplitude[5][i]*13,amplitude[6][i]*15,amplitude[7][i]*17,amplitude[8][i]*19),-0.9,0.9);

  }

  for(int i=0;i<16;i++){
    can_legen -> cd(i+1);
    legen[i] -> Draw();
    //diff[i] -> Draw();
  }

  for(int i=0;i<16;i++){
  std::cout<<"  "<<std::endl;
    for(int j=0;j<9;j++){
      std::cout<<"amplitude : "<<i <<" "<<amplitude[i][j]<<std::endl;
    }
  }
*/
}//end
