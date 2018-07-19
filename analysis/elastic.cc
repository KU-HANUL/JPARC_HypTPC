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
#include "hodo_elastic.h"

using namespace std;

void elastic()
{
  //gStyle -> SetOptFit(1000);

  bool pdf_saving = true;
  //bool pdf_saving = false;

  hodo_elastic hodo_class;

  TString pi_plus_pip="../rootfile/e45/elastic/pi_plus/pipip/";
  TString pi_minus_pip="../rootfile/e45/elastic/pi_minus/pipip/";

  TString filename[8]={"p635_phsp.root","p835_phsp.root","p1035_phsp.root","p1235_phsp.root","p1435_phsp.root","p1635_phsp.root","p1835_phsp.root","p2000_phsp.root" };
  double momentum[8]={0.635,0.835,1.035,1.235,1.435,1.635,1.835,2.000};

  TFile *Delta_pip[8];
  TFile *N_pip[8];
  for(int i=7;i<8;i++){
    Delta_pip[i] = new TFile(Form("%s%s",pi_plus_pip.Data(),filename[i].Data()),"READ");
    N_pip[i] = new TFile(Form("%s%s",pi_minus_pip.Data(),filename[i].Data()),"READ");
  }


  double ecut=1.0;
  //double ecut=0;
  double test;
  hodo_class.pip(Delta_pip[7],ecut,false,test,true);
  //hodo_class.pip(Delta_pip[7],ecut,true,test,true);
  //hodo_class.pip(N_pip[7],ecut,false,test,true);
  //hodo_class.monitor(Delta_pip[7],0.1,false);

  /*
  double prob_Delta_pip_wo[8];
  double prob_Delta_pip_w[8];
  for(int i=0;i<8;i++){
  hodo_class.pip(Delta_pip[i],ecut,false,prob_Delta_pip_wo[i],false);
    hodo_class.pip(Delta_pip[i],ecut,true,prob_Delta_pip_w[i],false);
  }

  double prob_N_pip_wo[8];
  double prob_N_pip_w[8];
  for(int i=0;i<8;i++){
    hodo_class.pip(N_pip[i],ecut,false,prob_N_pip_wo[i],false);
    hodo_class.pip(N_pip[i],ecut,true,prob_N_pip_w[i],false);
  }

  for(int i=0;i<8;i++){
    //
  }

  TGraphErrors *gr_N_pip_wo = new TGraphErrors(8,momentum,prob_N_pip_wo);
  gr_N_pip_wo -> SetMarkerColor(2);
  gr_N_pip_wo -> SetMarkerStyle(20);
  gr_N_pip_wo -> SetMarkerSize(1.5);

  TGraphErrors *gr_N_pip_w = new TGraphErrors(8,momentum,prob_N_pip_w);
  gr_N_pip_w -> SetMarkerColor(2);
  gr_N_pip_w -> SetMarkerStyle(21);
  gr_N_pip_w -> SetMarkerSize(1.5);

  TGraphErrors *gr_Delta_pip_wo = new TGraphErrors(8,momentum,prob_Delta_pip_wo);
  gr_Delta_pip_wo -> SetMarkerColor(4);
  gr_Delta_pip_wo -> SetMarkerStyle(20);
  gr_Delta_pip_wo -> SetMarkerSize(1.5);

  TGraphErrors *gr_Delta_pip_w = new TGraphErrors(8,momentum,prob_Delta_pip_w);
  gr_Delta_pip_w -> SetMarkerColor(4);
  gr_Delta_pip_w -> SetMarkerStyle(21);
  gr_Delta_pip_w -> SetMarkerSize(1.5);

  TCanvas *c1 = new TCanvas("c1","c1",1200, 1200);
  TMultiGraph *mg_pi_minus = new TMultiGraph();
  mg_pi_minus -> SetTitle("M-1 coin./M-1 events; Beam Momentum (GeV/c) ; %");
  mg_pi_minus -> Add(gr_N_pip_wo);
  mg_pi_minus -> Add(gr_N_pip_w);

  TLegend *legend_pi_minus = new TLegend(0.6,0.6,0.8,0.9);
  gStyle -> SetLegendFont(42);
  gStyle -> SetLegendTextSize(0.03);
  legend_pi_minus -> AddEntry((TObject*)0,"Without window","");
  legend_pi_minus -> AddEntry(gr_N_pip_wo,"#pi^{-}p #rightarrow #pi^{-}p","p");
  legend_pi_minus -> AddEntry((TObject*)0,"With window","");
  legend_pi_minus -> AddEntry(gr_N_pip_w,"#pi^{-}p #rightarrow #pi^{-}p","p");

  c1 -> cd();
  mg_pi_minus -> Draw("AP");
  //mg_pi_minus -> GetHistogram() -> GetYaxis() -> SetRangeUser(0.2,0.8);
  legend_pi_minus -> Draw("same");


  TCanvas *c2 = new TCanvas("c2","c2",1200, 1200);
  TMultiGraph *mg_pi_plus = new TMultiGraph();
  mg_pi_plus -> SetTitle("M-1 coin./M-1 events; Beam Momentum (GeV/c) ; %");
  mg_pi_plus -> Add(gr_Delta_pip_wo);
  mg_pi_plus -> Add(gr_Delta_pip_w);

  TLegend *legend_pi_plus = new TLegend(0.6,0.6,0.8,0.9);
  gStyle -> SetLegendFont(42);
  gStyle -> SetLegendTextSize(0.03);
  legend_pi_plus -> AddEntry((TObject*)0,"Without window","");
  legend_pi_plus -> AddEntry(gr_Delta_pip_wo,"#pi^{+}p #rightarrow #pi^{+}p","p");
  legend_pi_plus -> AddEntry((TObject*)0,"With window","");
  legend_pi_plus -> AddEntry(gr_Delta_pip_w,"#pi^{+}p #rightarrow #pi^{+}p","p");

  c2 -> cd();
  mg_pi_plus -> Draw("AP");
  //mg_pi_plus -> GetHistogram() -> GetYaxis() -> SetRangeUser(0.2,0.8);
  legend_pi_plus -> Draw("same");
  */
}//end
