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
#include "hodo.h"

using namespace std;

void fifc()
{
  //gStyle -> SetOptFit(1000);

  bool pdf_saving = true;
  //bool pdf_saving = false;

  hodo hodo_class;

  TString pi_plus_pipip="../rootfile/e45/phsp/pi_plus/pipip/";
  TString pi_plus_pipin="../rootfile/e45/phsp/pi_plus/pipin/";
  TString pi_minus_pipip="../rootfile/e45/phsp/pi_minus/pipip/";
  TString pi_minus_pipin="../rootfile/e45/phsp/pi_minus/pipin/";
  TString filename[8]={"p635_phsp.root","p835_phsp.root","p1035_phsp.root","p1235_phsp.root","p1435_phsp.root","p1635_phsp.root","p1835_phsp.root","p2000_phsp.root" };
  double momentum[8]={0.635,0.835,1.035,1.235,1.435,1.635,1.835,2.000};

  TFile *Delta_pipip[8];
  TFile *Delta_pipin[8];
  TFile *N_pipip[8];
  TFile *N_pipin[8];
  for(int i=0;i<8;i++){
    //Delta_pipip[i] = new TFile(Form("%s%s",pi_plus_pipip.Data(),filename[i].Data()),"READ");
    //Delta_pipin[i] = new TFile(Form("%s%s",pi_plus_pipin.Data(),filename[i].Data()),"READ");
    N_pipip[i] = new TFile(Form("%s%s",pi_minus_pipip.Data(),filename[i].Data()),"READ");
    N_pipin[i] = new TFile(Form("%s%s",pi_minus_pipin.Data(),filename[i].Data()),"READ");
  }

  double accep_N_pipip_wo[8];
  double accep_N_pipin_wo[8];
  double accep_N_pipip_w[8];
  double accep_N_pipin_w[8];
  for(int i=0;i<8;i++){
    hodo_class.pipip(N_pipip[i],1.0,false,accep_N_pipip_wo[i]);
    hodo_class.pipip(N_pipip[i],1.0,true,accep_N_pipip_w[i]);

    hodo_class.pipin(N_pipin[i],1.0,false,accep_N_pipin_wo[i]);
    hodo_class.pipin(N_pipin[i],1.0,true,accep_N_pipin_w[i]);
  }

  /*
  double accep_Delta_pipip_wo[8];
  double accep_Delta_pipin_wo[8];
  double accep_Delta_pipip_w[8];
  double accep_Delta_pipin_w[8];
  for(int i=0;i<8;i++){
    hodo_class.pipip(Delta_pipip[i],1.0,false,accep_Delta_pipip_wo[i]);
    hodo_class.pipip(Delta_pipip[i],1.0,true,accep_Delta_pipip_w[i]);

    hodo_class.pipin(Delta_pipin[i],1.0,false,accep_Delta_pipin_wo[i]);
    hodo_class.pipin(Delta_pipin[i],1.0,true,accep_Delta_pipin_w[i]);
  }
  */

  TCanvas *c1 = new TCanvas("c1","c1",1200, 1200);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1 -> SetTitle("Acceptance; Beam Momentum (GeV/c) ; Acceptance (Arb.)");
  TGraphErrors *gr_N_pipip_wo = new TGraphErrors(8,momentum,accep_N_pipip_wo);
  gr_N_pipip_wo -> SetMarkerColor(2);
  gr_N_pipip_wo -> SetMarkerStyle(20);
  mg1 -> Add(gr_N_pipip_wo);
  TGraphErrors *gr_N_pipip_w = new TGraphErrors(8,momentum,accep_N_pipip_w);
  gr_N_pipip_w -> SetMarkerColor(3);
  gr_N_pipip_w -> SetMarkerStyle(21);
  mg1 -> Add(gr_N_pipip_w);
  TGraphErrors *gr_N_pipin_wo = new TGraphErrors(8,momentum,accep_N_pipin_wo);
  gr_N_pipin_wo -> SetMarkerColor(4);
  gr_N_pipin_wo -> SetMarkerStyle(22);
  mg1 -> Add(gr_N_pipin_wo);
  TGraphErrors *gr_N_pipin_w = new TGraphErrors(8,momentum,accep_N_pipin_w);
  gr_N_pipin_w -> SetMarkerColor(5);
  gr_N_pipin_w -> SetMarkerStyle(23);
  mg1 -> Add(gr_N_pipin_w);

  c1 -> cd();
  mg1 -> Draw("AP");
  TLegend *legend1 = new TLegend(0.6,0.6,0.9,0.9);
  gStyle -> SetLegendFont(42);
  gStyle -> SetLegendTextSize(0.03);

  legend1 -> AddEntry(gr_N_pipip_wo,"","ep");
  legend1 -> AddEntry((TObject*)0,"MPV value corresconding to","");
  legend1 -> Draw("same");

  /*
    TGraphErrors *gr_Delta_pipip_wo = new TGraphErrors();
    TGraphErrors *gr_Delta_pipip_w = new TGraphErrors();
    TGraphErrors *gr_Delta_pipin_wo = new TGraphErrors();
    TGraphErrors *gr_Delta_pipin_w = new TGraphErrors();
  */
}//end
