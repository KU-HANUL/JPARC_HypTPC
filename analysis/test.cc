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
#include "hodo_3body.h"

using namespace std;

void test()
{
  //gStyle -> SetOptFit(1000);

  bool pdf_saving = true;
  //bool pdf_saving = false;

  hodo_3body hodo_class;

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
    Delta_pipip[i] = new TFile(Form("%s%s",pi_plus_pipip.Data(),filename[i].Data()),"READ");
    Delta_pipin[i] = new TFile(Form("%s%s",pi_plus_pipin.Data(),filename[i].Data()),"READ");
    N_pipip[i] = new TFile(Form("%s%s",pi_minus_pipip.Data(),filename[i].Data()),"READ");
    N_pipin[i] = new TFile(Form("%s%s",pi_minus_pipin.Data(),filename[i].Data()),"READ");
  }

  double accep_N_pipip_wo[8];
  double accep_N_pipin_wo[8];
  double accep_N_pipip_w[8];
  double accep_N_pipin_w[8];
  for(int i=7;i<8;i++){
    hodo_class.pipip(N_pipip[i],1.0,false,accep_N_pipip_wo[i],true);
    //hodo_class.pipip(N_pipip[i],1.0,true,accep_N_pipip_w[i],true);

    //hodo_class.pipin(N_pipin[i],1.0,false,accep_N_pipin_wo[i],true);
    //hodo_class.pipin(N_pipin[i],1.0,true,accep_N_pipin_w[i],true);
  }


  double accep_Delta_pipip_wo[8];
  double accep_Delta_pipin_wo[8];
  double accep_Delta_pipip_w[8];
  double accep_Delta_pipin_w[8];
  for(int i=7;i<8;i++){
    //hodo_class.pipip(Delta_pipip[i],1.0,false,accep_Delta_pipip_wo[i],true);
    //hodo_class.pipip(Delta_pipip[i],1.0,true,accep_Delta_pipip_w[i],true);

    //hodo_class.pipin2(Delta_pipin[i],1.0,false,accep_Delta_pipin_wo[i],true);
    //hodo_class.pipin2(Delta_pipin[i],1.0,true,accep_Delta_pipin_w[i]);
  }


}//end
