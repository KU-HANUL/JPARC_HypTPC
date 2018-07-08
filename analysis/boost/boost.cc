#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <fstream>
#include <math.h>
#include <iostream>

void boost(){

  //gStyle->SetOptStat(0);

  //TFile *file = new TFile("~/Desktop/hanulgit/JPARC_HypTPC/rootfile/e45/phsp/pi_plus/pipip/p2000_phsp.root","READ");
  TFile *file = new TFile("~/Desktop/p2000_phsp.root","READ");
  TTree *tree = (TTree*)file->Get("tree");

  int event;
  int nEvt;
  int evtpid[100];

  double evtpx[100];
  double evtpy[100];
  double evtpz[100];
  double evtvx[100];
  double evtvy[100];
  double evtvz[100];

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nEvt",&nEvt);

  tree->SetBranchAddress("evtpid",evtpid);
  tree->SetBranchAddress("evtpx",evtpx);
  tree->SetBranchAddress("evtpy",evtpy);
  tree->SetBranchAddress("evtpz",evtpz);
  tree->SetBranchAddress("evtvx",evtvx);
  tree->SetBranchAddress("evtvy",evtvy);
  tree->SetBranchAddress("evtvz",evtvz);

  double c = 0.299792458;
  double m_p = 0.938272;
  double m_pi = 0.139570;
  double m_pi0 = 0.134977;

  double p_beam = 2.0;
  double e_beam = sqrt(p_beam*p_beam + m_pi*m_pi);

  TVector3 beta_CM(0,0,-p_beam/(e_beam+m_p));

  double pT_pi, pL_pi, pT_p, pL_p, pT_pi0, pL_pi0;
  double pT_lab, pL_lab, pT_CM, pL_CM;
  double px_pi, py_pi, pz_pi, px_pi0, py_pi0, pz_pi0, px_p, py_p, pz_p;
  double E_pi, E_p, E_pi0;

  TH1D* hist_labpT_pi = new TH1D("hist_labpT_pi","Lab p_{T} #pi^{+}", 200, -2.2, 2.2);
  TH1D* hist_labpL_pi = new TH1D("hist_labpL_pi","Lab p_{L} #pi^{+}", 200, -2.2, 2.2);
  TH1D* hist_labpT_p = new TH1D("hist_labpT_p","Lab p_{T} p", 200, -2.2, 2.2);
  TH1D* hist_labpL_p = new TH1D("hist_labpL_p","Lab p_{L} p", 200, -2.2, 2.2);
  TH1D* hist_labpT_pi0 = new TH1D("hist_labpT_pi0","Lab p_{T} #pi^{0}", 200, -2.2, 2.2);
  TH1D* hist_labpL_pi0 = new TH1D("hist_labpL_pi0","Lab p_{L} #pi^{0}", 200, -2.2, 2.2);

  TH1D* hist_labpT = new TH1D("hist_labpT","Lab P_{T}", 200, -2.2, 2.2);
  TH1D* hist_labpL = new TH1D("hist_labpL","Lab P_{L}", 200, -2.2, 2.2);

  TH1D* hist_CMpT_pi = new TH1D("hist_CMpT_pi","C.M. p_{T} #pi^{+}", 200, -2.2, 2.2);
  TH1D* hist_CMpL_pi = new TH1D("hist_CMpL_pi","C.M. p_{L} #pi^{+}", 200, -2.2, 2.2);
  TH1D* hist_CMpT_p = new TH1D("hist_CMpT_p","C.M. p_{T} p", 200, -2.2, 2.2);
  TH1D* hist_CMpL_p = new TH1D("hist_CMpL_p","C.M. p_{L} p", 200, -2.2, 2.2);
  TH1D* hist_CMpT_pi0 = new TH1D("hist_CMpT_pi0","C.M. p_{T} #pi^{0}", 200, -2.2, 2.2);
  TH1D* hist_CMpL_pi0 = new TH1D("hist_CMpL_pi0","C.M. p_{L} #pi^{0}", 200, -2.2, 2.2);

  TH1D* hist_CM_E = new TH1D("hist_CM_E","C.M. Energy", 200, 0, 2.5);
  TH1D* hist_CMpT = new TH1D("hist_CMpT","C.M. P_{T}", 200, -2.2, 2.2);
  TH1D* hist_CMpL = new TH1D("hist_CMpL","C.M. P_{L}", 200, -2.2, 2.2);

  TCanvas *can_labpT = new TCanvas("can_labpT","",1200,800);
  can_labpT -> Divide(2,2);
  TCanvas *can_labpL = new TCanvas("can_labpL","",1200,800);
  can_labpL -> Divide(2,2);
  TCanvas *can_CMpT = new TCanvas("can_CMpT","",1200,800);
  can_CMpT -> Divide(2,2);
  TCanvas *can_CMpL = new TCanvas("can_CMpL","",1200,800);
  can_CMpL -> Divide(2,2);
  TCanvas *can_CM_E = new TCanvas("can_CM_E","",1200,800);

  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    for(int j=0;j<nEvt;j++){
      if(TMath::Abs(evtpid[j])==211){ //for selecting events w/ pi+
	pT_pi=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_pi=TMath::Sqrt(evtpz[j]*evtpz[j]);
	E_pi = sqrt(pT_pi*pT_pi+pL_pi*pL_pi + m_pi*m_pi);

	px_pi=evtpx[j];
	py_pi=evtpy[j];
	pz_pi=evtpz[j];

	hist_labpT_pi -> Fill(pT_pi);
	hist_labpL_pi -> Fill(pL_pi);
      }
      else if(evtpid[j]==111){ //for selecting events w/ pi0
	pT_pi0=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_pi0=TMath::Sqrt(evtpz[j]*evtpz[j]);
	E_pi0 = sqrt(pT_pi0*pT_pi0+pL_pi0*pL_pi0 + m_pi0*m_pi0);

	px_pi0=evtpx[j];
	py_pi0=evtpy[j];
	pz_pi0=evtpz[j];

	hist_labpT_pi0 -> Fill(pT_pi0);
	hist_labpL_pi0 -> Fill(pL_pi0);
      }
      else if(evtpid[j]==2212){ //for selecting events w/ p
    	pT_p=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_p=TMath::Sqrt(evtpz[j]*evtpz[j]);
	E_p = sqrt(pT_p*pT_p+pL_p*pL_p + m_p*m_p);

	px_p=evtpx[j];
	py_p=evtpy[j];
	pz_p=evtpz[j];

	hist_labpT_p -> Fill(pT_p);
	hist_labpL_p -> Fill(pL_p);
      }
    }
    pT_lab = TMath::Sqrt((px_pi+px_pi0+px_p)*(px_pi+px_pi0+px_p)+(py_pi+py_pi0+py_p)*(py_pi+py_pi0+py_p));
    pL_lab = pz_pi+pz_pi0+pz_p;

    hist_labpL -> Fill(pL_lab);
    hist_labpT -> Fill(pT_lab);

    TLorentzVector lv_pi(px_pi,py_pi,pz_pi,E_pi);
    lv_pi.Boost(beta_CM);
    hist_CMpT_pi -> Fill(lv_pi.Pt());
    hist_CMpL_pi -> Fill(lv_pi.Pz());

    TLorentzVector lv_pi0(px_pi0,py_pi0,pz_pi0,E_pi0);
    lv_pi0.Boost(beta_CM);
    hist_CMpT_pi0 -> Fill(lv_pi0.Pt());
    hist_CMpL_pi0 -> Fill(lv_pi0.Pz());

    TLorentzVector lv_p(px_p,py_p,pz_p,E_p);
    lv_p.Boost(beta_CM);
    hist_CMpT_p -> Fill(lv_p.Pt());
    hist_CMpL_p -> Fill(lv_p.Pz());

    TLorentzVector lv_CM;
    lv_CM=lv_p+lv_pi+lv_pi0;
    hist_CM_E -> Fill(lv_CM.E());
    hist_CMpT -> Fill(lv_CM.Pt());
    hist_CMpL -> Fill(lv_CM.Pz());
  }

  can_labpT -> cd(1);
  hist_labpT_pi -> Draw();
  can_labpT -> cd(2);
  hist_labpT_pi0 -> Draw();
  can_labpT -> cd(3);
  hist_labpT_p -> Draw();
  can_labpT -> cd(4);
  hist_labpT -> Draw();

  can_labpL -> cd(1);
  hist_labpL_pi -> Draw();
  can_labpL -> cd(2);
  hist_labpL_pi0 -> Draw();
  can_labpL -> cd(3);
  hist_labpL_p -> Draw();
  can_labpL -> cd(4);
  hist_labpL -> Draw();

  can_CMpL -> cd(1);
  hist_CMpL_pi -> Draw();
  can_CMpL -> cd(2);
  hist_CMpL_pi0 -> Draw();
  can_CMpL -> cd(3);
  hist_CMpL_p -> Draw();
  can_CMpL -> cd(4);
  hist_CMpL -> Draw();

  can_CMpT -> cd(1);
  hist_CMpT_pi -> Draw();
  can_CMpT -> cd(2);
  hist_CMpT_pi0 -> Draw();
  can_CMpT -> cd(3);
  hist_CMpT_p -> Draw();
  can_CMpT -> cd(4);
  hist_CMpT -> Draw();

  can_CM_E -> cd();
  hist_CM_E -> Draw();
}
