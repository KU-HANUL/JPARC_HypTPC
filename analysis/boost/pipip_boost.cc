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

double GetCos(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Px()*lv2.Px() + lv1.Py()*lv2.Py() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / (lv1.P() * lv2.P());
  return cos;
}

double GetCos_xz(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Px()*lv2.Px() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / TMath::Sqrt((lv1.Px()*lv1.Px()+lv1.Pz()*lv1.Pz()) * (lv2.Px()*lv2.Px()+lv2.Pz()*lv2.Pz()));
  return cos;
}

double GetCos_yz(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Py()*lv2.Py() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / TMath::Sqrt((lv1.Py()*lv1.Py()+lv1.Pz()*lv1.Pz()) * (lv2.Py()*lv2.Py()+lv2.Pz()*lv2.Pz()));
  return cos;
}

double GetAng(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(GetCos(lv1,lv2)))/TMath::Pi();
  return ang;
}

double GetAng_xz(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(GetCos_xz(lv1,lv2)))/TMath::Pi();
  return ang;
}

double GetAng_yz(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(GetCos_yz(lv1,lv2)))/TMath::Pi();
  return ang;
}

void pipip_boost(){

  //gStyle->SetOptStat(0);

  TFile *file = new TFile("~/Desktop/hanul_git/JPARC_HypTPC/rootfile/e45/phsp/pi_plus/pipip/p2000_phsp.root","READ");
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
  double theta_CM_p, cos_CM_p, theta_CM_pi, cos_CM_pi, theta_lab_pip, cos_lab_pip, theta_lab_pip_xz,theta_lab_pip_yz;

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
  TH1D* hist_lab_angle = new TH1D("hist_lab_angle","Lab frame, P, #pi^{+} Angle Distribution;#theta_{P, #pi^{+}};counts", 100, 0, 180);
  TH1D* hist_lab_angle_xz = new TH1D("hist_lab_angle_xz","Lab frame, XZ plane projection, P, #pi^{+} Angle Distribution;#theta_{P, #pi^{+}};counts", 100, 0, 180);
  TH1D* hist_lab_angle_xz_selected = new TH1D("hist_lab_angle_xz_selected","Lab frame, XZ plane projection, P, #pi^{+} Angle Distribution;#theta_{P, #pi^{+}};counts", 100, 0, 180);
  TH1D* hist_lab_angle_yz = new TH1D("hist_lab_angle_yz","Lab frame, YZ plane projection, P, #pi^{+} Angle Distribution;#theta_{P, #pi^{+}};counts", 100, 0, 180);
  TH1D* hist_lab_cos = new TH1D("hist_lab_cos","Lab frame, P, #pi^{+} Angle Distribution;Cos(#theta_{P, #pi^{+}});counts", 100, -1, 1);
  TH1D* hist_CM_angle_p = new TH1D("hist_CM_angle_p","Proton Angle Distribution;#theta_{C.M. Proton};counts", 100, 0, 180);
  TH1D* hist_CM_cos_p = new TH1D("hist_CM_cos_p","Proton Angle Distribution;Cos(#theta_{C.M. Proton});counts", 100, -1, 1);
  TH1D* hist_CM_angle_pi = new TH1D("hist_CM_angle_pi","#pi^{+} Angle Distribution;#theta_{C.M. #pi^{+}};counts", 100, 0, 180);
  TH1D* hist_CM_cos_pi = new TH1D("hist_CM_cos_pi","#pi^{+} Angle Distribution;Cos(#theta_{C.M. #pi^{+}});counts", 100, -1, 1);

  TCanvas *can_labpT = new TCanvas("can_labpT","",1200,800);
  can_labpT -> Divide(2,2);
  TCanvas *can_labpL = new TCanvas("can_labpL","",1200,800);
  can_labpL -> Divide(2,2);
  TCanvas *can_CMpT = new TCanvas("can_CMpT","",1200,800);
  can_CMpT -> Divide(2,2);
  TCanvas *can_CMpL = new TCanvas("can_CMpL","",1200,800);
  can_CMpL -> Divide(2,2);
  TCanvas *can_CM_E = new TCanvas("can_CM_E","",1200,800);
  TCanvas *can_lab_angle = new TCanvas("can_lab_angle","",1200,800);
  can_lab_angle -> Divide(2,2);
  TCanvas *can_CM_angle = new TCanvas("can_CM_angle","",1200,800);
  can_CM_angle -> Divide(2,2);

  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    for(int j=0;j<nEvt;j++){
      if(TMath::Abs(evtpid[j])==211){ //for selecting events w/ pi+
	pT_pi=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_pi=evtpz[j];
	E_pi = sqrt(pT_pi*pT_pi+pL_pi*pL_pi + m_pi*m_pi);

	px_pi=evtpx[j];
	py_pi=evtpy[j];
	pz_pi=evtpz[j];

	hist_labpT_pi -> Fill(pT_pi);
	hist_labpL_pi -> Fill(pL_pi);
      }
      else if(evtpid[j]==111){ //for selecting events w/ pi0
	pT_pi0=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_pi0=evtpz[j];
	E_pi0 = sqrt(pT_pi0*pT_pi0+pL_pi0*pL_pi0 + m_pi0*m_pi0);

	px_pi0=evtpx[j];
	py_pi0=evtpy[j];
	pz_pi0=evtpz[j];

	hist_labpT_pi0 -> Fill(pT_pi0);
	hist_labpL_pi0 -> Fill(pL_pi0);
      }
      else if(evtpid[j]==2212){ //for selecting events w/ p
    	pT_p=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_p=evtpz[j];
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

    //lab frame
    TLorentzVector lv_pi(px_pi,py_pi,pz_pi,E_pi);
    TLorentzVector lv_pi0(px_pi0,py_pi0,pz_pi0,E_pi0);
    TLorentzVector lv_p(px_p,py_p,pz_p,E_p);
    TLorentzVector lv_beam(0,0,p_beam,e_beam);

    theta_lab_pip = GetAng(lv_p,lv_pi);
    cos_lab_pip = GetCos(lv_p,lv_pi);
    theta_lab_pip_xz = GetAng_xz(lv_p,lv_pi);
    theta_lab_pip_yz = GetAng_yz(lv_p,lv_pi);

    hist_lab_angle -> Fill(theta_lab_pip);
    hist_lab_cos -> Fill(cos_lab_pip);
    hist_lab_angle_xz -> Fill(theta_lab_pip_xz);
    if(theta_lab_pip_xz<10){
      hist_lab_angle_xz_selected -> Fill(theta_lab_pip_xz);
      hist_lab_angle_yz -> Fill(theta_lab_pip_yz);
    }

    //C.M. frame
    lv_pi.Boost(beta_CM);
    lv_pi0.Boost(beta_CM);
    lv_p.Boost(beta_CM);

    hist_CMpT_pi -> Fill(lv_pi.Pt());
    hist_CMpL_pi -> Fill(lv_pi.Pz());
    hist_CMpT_pi0 -> Fill(lv_pi0.Pt());
    hist_CMpL_pi0 -> Fill(lv_pi0.Pz());
    hist_CMpT_p -> Fill(lv_p.Pt());
    hist_CMpL_p -> Fill(lv_p.Pz());

    theta_CM_pi = GetAng(lv_beam,lv_pi);
    cos_CM_pi = GetCos(lv_beam,lv_pi);
    theta_CM_p = GetAng(lv_beam,lv_p);
    cos_CM_p = GetCos(lv_beam,lv_p);

    hist_CM_angle_pi -> Fill(theta_CM_pi);
    hist_CM_cos_pi -> Fill(cos_CM_pi);
    hist_CM_angle_p -> Fill(theta_CM_p);
    hist_CM_cos_p -> Fill(cos_CM_p);

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

  can_lab_angle -> cd(1);
  hist_lab_angle -> Draw();
  can_lab_angle -> cd(2);
  hist_lab_cos -> Draw();
  can_lab_angle -> cd(3);
  hist_lab_angle_xz -> Draw();
  hist_lab_angle_xz_selected -> Draw("same");
  hist_lab_angle_xz_selected -> SetLineColor(2);
  can_lab_angle -> cd(4);
  hist_lab_angle_yz -> Draw();

  can_CM_angle -> cd(1);
  hist_CM_angle_pi -> Draw();
  can_CM_angle -> cd(2);
  hist_CM_cos_pi -> Draw();
  hist_CM_cos_pi -> GetYaxis() -> SetRangeUser(0,150);
  can_CM_angle -> cd(3);
  hist_CM_angle_p -> Draw();
  can_CM_angle -> cd(4);
  hist_CM_cos_p -> Draw();
  hist_CM_cos_p -> GetYaxis() -> SetRangeUser(0,150);

}
