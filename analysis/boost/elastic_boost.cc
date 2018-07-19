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

void elastic_boost(){

  //gStyle->SetOptStat(0);

  //TFile *file = new TFile("~/Desktop/hanul_git/JPARC_HypTPC/rootfile/e45/elastic/pi_plus/pipip/p2000_phsp.root","READ");
  //TFile *file = new TFile("~/Desktop/hanul_git/JPARC_HypTPC/rootfile/e45/elastic/pi_minus/pipip/p2000_phsp.root","READ");
  //TFile *file = new TFile("~/Desktop/test.root","READ");
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

  TF1 *legen[10];
  double coefficient[10];
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

  TF1 *func = new TF1("func","legen[0]*[0]+legen[1]*[1]+legen[2]*[2]+legen[3]*[3]+legen[4]*[4]+legen[5]*[5]+legen[6]*[6]+legen[7]*[7]+legen[8]*[8]+legen[9]*[9]",-1.0,1.0);

  double c = 0.299792458;
  double m_p = 0.938272;
  double m_pi = 0.139570;

  double p_beam = 2.0;
  double e_beam = sqrt(p_beam*p_beam + m_pi*m_pi);

  TVector3 beta_CM(0,0,-p_beam/(e_beam+m_p));

  double pT_pi, pL_pi, pT_p, pL_p;
  double pT_lab, pL_lab, pT_CM, pL_CM;
  double px_pi, py_pi, pz_pi, px_p, py_p, pz_p;
  double E_pi, E_p;
  double theta_CM, cos_CM, theta_lab_pip, cos_lab_pip, theta_lab_pip_xz, theta_lab_pip_yz,theta_lab_p,theta_lab_pi;

  TGraph *gr = new TGraph();

  TH1D* hist_labpT_pi = new TH1D("hist_labpT_pi","Lab p_{T} #pi^{+}", 500, -2.5, 2.5);
  TH1D* hist_labpL_pi = new TH1D("hist_labpL_pi","Lab p_{L} #pi^{+}", 500, -2.5, 2.5);
  TH1D* hist_labpT_p = new TH1D("hist_labpT_p","Lab p_{T} p", 500, -2.5, 2.5);
  TH1D* hist_labpL_p = new TH1D("hist_labpL_p","Lab p_{L} p", 500, -2.5, 2.5);

  TH1D* hist_labpT = new TH1D("hist_labpT","Lab P_{T}", 500, -2.5, 2.5);
  TH1D* hist_labpL = new TH1D("hist_labpL","Lab P_{L}", 500, -2.5, 2.5);

  TH1D* hist_CMpT_pi = new TH1D("hist_CMpT_pi","C.M. p_{T} #pi^{+}", 500, -2.5, 2.5);
  TH1D* hist_CMpL_pi = new TH1D("hist_CMpL_pi","C.M. p_{L} #pi^{+}", 500, -2.5, 2.5);
  TH1D* hist_CMpT_p = new TH1D("hist_CMpT_p","C.M. p_{T} p", 500, -2.5, 2.5);
  TH1D* hist_CMpL_p = new TH1D("hist_CMpL_p","C.M. p_{L} p", 500, -2.5, 2.5);

  TH1D* hist_CM_E = new TH1D("hist_CM_E","C.M. Energy", 500, 0, 2.5);
  TH1D* hist_CMpT = new TH1D("hist_CMpT","C.M. P_{T}", 500, -2.5, 2.5);
  TH1D* hist_CMpL = new TH1D("hist_CMpL","C.M. P_{L}", 500, -2.5, 2.5);

  TH1D* hist_lab_angle = new TH1D("hist_lab_angle","Lab frame, P, #pi^{+} Angle Distribution;#theta_{Lab, P, #pi^{+}};counts", 500, 0, 180);
  TH1D* hist_lab_angle_xz = new TH1D("hist_lab_angle_xz","Lab frame, XZ plane projection, P, #pi^{+} Angle Distribution;#theta_{Lab, P, #pi^{+}};counts", 500, 0, 180);
  TH1D* hist_lab_angle_xz_selected = new TH1D("hist_lab_angle_xz_selected","Lab frame, XZ plane projection, P, #pi^{+} Angle Distribution;#theta_{Lab, P, #pi^{+}};counts", 500, 0, 180);
  TH1D* hist_lab_angle_yz = new TH1D("hist_lab_angle_yz","Lab frame, YZ plane projection, P, #pi^{+} Angle Distribution;#theta_{Lab, P, #pi^{+}};counts", 500, 0, 180);
  TH1D* hist_lab_cos = new TH1D("hist_lab_cos","Lab frame, P, #pi^{+} Angle Distribution;Cos(#theta_{P, #pi^{+}});counts", 500, -1, 1);
  TH1D* hist_CM_angle = new TH1D("hist_CM_angle","Scattering Angle Distribution;#theta_{C.M.};counts", 500, 0, 180);
  TH1D* hist_CM_cos = new TH1D("hist_CM_cos","Scattering Angle Distribution;Cos(#theta_{C.M.});counts", 500, -1, 1);

  TH2D* hist_dummy = new TH2D("hist_dummy","Scatter plot of scattering angles of p, #pi;Scattering angle, #theta_{p};Scattering angle, #theta_{pi}", 100,0,180,100,0,180);

  TCanvas *can_scatterplot = new TCanvas("can_scatterplot","",800,800);
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

  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    for(int j=0;j<nEvt;j++){
      if(TMath::Abs(evtpid[j])==211){ //for selecting events w/ pi+
	pT_pi=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_pi=evtpz[j];
	E_pi =TMath::Sqrt(pT_pi*pT_pi+pL_pi*pL_pi + m_pi*m_pi);

	px_pi=evtpx[j];
	py_pi=evtpy[j];
	pz_pi=evtpz[j];

	hist_labpT_pi -> Fill(pT_pi);
	hist_labpL_pi -> Fill(pL_pi);
      }
      else if(evtpid[j]==2212){ //for selecting events w/ p
    	pT_p=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]);
	pL_p=evtpz[j];
	E_p = TMath::Sqrt(pT_p*pT_p+pL_p*pL_p + m_p*m_p);

	px_p=evtpx[j];
	py_p=evtpy[j];
	pz_p=evtpz[j];

	hist_labpT_p -> Fill(pT_p);
	hist_labpL_p -> Fill(pL_p);
      }
    }
    pT_lab = TMath::Sqrt((px_pi+px_p)*(px_pi+px_p)+(py_pi+py_p)*(py_pi+py_p));
    pL_lab = pz_pi+pz_p;

    hist_labpL -> Fill(pL_lab);
    hist_labpT -> Fill(pT_lab);

    //lab frame
    TLorentzVector lv_pi(px_pi,py_pi,pz_pi,E_pi); //p3
    TLorentzVector lv_p(px_p,py_p,pz_p,E_p); //p4
    //TLorentzVector lv_target(0,0,0,m_p); //p2
    TLorentzVector lv_beam(0,0,p_beam,e_beam); //p1

    theta_lab_pi = GetAng(lv_beam,lv_pi);
    theta_lab_p = GetAng(lv_beam,lv_p);
    gr -> SetPoint(i,theta_lab_pi,theta_lab_p);

    theta_lab_pip = GetAng(lv_p,lv_pi);
    cos_lab_pip = GetCos(lv_p,lv_pi);
    theta_lab_pip_xz = GetAng_xz(lv_p,lv_pi);
    theta_lab_pip_yz = GetAng_yz(lv_p,lv_pi);

    hist_lab_angle -> Fill(theta_lab_pip);
    hist_lab_cos -> Fill(cos_lab_pip);
    hist_lab_angle_xz -> Fill(theta_lab_pip_xz);
    if(theta_lab_pip_xz<20){
      hist_lab_angle_xz_selected -> Fill(theta_lab_pip_xz);
      hist_lab_angle_yz -> Fill(theta_lab_pip_yz);
    }

    //C.M. frame
    lv_pi.Boost(beta_CM); //p3'
    lv_p.Boost(beta_CM); //p4'
    //lv_target.Boost(beta_CM); //p2'
    lv_beam.Boost(beta_CM); //p1'

    hist_CMpT_pi -> Fill(lv_pi.Pt());
    hist_CMpL_pi -> Fill(lv_pi.Pz());
    hist_CMpT_p -> Fill(lv_p.Pt());
    hist_CMpL_p -> Fill(lv_p.Pz());

    theta_CM = GetAng(lv_beam,lv_pi);
    cos_CM = GetCos(lv_beam,lv_pi);
    hist_CM_angle -> Fill(theta_CM);
    hist_CM_cos -> Fill(cos_CM);

    TLorentzVector lv_CM;
    lv_CM=lv_p+lv_pi;
    hist_CM_E -> Fill(lv_CM.E());
    hist_CMpT -> Fill(lv_CM.Pt());
    hist_CMpL -> Fill(lv_CM.Pz());
  }

  can_scatterplot -> cd();
  //hist_dummy -> Draw();
  gr -> Draw("AP");
  gr -> SetTitle("Scatter plot for scattering angles of p, #pi");
  gr -> GetHistogram() -> GetXaxis()-> SetTitle("Scattering angle, #theta_{#pi}");
  gr -> GetHistogram() -> GetXaxis()-> SetRangeUser(0,180);
  gr -> GetHistogram() -> GetYaxis()-> SetTitle("Scattering angle, #theta_{p}");
  gr -> GetHistogram() -> GetYaxis()-> SetRangeUser(0,180);

  can_labpT -> cd(1);
  hist_labpT_pi -> Draw();
  can_labpT -> cd(2);
  hist_labpT_p -> Draw();
  can_labpT -> cd(3);
  hist_labpT -> Draw();

  can_labpL -> cd(1);
  hist_labpL_pi -> Draw();
  can_labpL -> cd(2);
  hist_labpL_p -> Draw();
  can_labpL -> cd(3);
  hist_labpL -> Draw();

  can_CMpL -> cd(1);
  hist_CMpL_pi -> Draw();
  can_CMpL -> cd(2);
  hist_CMpL_p -> Draw();
  can_CMpL -> cd(3);
  hist_CMpL -> Draw();

  can_CMpT -> cd(1);
  hist_CMpT_pi -> Draw();
  can_CMpT -> cd(2);
  hist_CMpT_p -> Draw();
  can_CMpT -> cd(3);
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

  can_CM_angle -> cd();
  hist_CM_cos -> Draw();
  //hist_CM_cos -> GetYaxis() -> SetRangeUser(0,1200);
  hist_CM_cos -> Fit(func);

  func -> GetParameters(coefficient);
  for(int i=0;i<10;i++){
    std::cout<<"coefficient "<<i<<" th : "<<0.676241*coefficient[i]/coefficient[0]<<std::endl;
  }
}
