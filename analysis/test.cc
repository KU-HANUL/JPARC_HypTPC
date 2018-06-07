#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <TPad.h>

using namespace std;

void test()
{
  bool pdf=false;
  //bool pdf=true;
  //bool printout=false;
  bool printout=true;

  gStyle->SetOptStat(0);
  //data load
  std::cout<<"data chain load"<<std::endl;
  std::cout<<std::endl;

  //TFile *file = new TFile("../rootfile/e45/p1035_pipip.root","READ");
  TFile *file = new TFile("../rootfile/e45/phsp/pi_plus/pipip/p1035_phsp.root","READ");
  TTree *tree = (TTree*)file->Get("tree");
  if(tree == 0) std::cout << "data open error : " << std::endl;

  int nhTof;
  int event;
  int nEvt;
  int tofpid[100];
  int toftrid[100];
  int tofseg[100];
  double toftime[100];
  double tofedep[100];
  double tofmomx[100];
  double tofmomy[100];
  double tofmomz[100];
  double tofpath[100];
  /*
    int tofparentpid1[nhTof];
    int tofparentpid2[nhTof];
    int tofparentpid3[nhTof];
    double tofposx[nhTof];
    double tofposy[nhTof];
    double tofposz[nhTof];
  */
  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("toftime",toftime);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("tofmomx",tofmomx);
  tree->SetBranchAddress("tofmomy",tofmomy);
  tree->SetBranchAddress("tofmomz",tofmomz);
  tree->SetBranchAddress("tofpath",tofpath);
  /*
    tree->SetBranchAddress("tofposx",tofposx);
    tree->SetBranchAddress("tofposy",tofposy);
    tree->SetBranchAddress("tofposz",tofposz);
    tree->SetBranchAddress("tofparentpid1",tofparentpid1);
    tree->SetBranchAddress("tofparentpid2",tofparentpid2);
    tree->SetBranchAddress("tofparentpid3",tofparentpid3);
  */

  TCanvas *can_tof = new TCanvas("can_tof","",1200,1200);
  TCanvas *can_ran_tof[6];
  for(int i=0;i<6;i++){
    can_ran_tof[i]=new TCanvas(Form("can_ran_tof%d",i),Form("can_ran_tof%d",i),1200,1200);
  }
  int resol[6]={150,160,170,180,190,200}; //ps
  double det_resol[6]={0.150,0.160,0.170,0.180,0.190,0.200};

  TH1D *hist_tof_p=new TH1D("hist_tof_p","Time of flight",1000,0,5);
  TH1D *hist_tof_pi=new TH1D("hist_tof_pi","Time of flight",1000,0,5);

  TH1D *hist_ran_tof_p[6];
  TH1D *hist_ran_tof_pi[6];
  for(int i=0;i<6;i++){
    hist_ran_tof_p[i]=new TH1D(Form("hist_ran_tof_p%d",i),Form("Time of flight /w #sigma = %d ps",resol[i]),1000,0,5);
    hist_ran_tof_pi[i]=new TH1D(Form("hist_ran_tof_pi%d",i),Form("Time of flight /w #sigma = %d ps",resol[i]),1000,0,5);
  }
  double c=0.299792458;
  int total_count,count_p,count_pi,count_mu,count_e,count_gamma; //flag for counting particles
  double tof_p,tof_pi; // time of flight of p,pi w/o detector resolution
  double tof_ran_p,tof_ran_pi; // time of flight of p,pi w/ detector resolution
  TRandom *eventgen = new TRandom();
  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    total_count=0;
    count_p=0;
    count_pi=0;
    count_e=0;
    count_mu=0;
    count_gamma=0;
    for(int j=0;j<nhTof;j++){
      if(TMath::Abs(tofpid[j])==211){ //for selecting events w/ pi
	total_count++;
	count_pi++;
	tof_pi=toftime[j];
      }
      if(tofpid[j]==2212){ //for selecting events w/ p
	total_count++;
	count_p++;
	tof_p=toftime[j];
      }
      if(tofpid[j]==11)	count_e++;
      if(tofpid[j]==13)	count_mu++;
      if(tofpid[j]==22)	count_gamma++;
    }
    if(count_pi>1){
      std::cout<<""<<std::endl;
      for(int j=0;j<nhTof;j++){
	if(TMath::Abs(tofpid[j])==211){
	  std::cout<<"seg : "<<tofseg[j]<<" id : "<<toftrid[j]<<std::endl;
	}
      }
    }
  }
}//end
