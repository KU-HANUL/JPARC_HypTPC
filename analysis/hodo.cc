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

void hodo()
{
  //bool pdf=false;
  bool pdf=true;
  //bool printout=false;
  bool printout=true;

  //gStyle->SetOptStat(0);
  //data load
  std::cout<<"data chain load"<<std::endl;
  std::cout<<std::endl;

  TFile *file = new TFile("../rootfile/e45/p1035_pipip.root","READ");
  TTree *tree = (TTree*)file->Get("tree");
  if(tree == 0) std::cout << "data open error : " << std::endl;

  int nhTof;
  int event;
  int nEvt;
  int tofpid[100];
  double toftime[100];
  /*
  int tofseg[nhTof];
  double tofedep[nhTof];
  double tofposx[nhTof];
  double tofposy[nhTof];
  double tofposz[nhTof];
  double tofmomx[nhTof];
  double tofmomy[nhTof];
  double tofmomz[nhTof];
  int toftrid[nhTof];
  double tofpath[nhTof];
  int tofparentpid1[nhTof];
  int tofparentpid2[nhTof];
  int tofparentpid3[nhTof];
  */
  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftime",toftime);
  /*
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("tofposx",tofposx);
  tree->SetBranchAddress("tofposy",tofposy);
  tree->SetBranchAddress("tofposz",tofposz);

  tree->SetBranchAddress("tofmomx",tofmomx);
  tree->SetBranchAddress("tofmomy",tofmomy);
  tree->SetBranchAddress("tofmomz",tofmomz);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofpath",tofpath);

  tree->SetBranchAddress("tofparentpid1",tofparentpid1);
  tree->SetBranchAddress("tofparentpid2",tofparentpid2);
  tree->SetBranchAddress("tofparentpid3",tofparentpid3);
  */

  TCanvas *can_multi = new TCanvas("can_multi","",1200,1200);
  TCanvas *can_pid = new TCanvas("can_pid","",1200,1200);
  TCanvas *can_tof = new TCanvas("can_tof","",1200,1200);
  can_tof -> Divide(2,1);

  //TCanvas *can_ran_tof = new TCanvas("can_ran_tof","",1200,1200);

  TCanvas *can_ran_tof[6];
  for(int i=0;i<6;i++){
    can_ran_tof[i]=new TCanvas(Form("can_ran_tof%d",i),Form("can_ran_tof%d",i),1200,1200);
    can_ran_tof[i] -> Divide(2,1);
  }


  double det_resol[6]={0.150,0.160,0.170,0.180,0.190,0.200};
  TH1D *hist_pid =new TH1D("hist_pid","Hodoscope PID",4600,-2300,2300);
  TH1D *hist_multi =new TH1D("hist_multi","Hodoscope Multiplicity",10,0,10);
  TH1D *hist_multi_pi=new TH1D("hist_multi_pi","Hodoscope Multiplicity",10,0,10);
  TH1D *hist_tof_p=new TH1D("hist_tof_p","Time of flight",100,0,5);
  TH1D *hist_tof_pi=new TH1D("hist_tof_pi","Time of flight",100,0,5);
  TH1D *hist_deltatof=new TH1D("hist_deltatof","#Delta Time of flight",100,0,5);
  TH1D *hist_ran_p[6];
  TH1D *hist_ran_pi[6];
  TH1D *hist_ran_deltatof[6];
  for(int i=0;i<6;i++){
    hist_ran_p[i]=new TH1D(Form("hist_ran_p%d",i),Form("Time of flight /w #sigma = %f ps",det_resol[i]),100,0,5);
    hist_ran_pi[i]=new TH1D(Form("hist_ran_pi%d",i),Form("Time of flight /w #sigma = %f ps",det_resol[i]),100,0,5);
    hist_ran_deltatof[i]=new TH1D(Form("hist_ran_deltatof%d",i),Form("#Delta Time of flight /w #sigma = %f ps",det_resol[i]*1000),100,0,5);
  }


  int dummy_2hits=0,dummy_2hits_wPpi=0;
  int dummy_p,dummy_n;
  double tof_p,tof_pi;
  double tof_ran_p,tof_ran_pi;
  TRandom *eventgen = new TRandom();

  can_pid -> cd();
  tree -> Draw("tofpid>>hist_pid");
  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    //for(int i=0;i<10;i++){
    tree -> GetEntry(i);
    if(nhTof>1) dummy_2hits++;
    dummy_p=0;
    dummy_n=0;
    hist_multi -> Fill(nhTof);
    for(int j=0;j<nhTof;j++){
      if(tofpid[j]==-211){
	dummy_n++;
	tof_pi=toftime[j];
      }
      if(tofpid[j]==2212){
	dummy_p++;
	tof_p=toftime[j];
      }
    }
    if(dummy_p==1&&dummy_n==1){
      hist_multi_pi -> Fill(nhTof);
      dummy_2hits_wPpi++;
      hist_tof_p -> Fill(tof_p);
      hist_tof_pi -> Fill(tof_pi);
      hist_deltatof -> Fill(tof_p-tof_pi);

      //tof_ran_p=eventgen->Gaus(tof_p,det_resol[0]);
      //tof_ran_p=eventgen->Gaus(0,1);
      //cout << "tof_ran_p  : " << tof_ran_p << endl;
      for(int j=0;j<6;j++){
	tof_ran_p=eventgen->Gaus(tof_p,det_resol[j]);
	tof_ran_pi=eventgen->Gaus(tof_pi,det_resol[j]);

	hist_ran_p[j] -> Fill(tof_ran_p);
	hist_ran_pi[j]-> Fill(tof_ran_pi);
	hist_ran_deltatof[j]-> Fill(tof_ran_p-tof_ran_pi);
      }
    }
  }
  can_multi->cd();
  hist_multi->Draw();
  hist_multi_pi->Draw("same");
  hist_multi_pi->SetLineColor(kRed);
  hist_multi->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi->GetYaxis()->SetTitle("Counts");

  can_tof->cd(1);
  hist_tof_pi->Draw();
  hist_tof_p->Draw("same");
  hist_tof_pi->SetLineColor(kRed);
  hist_tof_pi->GetXaxis()->SetTitle("TOF(ns) ");
  hist_tof_pi->GetYaxis()->SetTitle("Counts");

  can_tof->cd(2);
  hist_deltatof->Draw();
  hist_deltatof->GetXaxis()->SetTitle("TOF(ns) ");
  hist_deltatof->GetYaxis()->SetTitle("Counts");

  //can_ran_tof->cd();
  //hist_ran_pi[0]->Draw();


  for(int i=0;i<6;i++){
    can_ran_tof[i]->cd(1);
    hist_ran_pi[i]->Draw();
    hist_ran_p[i]->Draw("same");
    hist_ran_pi[i]->SetLineColor(kRed);
    hist_ran_pi[i]->GetXaxis()->SetTitle("TOF(ns) ");
    hist_ran_pi[i]->GetYaxis()->SetTitle("Counts");

    can_ran_tof[i]->cd(2);
    hist_ran_deltatof[i]->Draw();
    hist_ran_deltatof[i]->GetXaxis()->SetTitle("TOF(ns) ");
    hist_ran_deltatof[i]->GetYaxis()->SetTitle("Counts");

  }

  if(pdf){
    can_pid ->SaveAs("pdf/hodo_pid");
    can_multi ->SaveAs("pdf/hodo_multiplicity");
  }
  if(printout){
    std::cout<<"Event number  : "<<nevent<<std::endl;
    std::cout<<"2 hits        : "<<dummy_2hits<<std::endl;
    std::cout<<"2 hits w/ P pi: "<<dummy_2hits_wPpi<<std::endl;
    std::cout<<"Percentage    : "<<(double) dummy_2hits/nevent *100<<"("<<(double) dummy_2hits_wPpi/nevent *100<<") %"<<std::endl;
  }
}//end
