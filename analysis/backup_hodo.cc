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

  TCanvas *can_multi = new TCanvas("can_multi","",1200,1200);
  TCanvas *can_pid = new TCanvas("can_pid","",1200,1200);
  TCanvas *can_tof = new TCanvas("can_tof","",1200,1200);
  TCanvas *can_ran_tof[6];
  for(int i=0;i<6;i++){
    can_ran_tof[i]=new TCanvas(Form("can_ran_tof%d",i),Form("can_ran_tof%d",i),1200,1200);
  }
  TCanvas *can_pbeta = new TCanvas("can_pbeta","",1200,1200);
  TCanvas *can_ran_pbeta[6];
  for(int i=0;i<6;i++){
    can_ran_pbeta[i]=new TCanvas(Form("can_ran_pbeta%d",i),Form("can_ran_pbeta%d",i),1200,1200);
  }
  TCanvas *can_ptof = new TCanvas("can_ptof","",1200,1200);
  TCanvas *can_ran_ptof[6];
  for(int i=0;i<6;i++){
    can_ran_ptof[i]=new TCanvas(Form("can_ran_ptof%d",i),Form("can_ran_ptof%d",i),1200,1200);
  }

  double ecut=1.0; //energy cut
  int resol[6]={150,160,170,180,190,200}; //ps
  double det_resol[6]={0.150,0.160,0.170,0.180,0.190,0.200};
  TH1D *hist_pid =new TH1D("hist_pid","Hodoscope PID",4600,-2300,2300);
  TH1D *hist_multi =new TH1D("hist_multi","Hodoscope Multiplicity",7,0,7);
  TH1D *hist_multi_pi=new TH1D("hist_multi_pi","Hodoscope Multiplicity",7,0,7);
  TH1D *hist_tof_p=new TH1D("hist_tof_p","Time of flight",100,0,5);
  TH1D *hist_tof_pi=new TH1D("hist_tof_pi","Time of flight",100,0,5);
  TH1D *hist_deltatof=new TH1D("hist_deltatof","#Delta Time of flight",100,0,5);
  TH2D *hist_pbeta=new TH2D("hist_pbeta","#beta vs p/q",500,-1.5,1.5,500,0.2,1.2);
  TH2D *hist_ptof=new TH2D("hist_ptof","TOF vs p/q",500,-1.5,1.5,500,0,20);

  TH1D *hist_ran_p[6];
  TH1D *hist_ran_pi[6];
  TH1D *hist_ran_deltatof[6];
  TH2D *hist_ran_pbeta[6];
  TH2D *hist_ran_ptof[6];
  for(int i=0;i<6;i++){
    hist_ran_p[i]=new TH1D(Form("hist_ran_p%d",i),Form("Time of flight /w #sigma = %d ps",resol[i]),100,0,5);
    hist_ran_pi[i]=new TH1D(Form("hist_ran_pi%d",i),Form("Time of flight /w #sigma = %d ps",resol[i]),100,0,5);
    hist_ran_deltatof[i]=new TH1D(Form("hist_ran_deltatof%d",i),Form("#Delta Time of flight /w #sigma = %d ps",resol[i]),100,0,5);
    hist_ran_pbeta[i]=new TH2D(Form("hist_ran_pbeta%d",i),Form("#beta vs p/q /w #sigma = %d ps",resol[i]),500,-1.5,1.5,500,0.2,1.2);
    hist_ran_ptof[i]=new TH2D(Form("hist_ran_ptof%d",i),Form("TOF vs p/q /w #sigma = %d ps",resol[i]),500,-1.5,1.5,500,0,20);
  }
  double c=0.299792458;
  int event_w_2hits=0,event_w_2hits_wPpi=0; //flag for counting events with multi hits
  int total_count,count_all,count_p,count_pi,count_mu,count_e,count_gamma; //flag for counting particles
  int seg_p,seg_pi;
  double tof_p,tof_pi; // time of flight of p,pi w/o detector resolution
  double beta_p,beta_pi; // beta of p,pi w/o detector resolution
  double path_p,path_pi; // path length of p,pi
  double tof_ran_p,tof_ran_pi; // time of flight of p,pi w/ detector resolution
  double beta_ran_p,beta_ran_pi; // beta of p,pi w/ detector resolutiona
  double p_p,p_pi; //p,pi momentum

  TRandom *eventgen = new TRandom();

  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    total_count=0;
    count_all=0;
    count_p=0;
    count_pi=0;
    count_e=0;
    count_mu=0;
    count_gamma=0;
    for(int j=0;j<nhTof;j++){
      hist_pid -> Fill(tofpid[j]); //PID
      if(TMath::Abs(tofpid[j])==211&&tofedep[j]>ecut){ //for selecting events w/ pi
	total_count++;
	count_pi++;
	tof_pi=toftime[j];
	path_pi=tofpath[j]/1000; //pi pathlength(mm -> m)
	beta_pi=(path_pi/toftime[j])/c;
	p_pi=TMath::Sqrt(tofmomx[j]*tofmomx[j]+tofmomy[j]*tofmomy[j]+tofmomz[j]*tofmomz[j]);
	seg_pi=tofseg[j];
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	total_count++;
	count_p++;
	tof_p=toftime[j];
	path_p=tofpath[j]/1000; //p pathlength(mm -> m)
	beta_p=(path_p/toftime[j])/c;
	p_p=TMath::Sqrt(tofmomx[j]*tofmomx[j]+tofmomy[j]*tofmomy[j]+tofmomz[j]*tofmomz[j]);
	seg_p=tofseg[j];
      }
      else if(tofpid[j]==11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==-11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==13&&tofedep[j]>ecut) count_mu++;
      else if(tofpid[j]==-13&&tofedep[j]>ecut) count_mu++;
      else if(tofpid[j]==22&&tofedep[j]>ecut) count_gamma++;
    }
    count_all=total_count+count_e+count_mu+count_gamma;
    hist_multi -> Fill(count_all); //Hodoscope multiplicity
    //hist_multi -> Fill(nhTof); //Hodoscope multiplicity
    if(total_count>1) event_w_2hits++; //counting multi-hit events
    if(count_p>0&&count_pi>0){ // selecting events w/ p,pi
      //hist_multi_pi -> Fill(total_count); //counting multi-hit events w/ p,pi
      hist_multi_pi -> Fill(count_all); //counting multi-hit events w/ p,pi
      event_w_2hits_wPpi++; //counting multi-hit events w/ p,pi

      // w/o detector resolution
      hist_tof_p -> Fill(tof_p); //tof of p and pi
      hist_tof_pi -> Fill(tof_pi);
      hist_pbeta -> Fill(1.0*p_p,beta_p);
      //hist_pbeta -> Fill(-1.0*p_pi,beta_pi);
      hist_pbeta -> Fill(1.0*p_pi,beta_pi);
      hist_ptof -> Fill(1.0*p_p,tof_p);
      //hist_ptof -> Fill(-1.0*p_pi,tof_pi);
      hist_ptof -> Fill(1.0*p_pi,tof_pi);
      // w/ detector resolution
      for(int j=0;j<6;j++){
	tof_ran_p=eventgen->Gaus(tof_p,det_resol[j]);
	tof_ran_pi=eventgen->Gaus(tof_pi,det_resol[j]);
	beta_p=(path_p/tof_ran_p)/c;
	beta_pi=(path_pi/tof_ran_pi)/c;
	hist_ran_p[j] -> Fill(tof_ran_p); //tof of p and pi
	hist_ran_pi[j]-> Fill(tof_ran_pi);

	hist_ran_pbeta[j] -> Fill(1.0*p_p,beta_p);
	//hist_ran_pbeta[j] -> Fill(-1.0*p_pi,beta_pi);
	hist_ran_pbeta[j] -> Fill(1.0*p_pi,beta_pi);
	hist_ran_ptof[j] -> Fill(1.0*p_p,tof_ran_p);
	//hist_ran_ptof[j] -> Fill(-1.0*p_pi,tof_ran_pi);
	hist_ran_ptof[j] -> Fill(1.0*p_pi,tof_ran_pi);
      }
    }
  }
  can_pid -> cd();
  hist_pid->Draw();
  hist_pid->GetXaxis()->SetTitle("PID(PDG encording) ");
  hist_pid->GetYaxis()->SetTitle("Counts");

  can_multi->cd();
  hist_multi->Draw();
  hist_multi_pi->Draw("same");
  hist_multi_pi->SetLineColor(kRed);
  hist_multi->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi->GetYaxis()->SetTitle("Counts");

  can_tof->cd();
  hist_tof_pi->Draw();
  hist_tof_p->Draw("same");
  hist_tof_pi->SetLineColor(kRed);
  hist_tof_pi->GetXaxis()->SetTitle("TOF(ns) ");
  hist_tof_pi->GetYaxis()->SetTitle("Counts");

  can_pbeta->cd();
  hist_pbeta->Draw("colz");
  hist_pbeta->GetXaxis()->SetTitle("p/q(GeV/c/q)");
  hist_pbeta->GetYaxis()->SetTitle("#beta");
  hist_pbeta->GetZaxis()->SetRangeUser(0,15);

  can_ptof->cd();
  hist_ptof->Draw("colz");
  hist_ptof->GetXaxis()->SetTitle("p/q(GeV/c/q)");
  hist_ptof->GetYaxis()->SetTitle("tof(ns)");
  hist_ptof->GetZaxis()->SetRangeUser(0,15);

  for(int i=0;i<6;i++){
    can_ran_tof[i]->cd();
    hist_ran_pi[i]->Draw();
    hist_ran_p[i]->Draw("same");
    hist_ran_pi[i]->SetLineColor(kRed);
    hist_ran_pi[i]->GetXaxis()->SetTitle("TOF(ns) ");
    hist_ran_pi[i]->GetYaxis()->SetTitle("Counts");

    can_ran_pbeta[i]->cd();
    hist_ran_pbeta[i]->Draw("colz");
    hist_ran_pbeta[i]->GetXaxis()->SetTitle("p/q(GeV/c/q)");
    hist_ran_pbeta[i]->GetYaxis()->SetTitle("#beta");
    hist_ran_pbeta[i]->GetZaxis()->SetRangeUser(0,15);

    can_ran_ptof[i]->cd();
    hist_ran_ptof[i]->Draw("colz");
    hist_ran_ptof[i]->GetXaxis()->SetTitle("p/q(GeV/c/q)");
    hist_ran_ptof[i]->GetYaxis()->SetTitle("tof(ns)");
    hist_ran_ptof[i]->GetZaxis()->SetRangeUser(0,15);
  }

  if(pdf){
    can_pid ->SaveAs("pdf/hodo_pid");
    can_multi ->SaveAs("pdf/hodo_multiplicity");
    can_tof ->SaveAs("pdf/hodo_tof");
    for(int i=0;i<6;i++){
      can_ran_tof[i] ->SaveAs(Form("pdf/hodo_tof_%d",i));
    }
    can_pbeta ->SaveAs("pdf/hodo_pbeta");
    for(int i=0;i<6;i++){
      can_ran_pbeta[i] ->SaveAs(Form("pdf/hodo_pbeta_%d",i));
    }
  }

  if(printout){
    std::cout<<""<<std::endl;
    std::cout<<"Event number  : "<<nevent<<std::endl;
    std::cout<<"2 hits        : "<<event_w_2hits<<std::endl;
    std::cout<<"2 hits w/ P pi: "<<event_w_2hits_wPpi<<std::endl;
    std::cout<<"Percentage    : "<<(double) event_w_2hits/nevent *100<<"("<<(double) event_w_2hits_wPpi/nevent *100<<") %"<<std::endl;
  }
}//end
