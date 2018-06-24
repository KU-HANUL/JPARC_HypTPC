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

/*
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooFFTConvPdf.h>
#include <RooPlot.h>
#include <RooFormula.h>

using namespace RooFit;
*/

static int class_e45 = 0;
static int evt = 0;

class e45{

 private:

  int nhTof;
  int event;
  int nEvt;
  int evtpid[100];
  double evtpx[100];
  double evtpy[100];
  double evtpz[100];
  double evtvx[100];
  double evtvy[100];
  double evtvz[100];

  int tofpid[100];
  int toftrid[100];
  int tofseg[100];
  double toftime[100];
  double tofedep[100];
  double tofmomx[100];
  double tofmomy[100];
  double tofmomz[100];
  double tofpath[100];
  double tofvtxx[100];
  double tofvtxy[100];
  double tofvtxz[100];
  double tofposx[100];
  double tofposy[100];
  double tofposz[100];
  int tofparentpid1[100];
  int tofparentpid2[100];
  int tofparentpid3[100];

  TFile *file = new TFile("../rootfile/e45/phsp/pi_plus/pipip/p1035_phsp.root","READ");
  TTree *tree = (TTree*)file->Get("tree");

  double R=305, L=4*70;
  double c=0.299792458;

  double ecut;
  double p_p,p_pi; //p,pi momentum
  double edep_p,edep_pi,edep_mu; //p,pi,mu deposit energy
  double tof_p,tof_pi,tof_mu; // time of flight of p,pi w/o detector resolution
  double beta_p,beta_pi; // beta of p,pi w/o detector resolution
  double path_p,path_pi; // path length of p,pi
  double tof_ran_p,tof_ran_pi; // time of flight of p,pi w/ detector resolution
  double beta_ran_p,beta_ran_pi; // beta of p,pi w/ detector resolutiona
  double vtxx_mu, vtxz_mu; // pi -> mu vertex point
  int event_w_2hits=0,event_w_2hits_wPpi=0; //flag for counting events with multi hits
  int total_count,count_all,count_p,count_pi,count_mu,count_e,count_gamma; //flag for counting particle
  int seg_p,seg_pi,seg_mu;

  int resol[6]={150,160,170,180,190,200}; //ps
  double det_resol[6]={0.150,0.160,0.170,0.180,0.190,0.200};

 public:

  void display_draw();
  void draw_evtgen();
  void draw_multiplicity();
  void draw_edep();
  void draw_elastic();

  e45() = default;

};

void e45::display_draw(){

  bool pdf=false;
  //bool pdf=true;
  //bool printout=false;
  bool printout=true;

  //gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("nhTof",&nhTof);

  tree->SetBranchAddress("evtpid",evtpid);
  tree->SetBranchAddress("evtpx",evtpx);
  tree->SetBranchAddress("evtpy",evtpy);
  tree->SetBranchAddress("evtpz",evtpz);
  tree->SetBranchAddress("evtvx",evtvx);
  tree->SetBranchAddress("evtvy",evtvy);
  tree->SetBranchAddress("evtvz",evtvz);

  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("toftime",toftime);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("tofmomx",tofmomx);
  tree->SetBranchAddress("tofmomy",tofmomy);
  tree->SetBranchAddress("tofmomz",tofmomz);
  tree->SetBranchAddress("tofpath",tofpath);
  tree->SetBranchAddress("tofposx",tofposx);
  tree->SetBranchAddress("tofposy",tofposy);
  tree->SetBranchAddress("tofposz",tofposz);
  tree->SetBranchAddress("tofvtxx",tofvtxx);
  tree->SetBranchAddress("tofvtxy",tofvtxy);
  tree->SetBranchAddress("tofvtxz",tofvtxz);
  tree->SetBranchAddress("tofparentpid1",tofparentpid1);
  tree->SetBranchAddress("tofparentpid2",tofparentpid2);
  tree->SetBranchAddress("tofparentpid3",tofparentpid3);

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
  TCanvas *can_seg = new TCanvas("can_seg","",1200,1200);
  can_seg -> Divide(2,1);
  TCanvas *can_path = new TCanvas("can_path","",1200,1200);
  can_path -> Divide(2,1);
  TCanvas *can_momentum  = new TCanvas("can_momentum","",1200,1200);
  can_momentum -> Divide(2,1);
  TCanvas *can_tofseg  = new TCanvas("can_tofseg","",1200,1200);
  can_tofseg -> Divide(2,1);

  ecut=1.0; //energy cut
  TH1D *hist_tof_p=new TH1D("hist_tof_p","P Time of flight",100,0,5);
  TH1D *hist_tof_pi=new TH1D("hist_tof_pi","pi Time of flight",100,0,5);
  TH1D *hist_seg_p=new TH1D("hist_seg_p","P segment id",32,0,32);
  TH1D *hist_seg_pi=new TH1D("hist_seg_pi","pi segment id",32,0,32);
  TH1D *hist_path_p=new TH1D("hist_path_p","P path id",1000,0,2);
  TH1D *hist_path_pi=new TH1D("hist_path_pi","pi path id",1000,0,2);
  TH1D *hist_p_p=new TH1D("hist_p_p","P momentum",1000,0,2);
  TH1D *hist_p_pi=new TH1D("hist_p_pi","pi momentum",1000,0,2);

  TH1D *hist_deltatof=new TH1D("hist_deltatof","#Delta Time of flight",100,0,5);
  TH2D *hist_pbeta=new TH2D("hist_pbeta","#beta vs p/q",500,-1.5,1.5,500,0.2,1.2);
  TH2D *hist_ptof=new TH2D("hist_ptof","TOF vs p/q",500,-1.5,1.5,500,0,20);
  TH2D *hist_tofseg_p=new TH2D("hist_tofseg_p","P, TOF vs seg",32,0,32,500,0,20);
  TH2D *hist_tofseg_pi=new TH2D("hist_tofseg_pi","pi, TOF vs seg",32,0,32,500,0,20);

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
      if(TMath::Abs(tofpid[j])==211&&tofedep[j]>ecut){ //for selecting events w/ pi
	total_count++;
	count_pi++;
	if(count_pi==1){
	  tof_pi=toftime[j];
	  path_pi=tofpath[j]/1000; //pi pathlength(mm -> m)
	  beta_pi=(path_pi/toftime[j])/c;
	  p_pi=TMath::Sqrt(tofmomx[j]*tofmomx[j]+tofmomy[j]*tofmomy[j]+tofmomz[j]*tofmomz[j]);
	  seg_pi=tofseg[j];
	}
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	total_count++;
	count_p++;
	if(count_p==1){
	  tof_p=toftime[j];
	  path_p=tofpath[j]/1000; //p pathlength(mm -> m)
	  beta_p=(path_p/toftime[j])/c;
	  p_p=TMath::Sqrt(tofmomx[j]*tofmomx[j]+tofmomy[j]*tofmomy[j]+tofmomz[j]*tofmomz[j]);
	  seg_p=tofseg[j];
	}
      }
      else if(tofpid[j]==11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==-11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==13&&tofedep[j]>ecut) count_mu++;
      else if(tofpid[j]==-13&&tofedep[j]>ecut) count_mu++;
      else if(tofpid[j]==22&&tofedep[j]>ecut) count_gamma++;
    }
    if(count_p>0&&count_pi>0){ // selecting events w/ p,pi
      hist_seg_p -> Fill(seg_p);
      hist_seg_pi -> Fill(seg_pi);
      hist_path_p -> Fill(path_p);
      hist_path_pi -> Fill(path_pi);
      hist_p_p -> Fill(p_p);
      hist_p_pi -> Fill(p_pi);

      // w/o detector resolution
      hist_tof_p -> Fill(tof_p); //tof of p and pi
      hist_tof_pi -> Fill(tof_pi);
      hist_pbeta -> Fill(1.0*p_p,beta_p);
      //hist_pbeta -> Fill(-1.0*p_pi,beta_pi);
      hist_pbeta -> Fill(1.0*p_pi,beta_pi);
      hist_ptof -> Fill(1.0*p_p,tof_p);
      //hist_ptof -> Fill(-1.0*p_pi,tof_pi);
      hist_ptof -> Fill(1.0*p_pi,tof_pi);
      hist_tofseg_p -> Fill(seg_p,tof_p);
      hist_tofseg_pi -> Fill(seg_pi,tof_pi);

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

  can_seg -> cd(1);
  hist_seg_p -> Draw();
  can_seg -> cd(2);
  hist_seg_pi -> Draw();

  can_path -> cd(1);
  hist_path_p -> Draw();
  can_path -> cd(2);
  hist_path_pi -> Draw();

  can_momentum -> cd(1);
  hist_p_p -> Draw();
  can_momentum -> cd(2);
  hist_p_pi -> Draw();

  can_tofseg -> cd(1);
  hist_tofseg_p -> Draw("colz");
  can_tofseg -> cd(2);
  hist_tofseg_pi -> Draw("colz");

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
  }

  class_e45++;
}

void e45::draw_evtgen(){

  bool pdf=false;
  //bool pdf=true;
  //bool printout=false;
  bool printout=true;

  //gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("nhTof",&nhTof);

  tree->SetBranchAddress("evtpid",evtpid);
  tree->SetBranchAddress("evtpx",evtpx);
  tree->SetBranchAddress("evtpy",evtpy);
  tree->SetBranchAddress("evtpz",evtpz);
  tree->SetBranchAddress("evtvx",evtvx);
  tree->SetBranchAddress("evtvy",evtvy);
  tree->SetBranchAddress("evtvz",evtvz);

  TCanvas *can_pid = new TCanvas("can_pid","",1200,1200);
  TCanvas *can_p = new TCanvas("can_p","",1200,1200);
  can_p -> Divide(1,2);

  TH1D *hist_pid =new TH1D("hist_pid","Evtgen PID;PDG encording;counts",4600,-2300,2300);
  TH1D *hist_p_p=new TH1D("hist_p_p","P momentum;P(GeV/c);Counts",1000,0,2);
  TH1D *hist_p_pi=new TH1D("hist_p_pi","pi momentum;P(GeV/c);Counts",1000,0,2);

  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    for(int j=0;j<nEvt;j++){
      hist_pid -> Fill(evtpid[j]); //PID
      if(TMath::Abs(evtpid[j])==211){ //for selecting events w/ pi
	p_pi=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]);
	hist_p_pi -> Fill(p_pi);
      }
      else if(evtpid[j]==2212){ //for selecting events w/ p
	p_p=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]);
	hist_p_p -> Fill(p_p);
      }
    }
  }

  can_pid -> cd();
  hist_pid->Draw();
  can_p -> cd(1);
  hist_p_p->Draw();
  can_p -> cd(2);
  hist_p_pi->Draw();

  if(pdf){
    can_pid ->SaveAs("pdf/evt_pid");
    can_p ->SaveAs("pdf/evt_p");
  }

  class_e45++;
}



void e45::draw_multiplicity(){

  bool pdf=false;
  //bool pdf=true;
  //bool printout=false;
  bool printout=true;

  //gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nhTof",&nhTof);

  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("toftrid",toftrid);

  TCanvas *can_multi = new TCanvas("can_multi","",1200,1200);
  TCanvas *can_pid = new TCanvas("can_pid","",1200,1200);

  TH1D *hist_pid =new TH1D("hist_pid","Hodoscope PID",4600,-2300,2300);
  TH1D *hist_multi =new TH1D("hist_multi","Hodoscope Multiplicity",7,0,7);
  TH1D *hist_multi_pi=new TH1D("hist_multi_pi","Hodoscope Multiplicity",7,0,7);

  int nevent=tree->GetEntries();
  ecut=0.1;
  double evt_pi=0, evt_mu;
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
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	total_count++;
	count_p++;
      }
      else if(tofpid[j]==11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==-11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==13&&tofedep[j]>ecut) count_mu++;
      else if(tofpid[j]==-13&&tofedep[j]>ecut) count_mu++;
      else if(tofpid[j]==22&&tofedep[j]>ecut) count_gamma++;
    }
    count_all=total_count+count_e+count_mu+count_gamma;
    hist_multi -> Fill(count_all); //Hodoscope multiplicity
    if(total_count>1) event_w_2hits++; //counting multi-hit events
    if(count_p>0&&count_pi>0){ // selecting events w/ p,pi
      hist_multi_pi -> Fill(count_all); //counting multi-hit events w/ p,pi
      event_w_2hits_wPpi++; //counting multi-hit events w/ p,pi
      evt_pi++;
    }
    if(count_p>0&&count_mu>0){
      evt_mu++;
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

  if(pdf){

  }

  if(printout){
    std::cout<<""<<std::endl;
    std::cout<<"Event number  : "<<nevent<<std::endl;
    std::cout<<"2 hits        : "<<event_w_2hits<<std::endl;
    std::cout<<"2 hits w/ P pi: "<<event_w_2hits_wPpi<<std::endl;
    std::cout<<"Percentage    : "<<(double) event_w_2hits/nevent *100<<"("<<(double) event_w_2hits_wPpi/nevent *100<<") %"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"pion %"<< 100*evt_pi/nevent<<std::endl;
    std::cout<<"mu %"<< 100*evt_mu/nevent<<std::endl;
  }

  class_e45++;
}



void e45::draw_edep(){

  bool pdf=false;
  //bool pdf=true;
  //bool printout=false;
  bool printout=true;

  //gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nhTof",&nhTof);

  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("toftime",toftime);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("tofmomx",tofmomx);
  tree->SetBranchAddress("tofmomy",tofmomy);
  tree->SetBranchAddress("tofmomz",tofmomz);
  tree->SetBranchAddress("tofpath",tofpath);
  tree->SetBranchAddress("tofposx",tofposx);
  tree->SetBranchAddress("tofposy",tofposy);
  tree->SetBranchAddress("tofposz",tofposz);
  tree->SetBranchAddress("tofvtxx",tofvtxx);
  tree->SetBranchAddress("tofvtxy",tofvtxy);
  tree->SetBranchAddress("tofvtxz",tofvtxz);
  tree->SetBranchAddress("tofparentpid1",tofparentpid1);
  tree->SetBranchAddress("tofparentpid2",tofparentpid2);
  tree->SetBranchAddress("tofparentpid3",tofparentpid3);

  TCanvas *can_edep = new TCanvas("can_edep","",1200,1200);
  can_edep -> Divide(2,2);
  TCanvas *can_vtx = new TCanvas("can_vtx","",1200,1200);
  can_vtx -> Divide(2,2);

  TH1D *hist_edep_p=new TH1D("hist_edep_p","P deposit energy;E(MeV);Counts",100,0,20);
  TH1D *hist_edep_pi=new TH1D("hist_edep_pi","pi deposit energy;E(MeV);Counts",100,0,15);
  TH1D *hist_edep_mu=new TH1D("hist_edep_mu","mu deposit energy;E(MeV);Counts",100,0,15);
  TH1D *hist_edep_pimu=new TH1D("hist_edep_pimu","pi to mu deposit energy;E(MeV);Counts",100,0,15);

  TH2D *hist_vtx_p=new TH2D("hist_vtx_p","P hit position;x(mm);z(mm)",800,-400,400,800,-400,400);
  TH2D *hist_vtx_pi=new TH2D("hist_vtx_pi","pi hit position;x(mm);z(mm)",800,-400,400,800,-400,400);
  TH2D *hist_vtx_mu=new TH2D("hist_vtx_mu","mu vertex;x(mm);z(mm)",800,-400,400,800,-400,400);
  TH2D *hist_vtx_pimu=new TH2D("hist_vtx_pimu","mu vertex;x(mm);z(mm)",800,-400,400,800,-400,400);

  int nevent=tree->GetEntries();
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);
    total_count=0, count_all=0, count_p=0, count_pi=0, count_e=0, count_mu=0, count_gamma=0;
    edep_p=0, edep_pi=0, edep_mu=0;
    for(int j=0;j<nhTof;j++){
      if(TMath::Abs(tofpid[j])==211&&tofedep[j]>ecut){ //for selecting events w/ pi
	count_pi++;
	edep_pi=edep_pi+tofedep[j];
	hist_vtx_pi-> Fill(tofposx[j],tofposz[j]);
	if(count_pi==1){
	  tof_pi=toftime[j];
	  seg_pi=tofseg[j];
	}
	else if(tof_pi<toftime[j]){
	  tof_pi=toftime[j];
	  seg_pi=tofseg[j];
	}
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	count_p++;
	edep_p=edep_p+tofedep[j];
	hist_vtx_p-> Fill(tofposx[j],tofposz[j]);
	if(count_p==1){
	  tof_p=toftime[j];
	  seg_p=tofseg[j];
	}
	else if(tof_p<toftime[j]){
	  tof_p=toftime[j];
	  seg_p=tofseg[j];
	}
      }
      else if(tofpid[j]==11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==-11&&tofedep[j]>ecut) count_e++;
      else if(tofpid[j]==22&&tofedep[j]>ecut) count_gamma++;
      else if(TMath::Abs(tofpid[j])==13&&tofedep[j]>ecut){
	count_mu++;
	edep_mu=edep_mu+tofedep[j];
	hist_vtx_mu-> Fill(tofvtxx[j],tofvtxz[j]);
	if(count_mu==1){
	  tof_mu=toftime[j];
	  seg_mu=tofseg[j];
	  vtxx_mu=tofvtxx[j];
	  vtxz_mu=tofvtxz[j];
	}
	else if(tof_mu>toftime[j]){
	  tof_mu=toftime[j];
	  seg_mu=tofseg[j];
	  vtxx_mu=tofvtxx[j];
	  vtxz_mu=tofvtxz[j];
	}
      }
    }
    if(count_p>0&&count_pi>0){ // selecting events w/ p,pi
      hist_edep_p -> Fill(edep_p);
      hist_edep_pi -> Fill(edep_pi);
    }
    if(count_mu>0)  hist_edep_mu-> Fill(edep_mu);
    if(count_p>0&&count_pi>0&&count_mu>0){
      if(seg_mu==seg_pi){
	hist_edep_pimu-> Fill(edep_mu);
	hist_vtx_pimu-> Fill(vtxx_mu,vtxz_mu);
      }
    }
  }

  can_edep -> cd(1);
  hist_edep_p -> Draw();
  can_edep -> cd(2);
  hist_edep_pi -> Draw();
  can_edep -> cd(3);
  hist_edep_mu -> Draw();
  hist_edep_pimu -> Draw("same");
  hist_edep_pimu -> SetLineColor(2);
  can_edep -> cd(4);
  hist_edep_mu -> Draw();
  hist_edep_pimu -> Draw("same");
  hist_edep_pimu -> SetLineColor(2);
  gPad -> SetLogy();
  can_vtx -> cd(1);
  hist_vtx_p -> Draw("colz");
  can_vtx -> cd(2);
  hist_vtx_pi -> Draw("colz");
  can_vtx -> cd(3);
  hist_vtx_mu -> Draw("colz");
  can_vtx -> cd(4);
  hist_vtx_pimu -> Draw("colz");

  if(pdf){
  }

  class_e45++;
}

void e45::draw_elastic(){

  bool pdf=false;
  //bool pdf=true;
  //bool printout=false;
  bool printout=true;

  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree");
  chain -> Add("../rootfile/e45/elastic/p1035_elastic.root");

  chain->SetBranchAddress("event",&event);
  chain->SetBranchAddress("nhTof",&nhTof);

  chain->SetBranchAddress("tofpid",tofpid);
  chain->SetBranchAddress("toftrid",toftrid);
  chain->SetBranchAddress("tofseg",tofseg);
  chain->SetBranchAddress("toftime",toftime);
  chain->SetBranchAddress("tofedep",tofedep);
  chain->SetBranchAddress("tofmomx",tofmomx);
  chain->SetBranchAddress("tofmomy",tofmomy);
  chain->SetBranchAddress("tofmomz",tofmomz);
  chain->SetBranchAddress("tofpath",tofpath);
  chain->SetBranchAddress("tofposx",tofposx);
  chain->SetBranchAddress("tofposy",tofposy);
  chain->SetBranchAddress("tofposz",tofposz);
  chain->SetBranchAddress("tofvtxx",tofvtxx);
  chain->SetBranchAddress("tofvtxy",tofvtxy);
  chain->SetBranchAddress("tofvtxz",tofvtxz);
  chain->SetBranchAddress("tofparentpid1",tofparentpid1);
  chain->SetBranchAddress("tofparentpid2",tofparentpid2);
  chain->SetBranchAddress("tofparentpid3",tofparentpid3);

  TCanvas *can_beamseg = new TCanvas("can_beamseg","",1200,1200);
  TCanvas *can_beamhit  = new TCanvas("can_beamhit","",1200,1200);
  TCanvas *can_beampid = new TCanvas("can_beampid","",1200,1200);
  TCanvas *can_beamp = new TCanvas("can_beamp","",1200,1200);

  TH1D *hist_pid_beam = new TH1D("hist_pid_beam","Evtgen PID;PDG encording;counts",4600,-2300,2300);
  TH1D *hist_seg_beam = new TH1D("hist_seg_beam","segment;# of Seg;Counts",32,0,32);
  TH1D *hist_p_beam = new TH1D("hist_p_beam","pi momentum;GeV/c;Counts",100,0,2);

  TH2D *hist_hitpos_beam = new TH2D("hist_hitpos_beam","pi hit position;x(mm);z(mm)",800,-400,400,800,-400,400);

  int nevent=chain->GetEntries();
  for(int i=0;i<nevent;i++){
    chain->GetEntry(i);
    for(int j=0;j<nhTof;j++){
      hist_pid_beam -> Fill(tofpid[j]);
      if(TMath::Abs(tofpid[j])==211){ //for selecting events w/ pi
	hist_seg_beam -> Fill(tofseg[j]);
	hist_hitpos_beam -> Fill(tofposx[j],tofposz[j]);
	p_pi=TMath::Sqrt(tofmomx[j]*tofmomx[j]+tofmomy[j]*tofmomy[j]+tofmomz[j]*tofmomz[j]);
	hist_p_beam -> Fill(p_pi);
      }
    }
  }

  can_beampid -> cd();
  hist_pid_beam->Draw();
  can_beamseg -> cd();
  hist_seg_beam -> Draw();
  can_beamhit -> cd();
  hist_hitpos_beam ->Draw("colz");
  can_beamp -> cd();
  hist_p_beam ->Draw();

  if(pdf){
  }

  class_e45++;
}
