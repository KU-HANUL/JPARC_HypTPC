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

static int class_hodo = 0;

class hodo{

 private:

  double window_size = 40.0/2; //half of window size
  double mass_p=0.938272;
  double mass_pi=0.139570;

  int nhTof;
  int event;
  int nEvt;
  int evtpid[100];
  int evttrid[100];
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

  int nhTarget[100];
  int targettrid[100];
  double targetedep[100];

 public:

  void pipip(TFile *file, double ecut, bool window, double &acceptance, bool printout=false);
  void pipin(TFile *file, double ecut, bool window, double &acceptance, bool printout=false);
  hodo() = default;

  TCanvas *can_multi[100];
  TCanvas *can_pid[100];
  TCanvas *can_tan[100];
  TCanvas *can_pimu[100];
  TCanvas *can_ph[100];
  TCanvas *can_target[100];

  TH1D *hist_pid[100];
  TH1D *hist_multi[100];
  TH1D *hist_multi_pi[100];
  TH1D *hist_tan_p1[100];
  TH1D *hist_tan_pi1[100];
  TH1D *hist_tan_p2[100];
  TH1D *hist_tan_pi2[100];

  TH1D *hist_evtpid[100];
  TH1D *hist_mu_vtx[100];

  TH1D *hist_ph_pi1[100];
  TH1D *hist_ph_p1[100];
  TH1D *hist_ph_pi2[100];
  TH1D *hist_ph_p2[100];

  TH1D *hist_target_pi[100];
  TH1D *hist_target_p[100];

};

void hodo::pipip(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

  TTree *tree = (TTree*)file->Get("tree");

  bool pdf=false;
  //bool pdf=true;

  //gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);

  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("tofposx",tofposx);
  tree->SetBranchAddress("tofposy",tofposy);
  tree->SetBranchAddress("tofposz",tofposz);
  tree->SetBranchAddress("tofmomx",tofmomx);
  tree->SetBranchAddress("tofmomy",tofmomy);
  tree->SetBranchAddress("tofmomz",tofmomz);
  tree->SetBranchAddress("tofvtxx",tofvtxx);
  tree->SetBranchAddress("tofvtxy",tofvtxy);
  tree->SetBranchAddress("tofvtxz",tofvtxz);

  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("evtpid",evtpid);
  tree->SetBranchAddress("evttrid",evttrid);
  tree->SetBranchAddress("evtpx",evtpx);
  tree->SetBranchAddress("evtpy",evtpy);
  tree->SetBranchAddress("evtpz",evtpz);

  tree->SetBranchAddress("nhTarget",&nhTarget);
  tree->SetBranchAddress("targetedep",targetedep);
  tree->SetBranchAddress("targettrid",targettrid);

  can_multi[class_hodo] = new TCanvas(Form("can_multi_%d",class_hodo),"",1200,1200);
  can_pid[class_hodo] = new TCanvas(Form("can_pid_%d",class_hodo),"",1200,1200);
  can_pid[class_hodo] -> Divide(2,1);
  can_tan[class_hodo] = new TCanvas(Form("can_tan_%d",class_hodo),"",1200,1200);
  can_tan[class_hodo] -> Divide(2,1);
  can_pimu[class_hodo] = new TCanvas(Form("can_pimu_%d",class_hodo),"",1200,1200);
  can_ph[class_hodo] = new TCanvas(Form("can_ph_%d",class_hodo),"",1200,1200);
  can_ph[class_hodo] -> Divide(2,1);
  can_target[class_hodo] = new TCanvas(Form("can_target_%d",class_hodo),"",1200,1200);
  can_target[class_hodo]->Divide(2,1);

  hist_evtpid[class_hodo] =new TH1D(Form("hist_evtpid_%d",class_hodo),"EVTGEN PID",4600,-2300,2300);
  hist_pid[class_hodo] =new TH1D(Form("hist_pid_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);
  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_multi_pi[class_hodo]=new TH1D(Form("hist_multi_pi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_tan_pi1[class_hodo]=new TH1D(Form("hist_tan_pi1_%d",class_hodo),"tan#theta of pi",1000,-5,5);
  hist_tan_p1[class_hodo]=new TH1D(Form("hist_tan_p1_%d",class_hodo),"tan#theta of p",1000,-5,5);
  hist_tan_pi2[class_hodo]=new TH1D(Form("hist_tan_pi2_%d",class_hodo),"tan#theta of pi",1000,-5,5);
  hist_tan_p2[class_hodo]=new TH1D(Form("hist_tan_p2_%d",class_hodo),"tan#theta of p",1000,-5,5);
  hist_mu_vtx[class_hodo]=new TH1D(Form("hist_mu_vt_%d",class_hodo),"origin to #pi-> #mu vtx",5000,0,500);

  hist_ph_pi1[class_hodo]=new TH1D(Form("hist_ph_pi1_%d",class_hodo),"Ph, pi",1000,0,2.5);
  hist_ph_p1[class_hodo]=new TH1D(Form("hist_ph_p1_%d",class_hodo),"Ph, P",1000,0,2.5);
  hist_ph_pi2[class_hodo]=new TH1D(Form("hist_ph_pi2_%d",class_hodo),"Ph, pi",1000,0,2.5);
  hist_ph_p2[class_hodo]=new TH1D(Form("hist_ph_p2_%d",class_hodo),"Ph, P",1000,0,2.5);


  hist_target_pi[class_hodo]=new TH1D(Form("hist_target_pi_%d",class_hodo),"pi,kE - target edep",2500,0,2.5);
  hist_target_p[class_hodo]=new TH1D(Form("hist_target_p_%d",class_hodo),"p,kE - target edep",2500,0,2.5);
  int nevent=tree->GetEntries();

  double tan_p, tan_pi, mu_vtx;

  int total_count,count_p,count_pi,count_mu,count_e,count_gamma; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi, flag_seg_p;
  int hodo_seg[32]={0};
  int trid;
  int pimu=0;
  int count_m0=0;
  int count_m1=0, count_m1_coin=0, count_m1_p=0, count_m1_pi=0, count_m1_dummy=0;
  int count_m2=0, count_m2_coin=0, count_m2_multi=0, count_m2_p=0, count_m2_pi=0, count_m2_dummy=0, count_m2_tan=0;

  int dummy=0;
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);

    total_count=0;
    count_p=0;
    count_pi=0;
    count_e=0;
    count_mu=0;
    multiplicity=0;
    check_flag=0;
    tan_p=0, tan_pi=0;
    flag_seg_pi, flag_seg_p;

    for(int j=0;j<nEvt;j++){
      hist_evtpid[class_hodo] -> Fill(evtpid[j]);
      double pv=evtpy[j];
      double ph=TMath::Sqrt(evtpx[j]*evtpx[j]+evtpz[j]*evtpz[j]);
      if(TMath::Abs(evtpid[j])==211){
	double kE=TMath::Sqrt(pv*pv+ph*ph+mass_pi*mass_pi) - mass_pi;
	/*
	for(int k=0;k<nhTarget;k++){
	  if(targettrid[k]==trid){
	    hist_target_pi[class_hodo] -> Fill(kE-targetedep[k]);
	  }
	}
	*/
	hist_ph_pi1[class_hodo]-> Fill(ph);

	tan_pi=(double) pv/ph;
	hist_tan_pi1[class_hodo] -> Fill(tan_pi);
	trid=evttrid[j];

	for(int k=0;k<nhTof;k++){
	  if(toftrid[k]==trid){
	    hist_tan_pi2[class_hodo] -> Fill(tan_pi);
	  }
	}
      }
      else if(TMath::Abs(evtpid[j])==2212){
	double kE=TMath::Sqrt(pv*pv+ph*ph+mass_p*mass_p) - mass_p;
	/*
	for(int k=0;k<nhTarget;k++){
	  if(targettrid[k]==trid){
	    hist_target_p[class_hodo] -> Fill(kE-targetedep[k]);
	  }
	}
	*/
	hist_ph_p1[class_hodo]-> Fill(ph);

	tan_p=(double) pv/ph;
	hist_tan_p1[class_hodo] -> Fill(tan_p);
	trid=evttrid[j];
	for(int k=0;k<nhTof;k++){
	  if(toftrid[k]==trid){
	    hist_tan_p2[class_hodo] -> Fill(tan_p);
	  }
	}

      }
    }
    for(int j=0;j<nhTof;j++){
      hist_pid[class_hodo] -> Fill(tofpid[j]); //PID
      if(window==true&&tofseg[j]==0&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==1&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==2&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==3&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==17&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==18&&TMath::Abs(tofposy[j])<window_size) dummy++;

      else if(TMath::Abs(tofpid[j])==211&&tofedep[j]>ecut){ //for selecting events w/ pi
	total_count++;
	count_pi++;
	hodo_seg[tofseg[j]] += 1;
	flag_seg_pi=tofseg[j];
	trid=toftrid[j];
	for(int k=0;k<nEvt;k++){
	  if(evttrid[k]==trid){
	    double ph=TMath::Sqrt(evtpx[k]*evtpx[k]+evtpz[k]*evtpz[k]);
	    hist_ph_pi2[class_hodo]-> Fill(ph);
	  }
	}
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	total_count++;
	count_p++;
	hodo_seg[tofseg[j]] += 1;
	flag_seg_p=tofseg[j];
	trid=toftrid[j];
	for(int k=0;k<nEvt;k++){
	  if(evttrid[k]==trid){
	    double ph=TMath::Sqrt(evtpx[k]*evtpx[k]+evtpz[k]*evtpz[k]);
	    hist_ph_p2[class_hodo]-> Fill(ph);
	  }
	}
      }
      else if(TMath::Abs(tofpid[j])==11&&tofedep[j]>ecut){
	total_count++;
	count_e++;
	hodo_seg[tofseg[j]] += 1;
      }
      else if(TMath::Abs(tofpid[j])==13&&tofedep[j]>ecut){
	total_count++;
	count_mu++;
	hodo_seg[tofseg[j]] += 1;
	mu_vtx=TMath::Sqrt(tofvtxx[j]*tofvtxx[j]+tofvtxz[j]*tofvtxz[j]);
	hist_mu_vtx[class_hodo] -> Fill(mu_vtx);
      }
    }
    for(int i=0;i<32;i++){
      if(hodo_seg[i]>0) multiplicity++;
      check_flag+=hodo_seg[i];
      hodo_seg[i]=0;
    }
    hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
    if(multiplicity==1){
      if(count_p==1 && count_pi==1) count_m1_coin++;
      if(count_p==1 && count_pi==0) count_m1_p++;
      if(count_pi==1 && count_p==0) count_m1_pi++;
      if(count_pi==0 && count_p==0) count_m1_dummy++;
      count_m1++;
    }
    else if(multiplicity>1){
      if(count_p==1 && count_pi==1 && flag_seg_pi==flag_seg_p) count_m2_coin++;
      if(count_p>0 && count_pi>0) count_m2_multi++;
      if(count_p>0 && count_pi==0){
	count_m2_p++;
      }
      if(count_pi>0 && count_p==0){
	count_m2_pi++;
      }
      if(count_pi==0 && count_p==0) count_m2_dummy++;
      count_m2++;
    }
    else count_m0++;
  }
  can_pid[class_hodo] -> cd(1);
  hist_pid[class_hodo]->Draw();
  hist_pid[class_hodo]->GetXaxis()->SetTitle("Hodoscope PID(PDG encording) ");
  hist_pid[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_pid[class_hodo] -> cd(2);
  hist_evtpid[class_hodo]->Draw();
  hist_evtpid[class_hodo]->GetXaxis()->SetTitle("EVT PID(PDG encording) ");
  hist_evtpid[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_multi[class_hodo]->cd();
  hist_multi[class_hodo]->Draw();
  hist_multi[class_hodo]->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_tan[class_hodo]->cd(1);
  hist_tan_p1[class_hodo]->Draw();
  hist_tan_p1[class_hodo]->GetXaxis()->SetTitle("tan#theta of p ");
  hist_tan_p1[class_hodo]->GetYaxis()->SetTitle("Counts");
  hist_tan_p2[class_hodo]->Draw("same");
  hist_tan_p2[class_hodo]->SetLineColor(2);

  can_tan[class_hodo]->cd(2);
  hist_tan_pi1[class_hodo]->Draw();
  hist_tan_pi1[class_hodo]->GetXaxis()->SetTitle("tan#theta of pi ");
  hist_tan_pi1[class_hodo]->GetYaxis()->SetTitle("Counts");
  hist_tan_pi2[class_hodo]->Draw("same");
  hist_tan_pi2[class_hodo]->SetLineColor(2);

  can_pimu[class_hodo]->cd();
  hist_mu_vtx[class_hodo]->Draw();
  hist_mu_vtx[class_hodo]->GetXaxis()->SetTitle("length(mm) ");
  hist_mu_vtx[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_ph[class_hodo]->cd(1);
  hist_ph_pi1[class_hodo]->Draw();
  hist_ph_pi2[class_hodo]->Draw("same");
  hist_ph_pi2[class_hodo]->SetLineColor(2);

  can_ph[class_hodo]->cd(2);
  hist_ph_p1[class_hodo]->Draw();
  hist_ph_p2[class_hodo]->Draw("same");
  hist_ph_p2[class_hodo]->SetLineColor(2);

  can_target[class_hodo]->cd(1);
  hist_target_pi[class_hodo]->Draw();

  can_target[class_hodo]->cd(2);
  hist_target_p[class_hodo]->Draw();

  if(printout==false){
    can_pid[class_hodo] -> Close();
    can_multi[class_hodo] -> Close();
    can_tan[class_hodo] -> Close();
  }
  if(pdf){

  }
  acceptance= (double) count_m2_multi/nevent;

  if(printout){
    std::cout<<""<<std::endl;
    std::cout<<"Event number   : "<<nevent<<std::endl;
    std::cout<<"hit window     : "<<dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-0 number : "<<(double) count_m0<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-1 number : "<<count_m1<<std::endl;
    std::cout<<"multi-1 coin   : "<<count_m1_coin<<std::endl;
    std::cout<<"multi-1 P hits : "<<count_m1_p<<std::endl;
    std::cout<<"multi-1 pi hits: "<<count_m1_pi<<std::endl;
    std::cout<<"multi-1 dummy  : "<<count_m1_dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-2 number : "<<count_m2<<std::endl;
    std::cout<<"multi-2 coin   : "<<count_m2_coin<<std::endl;
    std::cout<<"multi-2 Ppi hits : "<<count_m2_multi<<std::endl;
    std::cout<<"multi-2 P hits : "<<count_m2_p<<std::endl;
    std::cout<<"multi-2 pi hits: "<<count_m2_pi<< " tan theta > 1.14 "<<count_m2_tan<<std::endl;

    std::cout<<"multi-2 dummy  : "<<count_m2_dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"acceptance     : "<<(double) acceptance*100<< " %"<<std::endl;
    std::cout<<"multi-0 ratio  : "<<(double) (count_m0+dummy)/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-1 ratio  : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-2 ratio  : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

    std::cout<<"m1 coin/m1     : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
    std::cout<<"m2 coin/m2     : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

    std::cout<<""<<std::endl;
  }

  class_hodo++;
}


void hodo::pipin(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

  TTree *tree = (TTree*)file->Get("tree");

  bool pdf=false;
  //bool pdf=true;

  //gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("tofposx",tofposx);
  tree->SetBranchAddress("tofposy",tofposy);
  tree->SetBranchAddress("tofposz",tofposz);
  //tree->SetBranchAddress("toftrid",toftrid);


  can_multi[class_hodo] = new TCanvas(Form("can_multi_%d",class_hodo),"",1200,1200);
  can_pid[class_hodo] = new TCanvas(Form("can_pid_%d",class_hodo),"",1200,1200);
  hist_pid[class_hodo] =new TH1D(Form("hist_pid_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);
  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_multi_pi[class_hodo]=new TH1D(Form("hist_multi_pi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  int nevent=tree->GetEntries();

  int total_count,count_pi_p,count_pi_m,count_mu,count_e,count_gamma; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi_p, flag_seg_pi_m;
  int hodo_seg[32]={0};

  int count_m0=0;
  int count_m1=0, count_m1_coin=0, count_m1_pi_p=0, count_m1_pi_m=0, count_m1_dummy=0;
  int count_m2=0, count_m2_coin=0, count_m2_multi=0, count_m2_pi_p=0, count_m2_pi_m=0, count_m2_dummy=0;

  int dummy=0;
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);

    total_count=0;
    count_pi_p=0;
    count_pi_m=0;
    count_e=0;
    count_mu=0;
    multiplicity=0;
    check_flag=0;

    flag_seg_pi_p, flag_seg_pi_m;

    for(int j=0;j<nhTof;j++){
      hist_pid[class_hodo] -> Fill(tofpid[j]); //PID
      if(window==true&&tofseg[j]==0&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==1&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==2&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==3&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==17&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==18&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(tofpid[j]==211&&tofedep[j]>ecut){ //for selecting events w/ pi+
	total_count++;
	count_pi_p++;
	hodo_seg[tofseg[j]] += 1;
	flag_seg_pi_p=tofseg[j];
      }
      else if(tofpid[j]==-211&&tofedep[j]>ecut){ //for selecting events w/ pi-
	total_count++;
	count_pi_m++;
	hodo_seg[tofseg[j]] += 1;
	flag_seg_pi_m=tofseg[j];
      }
      else if(TMath::Abs(tofpid[j])==11&&tofedep[j]>ecut){
	total_count++;
	count_e++;
	hodo_seg[tofseg[j]] += 1;
      }
      else if(TMath::Abs(tofpid[j])==13&&tofedep[j]>ecut){
	total_count++;
	count_mu++;
	hodo_seg[tofseg[j]] += 1;
      }
    }
    for(int i=0;i<32;i++){
      if(hodo_seg[i]>0) multiplicity++;
      check_flag+=hodo_seg[i];
      hodo_seg[i]=0;
    }
    hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
    if(multiplicity==1){
      if(count_pi_p==1 && count_pi_m==1) count_m1_coin++;
      if(count_pi_p==1 && count_pi_m==0) count_m1_pi_p++;
      if(count_pi_m==1 && count_pi_p==0) count_m1_pi_m++;
      if(count_pi_m==0 && count_pi_p==0) count_m1_dummy++;
      count_m1++;
    }
    else if(multiplicity>1){
      if(count_pi_p==1 && count_pi_m==1 && flag_seg_pi_p==flag_seg_pi_m) count_m2_coin++;
      if(count_pi_p>0 && count_pi_m>0) count_m2_multi++;
      if(count_pi_p>0 && count_pi_m==0) count_m2_pi_p++;
      if(count_pi_m>0 && count_pi_p==0) count_m2_pi_m++;
      if(count_pi_m==0 && count_pi_p==0) count_m2_dummy++;
      count_m2++;
    }
    else count_m0++;
  }
  can_pid[class_hodo] -> cd();
  hist_pid[class_hodo]->Draw();
  hist_pid[class_hodo]->GetXaxis()->SetTitle("PID(PDG encording) ");
  hist_pid[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_multi[class_hodo]->cd();
  hist_multi[class_hodo]->Draw();
  hist_multi[class_hodo]->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi[class_hodo]->GetYaxis()->SetTitle("Counts");

  if(printout==false){
    can_pid[class_hodo] -> Close();
    can_multi[class_hodo]->Close();
  }
  if(pdf){

  }
  acceptance= (double) count_m2_multi/nevent;

    if(printout){
    std::cout<<""<<std::endl;
    std::cout<<"Event number   : "<<nevent<<std::endl;
    std::cout<<"hit window     : "<<dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-0 number : "<<(double) count_m0<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-1 number : "<<count_m1<<std::endl;
    std::cout<<"multi-1 coin   : "<<count_m1_coin<<std::endl;
    std::cout<<"multi-1 pi+ hits : "<<count_m1_pi_p<<std::endl;
    std::cout<<"multi-1 pi- hits: "<<count_m1_pi_m<<std::endl;
    std::cout<<"multi-1 dummy  : "<<count_m1_dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-2 number : "<<count_m2<<std::endl;
    std::cout<<"multi-2 coin   : "<<count_m2_coin<<std::endl;
    std::cout<<"multi-2 pipi hits : "<<count_m2_multi<<std::endl;
    std::cout<<"multi-2 pi+ hits : "<<count_m2_pi_p<<std::endl;
    std::cout<<"multi-2 pi- hits: "<<count_m2_pi_m<<std::endl;
    std::cout<<"multi-2 dummy  : "<<count_m2_dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"acceptance     : "<<(double) acceptance*100<< " %"<<std::endl;
    std::cout<<"multi-0 ratio  : "<<(double) (count_m0+dummy)/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-1 ratio  : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-2 ratio  : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

    std::cout<<"m1 coin/m1     : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
    std::cout<<"m2 coin/m2     : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

    std::cout<<""<<std::endl;
  }

  class_hodo++;
}
