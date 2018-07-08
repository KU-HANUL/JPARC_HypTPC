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

class hodo3{

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
  void pipin2(TFile *file, double ecut, bool window, double &acceptance, bool printout=false);
  hodo3() = default;

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

void hodo3::pipip(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

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
  hist_pid[class_hodo] =new TH1D(Form("hist_pid_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);
  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_multi_pi[class_hodo]=new TH1D(Form("hist_multi_pi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);

  int nevent=tree->GetEntries();

  double tan_p, tan_pi, mu_vtx;

  int total_count,count_p,count_pi,count_mu,count_e,count_gamma; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi, flag_seg_p;
  int hodo_seg[38]={0};
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
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  flag_seg_pi=32+0;
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  flag_seg_pi=32+1;
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 2;
	  flag_seg_pi=32+2;
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 3;
	  flag_seg_pi=32+3;
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 4;
	  flag_seg_pi=32+4;
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 5;
	  flag_seg_pi=32+5;
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  flag_seg_pi=tofseg[j];
	}
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	total_count++;
	count_p++;
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  flag_seg_p=32+0;
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  flag_seg_p=32+1;
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 2;
	  flag_seg_p=32+2;
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 3;
	  flag_seg_p=32+3;
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 4;
	  flag_seg_p=32+4;
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 5;
	  flag_seg_p=32+5;
	}
	else{hodo_seg[tofseg[j]] += 1;
	  flag_seg_p=tofseg[j];
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
      }
    }
    for(int i=0;i<38;i++){
      if(hodo_seg[i]>0) multiplicity++;
      check_flag+=hodo_seg[i];
      hodo_seg[i]=0;
    }
    hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
    if(multiplicity==1){
      if(count_p==1 && count_pi==1 && flag_seg_pi==flag_seg_p) count_m1_coin++;
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
    can_multi[class_hodo] -> Close();
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
    std::cout<<"multi-0 ratio  : "<<(double) (count_m0)/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-1 ratio  : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-2 ratio  : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

    std::cout<<"m1 coin/m1     : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
    std::cout<<"m2 coin/m2     : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

    std::cout<<""<<std::endl;
  }

  class_hodo++;
}


void hodo3::pipin(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

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
  int hodo_seg[38]={0};

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
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  flag_seg_pi_p=32+0;
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  flag_seg_pi_p=32+1;
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 2;
	  flag_seg_pi_p=32+2;
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 3;
	  flag_seg_pi_p=32+3;
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 4;
	  flag_seg_pi_p=32+4;
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 5;
	  flag_seg_pi_p=32+5;
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  flag_seg_pi_p=tofseg[j];
	}
      }
      else if(tofpid[j]==-211&&tofedep[j]>ecut){ //for selecting events w/ pi-
	total_count++;
	count_pi_m++;
	hodo_seg[tofseg[j]] += 1;
	flag_seg_pi_m=tofseg[j];
	total_count++;
	count_pi_m++;
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  flag_seg_pi_m=32+0;
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  flag_seg_pi_m=32+1;
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 2;
	  flag_seg_pi_m=32+2;
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 3;
	  flag_seg_pi_m=32+3;
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 4;
	  flag_seg_pi_m=32+4;
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 5;
	  flag_seg_pi_m=32+5;
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  flag_seg_pi_m=tofseg[j];
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
      }
    }
    for(int i=0;i<38;i++){
      if(hodo_seg[i]>0) multiplicity++;
      check_flag+=hodo_seg[i];
      hodo_seg[i]=0;
    }
    hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
    if(multiplicity==1){
      if(count_pi_p==1 && count_pi_m==1 && flag_seg_pi_p==flag_seg_pi_m) count_m1_coin++;
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
    std::cout<<"multi-0 ratio  : "<<(double) (count_m0)/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-1 ratio  : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-2 ratio  : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

    std::cout<<"m1 coin/m1     : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
    std::cout<<"m2 coin/m2     : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

    std::cout<<""<<std::endl;
  }

  class_hodo++;
}



void hodo3::pipin2(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

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
  tree->SetBranchAddress("toftrid",toftrid);

  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("evtpid",evtpid);
  tree->SetBranchAddress("evttrid",evttrid);

  can_multi[class_hodo] = new TCanvas(Form("can_multi_%d",class_hodo),"",1200,1200);
  can_pid[class_hodo] = new TCanvas(Form("can_pid_%d",class_hodo),"",1200,1200);
  hist_pid[class_hodo] =new TH1D(Form("hist_pid_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);
  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_multi_pi[class_hodo]=new TH1D(Form("hist_multi_pi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  int nevent=tree->GetEntries();

  int total_count,count_pi1,count_pi2,count_mu,count_e,count_gamma; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi1, flag_seg_pi2;
  int hodo_seg[38]={0};

  int count_m0=0;
  int count_m1=0, count_m1_coin=0, count_m1_pi1=0, count_m1_pi2=0, count_m1_dummy=0;
  int count_m2=0, count_m2_coin=0, count_m2_multi=0, count_m2_pi1=0, count_m2_pi2=0, count_m2_dummy=0;

  int dummy=0;
  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);

    total_count=0;
    count_pi1=0;
    count_pi2=0;
    count_e=0;
    count_mu=0;
    multiplicity=0;
    check_flag=0;

    flag_seg_pi1, flag_seg_pi2;

    for(int j=0;j<nhTof;j++){
      hist_pid[class_hodo] -> Fill(tofpid[j]); //PID
      if(window==true&&tofseg[j]==0&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==1&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==2&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==3&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==17&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==18&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(toftrid[j]==1&&tofpid[j]==211&&tofedep[j]>ecut){ //for selecting events w/ 1st pi+
	total_count++;
	count_pi1++;
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  flag_seg_pi1=32+0;
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  flag_seg_pi1=32+1;
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 2;
	  flag_seg_pi1=32+2;
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 3;
	  flag_seg_pi1=32+3;
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 4;
	  flag_seg_pi1=32+4;
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 5;
	  flag_seg_pi1=32+5;
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  flag_seg_pi1=tofseg[j];
	}
      }
      else if(toftrid[j]==2&&tofpid[j]==211&&tofedep[j]>ecut){ //for selecting events w/ 2nd pi+
	total_count++;
	count_pi2++;
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  flag_seg_pi2=32+0;
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  flag_seg_pi2=32+1;
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 2;
	  flag_seg_pi2=32+2;
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 3;
	  flag_seg_pi2=32+3;
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 4;
	  flag_seg_pi2=32+4;
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 5;
	  flag_seg_pi2=32+5;
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  flag_seg_pi2=tofseg[j];
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
      }
    }
    for(int i=0;i<38;i++){
      if(hodo_seg[i]>0) multiplicity++;
      check_flag+=hodo_seg[i];
      hodo_seg[i]=0;
    }
    hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
    if(multiplicity==1){
      if(count_pi1==1 && count_pi2==1 && flag_seg_pi1==flag_seg_pi2) count_m1_coin++;
      if(count_pi1==1 && count_pi2==0) count_m1_pi1++;
      if(count_pi2==1 && count_pi1==0) count_m1_pi2++;
      if(count_pi2==0 && count_pi1==0) count_m1_dummy++;
      count_m1++;
    }
    else if(multiplicity>1){
      if(count_pi1==1 && count_pi2==1 && flag_seg_pi1==flag_seg_pi2) count_m2_coin++;
      if(count_pi1>0 && count_pi2>0) count_m2_multi++;
      if(count_pi1>0 && count_pi2==0) count_m2_pi1++;
      if(count_pi2>0 && count_pi1==0) count_m2_pi2++;
      if(count_pi2==0 && count_pi1==0) count_m2_dummy++;
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
    std::cout<<"multi-1 one pi+ hits : "<<count_m1_pi1+count_m1_pi2<<std::endl;
    std::cout<<"multi-1 dummy  : "<<count_m1_dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-2 number : "<<count_m2<<std::endl;
    std::cout<<"multi-2 coin   : "<<count_m2_coin<<std::endl;
    std::cout<<"multi-2 pipi hits : "<<count_m2_multi<<std::endl;
    std::cout<<"multi-2 one pi+ hits : "<<count_m2_pi1+count_m2_pi2<<std::endl;
    std::cout<<"multi-2 dummy  : "<<count_m2_dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"acceptance     : "<<(double) acceptance*100<< " %"<<std::endl;
    std::cout<<"multi-0 ratio  : "<<(double) (count_m0)/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-1 ratio  : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-2 ratio  : "<<(double) count_m2/nevent*100<< " %"<<std::endl;
    std::cout<<"m1 coin/m1     : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
    std::cout<<"m2 coin/m2     : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

    std::cout<<""<<std::endl;
  }

  class_hodo++;
}
