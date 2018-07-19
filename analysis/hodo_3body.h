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

class hodo_3body{

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

  double GetCos(TLorentzVector lv1, TLorentzVector lv2);
  double GetCos_xz(TLorentzVector lv1, TLorentzVector lv2);
  //double GetCos_yz(TLorentzVector lv1, TLorentzVector lv2);
  double GetCos_yz(TLorentzVector lv1);
  double GetAng(TLorentzVector lv1, TLorentzVector lv2);
  double GetAng_xz(TLorentzVector lv1, TLorentzVector lv2);
  //double GetAng_yz(TLorentzVector lv1, TLorentzVector lv2);
  double GetAng_yz(TLorentzVector lv1);

  void pipip(TFile *file, double ecut, bool window, double &acceptance, bool printout=false);
  void pipin(TFile *file, double ecut, bool window, double &acceptance, bool printout=false);
  void pipin2(TFile *file, double ecut, bool window, double &acceptance, bool printout=false);

  hodo_3body() = default;

  TCanvas *can_multi[100];
  TCanvas *can_pid[100];
  TCanvas *can_deltatof[100];
  TCanvas *can_hitpattern2D[100];
  TCanvas *can_segpattern2D[100];
  TCanvas *can_xzpattern2D[100];
  TCanvas *can_yzpattern2D[100];
  TCanvas *can_xzpattern2D2[100];
  TCanvas *can_yzpattern2D2[100];

  TH1D *hist_pid[100];
  TH1D *hist_deltatof[100];
  TH1D *hist_multi[100];
  TH1D *hist_evtpid[100];

  TH2D *hist_hitpattern2D_pi[100];
  TH2D *hist_hitpattern2D_p[100];

  TH2D *hist_segpattern2D[100];
  TH2D *hist_xzpattern2D[100];
  TH2D *hist_yzpattern2D[100];
  TH2D *hist_xzpattern2D2[100];
  TH2D *hist_yzpattern2D2[100];

};


double hodo_3body::GetCos(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Px()*lv2.Px() + lv1.Py()*lv2.Py() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / (lv1.P() * lv2.P());
  return cos;
}

double hodo_3body::GetCos_xz(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Px()*lv2.Px() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / TMath::Sqrt((lv1.Px()*lv1.Px()+lv1.Pz()*lv1.Pz()) * (lv2.Px()*lv2.Px()+lv2.Pz()*lv2.Pz()));
  return cos;
}

/*
  double hodo_3body::GetCos_yz(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Py()*lv2.Py() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / TMath::Sqrt((lv1.Py()*lv1.Py()+lv1.Pz()*lv1.Pz()) * (lv2.Py()*lv2.Py()+lv2.Pz()*lv2.Pz()));
  return cos;
  }

*/
double hodo_3body::GetCos_yz(TLorentzVector lv1) {
  double dot_product = lv1.Py();
  double cos = dot_product / lv1.P();
  return cos;
}
double hodo_3body::GetAng_yz(TLorentzVector lv1){
  //  double rad = 57.2958;
  double ang = 90-(180*TMath::ACos(GetCos_yz(lv1)))/TMath::Pi();
  return ang;
}
double hodo_3body::GetAng(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(GetCos(lv1,lv2)))/TMath::Pi();
  return ang;
}

double hodo_3body::GetAng_xz(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(GetCos_xz(lv1,lv2)))/TMath::Pi();
  return ang;
}
/*
  double hodo_3body::GetAng_yz(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(GetCos_yz(lv1,lv2)))/TMath::Pi();
  return ang;
  }
*/

void hodo_3body::pipip(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

  TTree *tree = (TTree*)file->Get("tree");

  //bool pdf=false;
  //bool pdf=true;

  gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);

  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("toftime",toftime);
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
  //tree->SetBranchAddress("targetpid",targetpid);

  can_multi[class_hodo] = new TCanvas(Form("can_multi_%d",class_hodo),"",1200,800);
  can_deltatof[class_hodo] = new TCanvas(Form("can_deltatof_%d",class_hodo),"",800,800);
  can_segpattern2D[class_hodo] = new TCanvas(Form("can_segpattern2D_%d",class_hodo),"",800,800);
  can_xzpattern2D[class_hodo] = new TCanvas(Form("can_xzpattern2D_%d",class_hodo),"",800,800);
  can_yzpattern2D[class_hodo] = new TCanvas(Form("can_yzpattern2D_%d",class_hodo),"",800,800);
  can_xzpattern2D2[class_hodo] = new TCanvas(Form("can_xzpattern2D2_%d",class_hodo),"",800,800);
  can_yzpattern2D2[class_hodo] = new TCanvas(Form("can_yzpattern2D2_%d",class_hodo),"",800,800);
  can_hitpattern2D[class_hodo] = new TCanvas(Form("can_hitpattern2D_%d",class_hodo),"",800,800);
  can_hitpattern2D[class_hodo] -> Divide(1,2);

  //can_pid[class_hodo] = new TCanvas(Form("can_pid_%d",class_hodo),"",800,800);
  //hist_pid[class_hodo] =new TH1D(Form("hist_pid_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);

  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_deltatof[class_hodo] =new TH1D(Form("hist_deltatof_%d",class_hodo),"Hodoscope #Delta tof",100,0,5);
  if(window==true) hist_segpattern2D[class_hodo] = new TH2D(Form("hist_segpattern2D_%d",class_hodo),"Hodoscope hitpattern;ID of #pi hit for TPC-Hodo;ID of p hit for TPC-Hodo;",38,0,38,38,0,38);
  if(window==false) hist_segpattern2D[class_hodo] = new TH2D(Form("hist_segpattern2D_%d",class_hodo),"Hodoscope hitpattern;ID of #pi hit for TPC-Hodo;ID of p hit for TPC-Hodo;",32,0,32,32,0,32);

  hist_xzpattern2D[class_hodo] = new TH2D(Form("hist_xzpattern2D_%d",class_hodo),"#theta_{xz} pattern;#theta_{xz} of #pi track for TPC-Hodo;#theta_{xz} of p track for TPC-Hodo",100,-180,180,100,-180,180);
  hist_yzpattern2D[class_hodo] = new TH2D(Form("hist_yzpattern2D_%d",class_hodo),"#theta_{yz} pattern;#theta_{yz} of #pi track for TPC-Hodo;#theta_{yz} of p track for TPC-Hodo",100,-90,90,100,-90,90);
  hist_xzpattern2D2[class_hodo] = new TH2D(Form("hist_xzpattern2D2_%d",class_hodo),"single hit, #theta_{xz} pattern;#theta_{xz} of #pi track for TPC-Hodo;#theta_{xz} of p track for TPC-Hodo",100,-180,180,100,-180,180);
  hist_yzpattern2D2[class_hodo] = new TH2D(Form("hist_yzpattern2D2_%d",class_hodo),"single hit, #theta_{yz} pattern;#theta_{yz} of #pi track for TPC-Hodo;#theta_{yz} of p track for TPC-Hodo",100,-90,90,100,-90,90);

  hist_hitpattern2D_pi[class_hodo] = new TH2D(Form("hist_hitpattern2D_pi_%d",class_hodo),"#pi Hit pattern near the window;X(mm);Y(mm)",1000,-140,140,1000,-50,50);
  hist_hitpattern2D_p[class_hodo] = new TH2D(Form("hist_hitpattern2D_p_%d",class_hodo),"p Hit pattern near the window;X(mm);Y(mm)",1000,-140,140,1000,-50,50);

  int nevent=tree->GetEntries();
  int count_p,count_pi,count_mu,count_e; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi, flag_seg_p;

  int hodo_seg[38]={0};
  int hodo_seg_p[38]={0};
  int hodo_seg_pi[38]={0};

  double flag_time_pi, flag_time_p;
  double hodo_hittime_p[38]={0};
  double hodo_hittime_pi[38]={0};

  double m_p = 0.938272;
  double m_pi = 0.139570;
  double p_beam = 2.0;
  double e_beam = sqrt(p_beam*p_beam + m_pi*m_pi);
  TLorentzVector lv_beam(0,0,p_beam,e_beam); //p1
  double theta_xz_p,theta_xz_pi;
  double theta_yz_p,theta_yz_pi;
  double E_p,E_pi;

  int count_m0=0;
  int count_m1=0, count_m1_coin=0, count_m1_p=0, count_m1_pi=0, count_m1_dummy=0;
  int count_m2=0, count_m2_coin=0, count_m2_multi=0, count_m2_p=0, count_m2_pi=0, count_m2_dummy=0;
  int dummy=0;

  for(int i=0;i<nevent;i++){
    tree->GetEntry(i);

    count_p=0;
    count_pi=0;
    count_e=0;
    count_mu=0;
    multiplicity=0;
    check_flag=0;
    flag_seg_pi, flag_seg_p;
    flag_time_pi=9999.;
    flag_time_p=9999.;
    theta_xz_p=9999;
    theta_xz_pi=9999;
    theta_yz_p=9999;
    theta_yz_pi=9999;

    for(int j=0;j<nEvt;j++){
      if(TMath::Abs(evtpid[j])==211){
	E_pi = TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]+m_pi*m_pi);
	TLorentzVector lv_pi(evtpx[j],evtpy[j],evtpz[j],E_pi);
	theta_xz_pi = (evtpx[j]/TMath::Abs(evtpx[j]))*GetAng_xz(lv_beam,lv_pi);
	//theta_yz_pi = (evtpy[j]/TMath::Abs(evtpy[j]))*GetAng_yz(lv_beam,lv_pi);
	theta_yz_pi = GetAng_yz(lv_pi);
      }
      if(evtpid[j]==2212){
	E_p = TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]+m_p*m_p);
	TLorentzVector lv_p(evtpx[j],evtpy[j],evtpz[j],E_p);

	theta_xz_p =(evtpx[j]/TMath::Abs(evtpx[j]))*GetAng_xz(lv_beam,lv_p);
	//theta_yz_p =(evtpy[j]/TMath::Abs(evtpy[j]))*GetAng_yz(lv_beam,lv_p);
	theta_yz_p =GetAng_yz(lv_p);
      }
    }

    for(int j=0;j<nhTof;j++){
      //hist_pid[class_hodo] -> Fill(tofpid[j]); //PID
      if(window==true&&tofseg[j]==0&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==1&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==2&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==3&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==17&&TMath::Abs(tofposy[j])<window_size) dummy++;
      else if(window==true&&tofseg[j]==18&&TMath::Abs(tofposy[j])<window_size) dummy++;

      else if(TMath::Abs(tofpid[j])==211&&tofedep[j]>ecut){ //for selecting events w/ pi
	count_pi++;
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  hodo_seg_pi[32+0] += 1;
	  hodo_hittime_pi[32+0]=toftime[j];
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  hodo_seg_pi[32+1] += 1;
	  hodo_hittime_pi[32+1]=toftime[j];
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 1;
	  hodo_seg_pi[32+2] += 1;
	  hodo_hittime_pi[32+2]=toftime[j];
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 1;
	  hodo_seg_pi[32+3] += 1;
	  hodo_hittime_pi[32+3]=toftime[j];
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 1;
	  hodo_seg_pi[32+4] += 1;
	  hodo_hittime_pi[32+4]=toftime[j];
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 1;
	  hodo_seg_pi[32+5] += 1;
	  hodo_hittime_pi[32+5]=toftime[j];
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  hodo_seg_pi[tofseg[j]] += 1;
	  hodo_hittime_pi[tofseg[j]]=toftime[j];
	}
      }
      else if(tofpid[j]==2212&&tofedep[j]>ecut){ //for selecting events w/ p
	count_p++;
	if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
	  hodo_seg[32+0] += 1;
	  hodo_seg_p[32+0] += 1;
	  hodo_hittime_p[32+0]=toftime[j];
	}
	else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
	  hodo_seg[32+1] += 1;
	  hodo_seg_p[32+1] += 1;
	  hodo_hittime_p[32+1]=toftime[j];
	}
	else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
	  hodo_seg[32+2] += 1;
	  hodo_seg_p[32+2] += 1;
	  hodo_hittime_p[32+2]=toftime[j];
	}
	else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
	  hodo_seg[32+3] += 1;
	  hodo_seg_p[32+3] += 1;
	  hodo_hittime_p[32+3]=toftime[j];
	}
	else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
	  hodo_seg[32+4] += 1;
	  hodo_seg_p[32+4] += 1;
	  hodo_hittime_p[32+4]=toftime[j];
	}
	else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
	  hodo_seg[32+5] += 1;
	  hodo_seg_p[32+5] += 1;
	  hodo_hittime_p[32+5]=toftime[j];
	}
	else{
	  hodo_seg[tofseg[j]] += 1;
	  hodo_seg_p[tofseg[j]] += 1;
	  hodo_hittime_p[tofseg[j]]=toftime[j];
	}
      }
      else if(TMath::Abs(tofpid[j])==11&&tofedep[j]>ecut){
	count_e++;
	hodo_seg[tofseg[j]] += 1;
      }
      else if(TMath::Abs(tofpid[j])==13&&tofedep[j]>ecut){
	count_mu++;
	hodo_seg[tofseg[j]] += 1;
      }
      //std::cout<<"event: "<<i<<std::endl;
    }

    if(count_p>0){
      for(int i=0;i<38;i++){
	if(hodo_seg_p[i]>0&&flag_time_p>hodo_hittime_p[i]){
	  flag_time_p=hodo_hittime_p[i];
	  flag_seg_p=i;
	  //std::cout<<"Proton hodo time: "<<flag_time_p<<" hodo seg: "<<flag_seg_p<<std::endl;
	}
      }
    }
    if(count_pi>0){
      for(int i=0;i<38;i++){
	if(hodo_seg_pi[i]>0&&flag_time_pi>hodo_hittime_pi[i]){
	  flag_time_pi=hodo_hittime_pi[i];
	  flag_seg_pi=i;
	  //std::cout<<"Pion hodo time: "<<flag_time_pi<<" hodo seg: "<<flag_seg_pi<<std::endl;
	}
      }
    }
    /*
      if(count_p>0 && count_pi>0){
      hist_deltatof[class_hodo] -> Fill(flag_time_p-flag_time_pi);
      if(flag_time_p < flag_time_pi) hist_segpattern2D[class_hodo] -> Fill(flag_seg_p, flag_seg_pi);
      else hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi, flag_seg_p);
      }
    */
    if(count_p>0 && count_pi==0){
      hist_xzpattern2D2[class_hodo] -> Fill(theta_xz_pi,theta_xz_p);
      hist_yzpattern2D2[class_hodo] -> Fill(theta_yz_pi,theta_yz_p);
    }
    if(count_pi>0 && count_p==0){
      hist_xzpattern2D2[class_hodo] -> Fill(theta_xz_pi,theta_xz_p);
      hist_yzpattern2D2[class_hodo] -> Fill(theta_yz_pi,theta_yz_p);
    }
    if(count_p>0 && count_pi>0){
      if(flag_seg_pi<32) flag_seg_pi=(flag_seg_pi+24)%32;
      if(flag_seg_p<32) flag_seg_p=(flag_seg_p+24)%32;

      hist_deltatof[class_hodo] -> Fill(flag_time_p-flag_time_pi);
      hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi, flag_seg_p);
      hist_xzpattern2D[class_hodo] -> Fill(theta_xz_pi,theta_xz_p);
      hist_yzpattern2D[class_hodo] -> Fill(theta_yz_pi,theta_yz_p);
    }
    for(int i=0;i<38;i++){
      if(hodo_seg[i]>0) multiplicity++;
      check_flag += hodo_seg[i];
      hodo_seg[i]=0;
      hodo_seg_p[i]=0;
      hodo_seg_pi[i]=0;
      hodo_hittime_p[i]=0;
      hodo_hittime_pi[i]=0;
    }

    //counting check
    if(check_flag!=count_p+count_pi+count_e+count_mu) std::cout<< "fatal error !!! "<<std::endl;

    //M-1
    if(multiplicity==1){
      if(count_p==1 && count_pi==1 && flag_seg_pi==flag_seg_p) count_m1_coin++;
      if(count_p==1 && count_pi==0) count_m1_p++;
      if(count_pi==1 && count_p==0) count_m1_pi++;
      if(count_pi==0 && count_p==0) count_m1_dummy++;
      count_m1++;
    }

    //M-2
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

    //M-0
    else count_m0++;
    hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
  }

  //can_pid[class_hodo] -> cd();
  //hist_pid[class_hodo]->Draw();
  //hist_pid[class_hodo]->GetXaxis()->SetTitle("PID(PDG encording) ");
  //hist_pid[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_multi[class_hodo]->cd();
  hist_multi[class_hodo]->Draw();
  hist_multi[class_hodo]->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_segpattern2D[class_hodo]->cd();
  hist_segpattern2D[class_hodo]->Draw("colz");
  hist_segpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,250);

  can_xzpattern2D[class_hodo]->cd();
  hist_xzpattern2D[class_hodo]->Draw("colz");
  hist_xzpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_yzpattern2D[class_hodo]->cd();
  hist_yzpattern2D[class_hodo]->Draw("colz");
  hist_yzpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_xzpattern2D2[class_hodo]->cd();
  hist_xzpattern2D2[class_hodo]->Draw("colz");
  hist_xzpattern2D2[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_yzpattern2D2[class_hodo] -> cd();
  hist_yzpattern2D2[class_hodo] -> Draw("colz");
  hist_yzpattern2D2[class_hodo] -> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_hitpattern2D[class_hodo]->cd(1);
  hist_hitpattern2D_pi[class_hodo] -> Draw("colz");
  hist_hitpattern2D_pi[class_hodo] -> GetZaxis() -> SetRangeUser(0,100);

  can_hitpattern2D[class_hodo]->cd(2);
  hist_hitpattern2D_p[class_hodo] -> Draw("colz");
  hist_hitpattern2D_p[class_hodo] -> GetZaxis() -> SetRangeUser(0,100);

  can_deltatof[class_hodo]->cd();
  hist_deltatof[class_hodo]->Draw();

  if(printout==false){
    can_multi[class_hodo] -> Close();
  }

  acceptance= (double) count_m2_multi/nevent;

  if(printout){
    std::cout<<""<<std::endl;
    std::cout<<"Event number   : "<<nevent<<std::endl;
    if(nevent!=count_m0+count_m1+count_m2){
      std::cout<<"event number error !!! "<<std::endl;
      std::cout<<"M-0,1,2 number : "<<count_m0+count_m1+count_m2<<std::endl;
    }
    //std::cout<<"hit window     : "<<dummy<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-0 number   : "<<(double) count_m0<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"multi-1 number   : "<<count_m1<<std::endl;
    std::cout<<"multi-1 coin     : "<<count_m1_coin<<std::endl;
    std::cout<<"multi-1 P hit    : "<<count_m1_p<<std::endl;
    std::cout<<"multi-1 pi hit   : "<<count_m1_pi<<std::endl;
    std::cout<<"multi-1 e,mu hit : "<<count_m1_dummy<<std::endl;

    if(count_m1!=count_m1_coin+count_m1_p+count_m1_pi+count_m1_dummy){
      std::cout<< "M-1 error !!! "<<std::endl;
      std::cout<<"multi-1 check    : "<<count_m1_coin+count_m1_p+count_m1_pi+count_m1_dummy<<std::endl;
    }

    std::cout<<""<<std::endl;
    std::cout<<"multi-2 number   : "<<count_m2<<std::endl;
    std::cout<<"multi-2 coin     : "<<count_m2_coin<<std::endl;
    std::cout<<"multi-2 Ppi hits : "<<count_m2_multi<<std::endl;
    std::cout<<"multi-2 P only   : "<<count_m2_p<<std::endl;
    std::cout<<"multi-2 pi only  : "<<count_m2_pi<<std::endl;
    std::cout<<"multi-2 e,mu hits: "<<count_m2_dummy<<std::endl;
    if(count_m2!=count_m2_multi+count_m2_p+count_m2_pi+count_m2_dummy){
      std::cout<<"M-2 error !!! "<<std::endl;
      std::cout<<"multi-2 check    : "<<count_m2_multi+count_m2_p+count_m2_pi+count_m2_dummy<<std::endl;
    }

    std::cout<<""<<std::endl;
    std::cout<<"acceptance     : "<<(double) acceptance*100<< " %"<<std::endl;
    std::cout<<"multi-0 ratio    : "<<(double) count_m0/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-1 ratio    : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
    std::cout<<"multi-2 ratio    : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

    std::cout<<"m1 coin/m1       : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
    std::cout<<"m2 coin/m2       : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

    std::cout<<""<<std::endl;
    std::cout<<"m-0,1,2 ratio    : "<<(double) (count_m0+count_m1+count_m2)/nevent*100<< " %"<<std::endl;
  }

  class_hodo++;
}



void hodo_3body::pipin(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

  TTree *tree = (TTree*)file->Get("tree");

  //bool pdf=false;
  //bool pdf=true;

  gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);

  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("toftime",toftime);
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
  //tree->SetBranchAddress("targetpid",targetpid);

  can_multi[class_hodo] = new TCanvas(Form("can_multi_%d",class_hodo),"",800,800);
  can_deltatof[class_hodo] = new TCanvas(Form("can_deltatof_%d",class_hodo),"",800,800);
  can_segpattern2D[class_hodo] = new TCanvas(Form("can_segpattern2D_%d",class_hodo),"",800,800);
  can_xzpattern2D[class_hodo] = new TCanvas(Form("can_xzpattern2D_%d",class_hodo),"",800,800);
  can_yzpattern2D[class_hodo] = new TCanvas(Form("can_yzpattern2D_%d",class_hodo),"",800,800);
  can_xzpattern2D2[class_hodo] = new TCanvas(Form("can_xzpattern2D2_%d",class_hodo),"",800,800);
  can_yzpattern2D2[class_hodo] = new TCanvas(Form("can_yzpattern2D2_%d",class_hodo),"",800,800);

  //can_pi_pd[class_hodo] = new TCanvas(Form("can_pi_pd_%d",class_hodo),"",800,800);
  //hist_pi_pd[class_hodo] =new TH1D(Form("hist_pi_pd_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);
  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_deltatof[class_hodo] =new TH1D(Form("hist_deltatof_%d",class_hodo),"Hodoscope #Delta tof",100,-5,5);

  if(window==true) hist_segpattern2D[class_hodo] = new TH2D(Form("hist_segpattern2D_%d",class_hodo),"Hodoscope hitpattern;ID of #pi^{+} hit for TPC-Hodo;ID of #pi^{-} hit for TPC-Hodo;",38,0,38,38,0,38);
  if(window==false) hist_segpattern2D[class_hodo] = new TH2D(Form("hist_segpattern2D_%d",class_hodo),"Hodoscope hitpattern;ID of #pi^{+} hit for TPC-Hodo;ID of #pi^{-} hit for TPC-Hodo;",32,0,32,32,0,32);

  hist_xzpattern2D[class_hodo] = new TH2D(Form("hist_xzpattern2D_%d",class_hodo),"#theta_{xz} pattern;#theta_{xz} of #pi^{+} track for TPC-Hodo;#theta_{xz} of #pi^{-} track for TPC-Hodo",100,-180,180,100,-180,180);
  hist_yzpattern2D[class_hodo] = new TH2D(Form("hist_yzpattern2D_%d",class_hodo),"#theta_{yz} pattern;#theta_{yz} of #pi^{+} track for TPC-Hodo;#theta_{yz} of #pi^{-} track for TPC-Hodo",100,-90,90,100,-90,90);
  hist_xzpattern2D2[class_hodo] = new TH2D(Form("hist_xzpattern2D2_%d",class_hodo),"single hit, #theta_{xz} pattern;#theta_{xz} of #pi^{+} track for TPC-Hodo;#theta_{xz} of #pi^{-} track for TPC-Hodo",100,-180,180,100,-180,180);
  hist_yzpattern2D2[class_hodo] = new TH2D(Form("hist_yzpattern2D2_%d",class_hodo),"single hit, #theta_{yz} pattern;#theta_{yz} of #pi^{+} track for TPC-Hodo;#theta_{yz} of #pi^{-} track for TPC-Hodo",100,-90,90,100,-90,90);

  int nevent=tree->GetEntries();
  int count_pi_m,count_pi_p,count_mu,count_e; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi_p, flag_seg_pi_m;

  int hodo_seg[38]={0};
  int hodo_seg_pi_m[38]={0};
  int hodo_seg_pi_p[38]={0};

  double flag_time_pi_p, flag_time_pi_m;
  double hodo_hittime_pi_m[38]={0};
  double hodo_hittime_pi_p[38]={0};

  double m_pi = 0.139570;
  double p_beam = 2.0;
  double e_beam = sqrt(p_beam*p_beam + m_pi*m_pi);
  TLorentzVector lv_beam(0,0,p_beam,e_beam); //p1
  double theta_xz_pi_p,theta_xz_pi_m;
  double theta_yz_pi_p,theta_yz_pi_m;
  double E_pi_p,E_pi_m;

  int count_m0=0;
  int count_m1=0, count_m1_coin=0, count_m1_pi_m=0, count_m1_pi_p=0, count_m1_dummy=0;
  int count_m2=0, count_m2_coin=0, count_m2_multi=0, count_m2_pi_m=0, count_m2_pi_p=0, count_m2_dummy=0;
  int dummy=0;

  for(int i=0;i<nevent;i++){
  tree->GetEntry(i);

  count_pi_m=0;
  count_pi_p=0;
  count_e=0;
  count_mu=0;
  multiplicity=0;
  check_flag=0;
  flag_seg_pi_p, flag_seg_pi_m;
  flag_time_pi_p=9999.;
  flag_time_pi_m=9999.;
  theta_xz_pi_p=9999;
  theta_xz_pi_m=9999;
  theta_yz_pi_p=9999;
  theta_yz_pi_m=9999;

  for(int j=0;j<nEvt;j++){
  if(evtpid[j]==211){
  E_pi_p = TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]+m_pi*m_pi);
  TLorentzVector lv_pi_p(evtpx[j],evtpy[j],evtpz[j],E_pi_p);
  theta_xz_pi_p = (evtpx[j]/TMath::Abs(evtpx[j]))*GetAng_xz(lv_beam,lv_pi_p);
  //theta_yz_pi_p = (evtpy[j]/TMath::Abs(evtpy[j]))*GetAng_yz(lv_beam,lv_pi_p);
  theta_yz_pi_p = GetAng_yz(lv_pi_p);
}
  if(evtpid[j]==-211){
  E_pi_m = TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]+m_pi*m_pi);
  TLorentzVector lv_pi_m(evtpx[j],evtpy[j],evtpz[j],E_pi_m);

  theta_xz_pi_m =(evtpx[j]/TMath::Abs(evtpx[j]))*GetAng_xz(lv_beam,lv_pi_m);
  //theta_yz_pi_m =(evtpy[j]/TMath::Abs(evtpy[j]))*GetAng_yz(lv_beam,lv_pi_m);
  theta_yz_pi_m =GetAng_yz(lv_pi_m);
}
}

  for(int j=0;j<nhTof;j++){
  //hist_pi_pd[class_hodo] -> Fill(tofpid[j]); //PID
  if(window==true&&tofseg[j]==0&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==1&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==2&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==3&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==17&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==18&&TMath::Abs(tofposy[j])<window_size) dummy++;

  else if(tofpid[j]==211&&tofedep[j]>ecut){ //for selecting events w/ pi
  count_pi_p++;
  if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
  hodo_seg[32+0] += 1;
  hodo_seg_pi_p[32+0] += 1;
  hodo_hittime_pi_p[32+0]=toftime[j];
}
  else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
  hodo_seg[32+1] += 1;
  hodo_seg_pi_p[32+1] += 1;
  hodo_hittime_pi_p[32+1]=toftime[j];
}
  else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
  hodo_seg[32+2] += 1;
  hodo_seg_pi_p[32+2] += 1;
  hodo_hittime_pi_p[32+2]=toftime[j];
}
  else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
  hodo_seg[32+3] += 1;
  hodo_seg_pi_p[32+3] += 1;
  hodo_hittime_pi_p[32+3]=toftime[j];
}
  else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
  hodo_seg[32+4] += 1;
  hodo_seg_pi_p[32+4] += 1;
  hodo_hittime_pi_p[32+4]=toftime[j];
}
  else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
  hodo_seg[32+5] += 1;
  hodo_seg_pi_p[32+5] += 1;
  hodo_hittime_pi_p[32+5]=toftime[j];
}
  else{
  hodo_seg[tofseg[j]] += 1;
  hodo_seg_pi_p[tofseg[j]] += 1;
  hodo_hittime_pi_p[tofseg[j]]=toftime[j];
}
}
  else if(tofpid[j]==-211&&tofedep[j]>ecut){ //for selecting events w/ p
  count_pi_m++;
  if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
  hodo_seg[32+0] += 1;
  hodo_seg_pi_m[32+0] += 1;
  hodo_hittime_pi_m[32+0]=toftime[j];
}
  else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
  hodo_seg[32+1] += 1;
  hodo_seg_pi_m[32+1] += 1;
  hodo_hittime_pi_m[32+1]=toftime[j];
}
  else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
  hodo_seg[32+2] += 1;
  hodo_seg_pi_m[32+2] += 1;
  hodo_hittime_pi_m[32+2]=toftime[j];
}
  else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
  hodo_seg[32+3] += 1;
  hodo_seg_pi_m[32+3] += 1;
  hodo_hittime_pi_m[32+3]=toftime[j];
}
  else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
  hodo_seg[32+4] += 1;
  hodo_seg_pi_m[32+4] += 1;
  hodo_hittime_pi_m[32+4]=toftime[j];
}
  else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
  hodo_seg[32+5] += 1;
  hodo_seg_pi_m[32+5] += 1;
  hodo_hittime_pi_m[32+5]=toftime[j];
}
  else{
  hodo_seg[tofseg[j]] += 1;
  hodo_seg_pi_m[tofseg[j]] += 1;
  hodo_hittime_pi_m[tofseg[j]]=toftime[j];
}
}
  else if(TMath::Abs(tofpid[j])==11&&tofedep[j]>ecut){
  count_e++;
  hodo_seg[tofseg[j]] += 1;
}
  else if(TMath::Abs(tofpid[j])==13&&tofedep[j]>ecut){
  count_mu++;
  hodo_seg[tofseg[j]] += 1;
}
  //std::cout<<"event: "<<i<<std::endl;
}

  if(count_pi_m>0){
  for(int i=0;i<38;i++){
  if(hodo_seg_pi_m[i]>0&&flag_time_pi_m>hodo_hittime_pi_m[i]){
  flag_time_pi_m=hodo_hittime_pi_m[i];
  flag_seg_pi_m=i;
  //std::cout<<"Proton hodo time: "<<flag_time_pi_m<<" hodo seg: "<<flag_seg_pi_m<<std::endl;
}
}
}
  if(count_pi_p>0){
  for(int i=0;i<38;i++){
  if(hodo_seg_pi_p[i]>0&&flag_time_pi_p>hodo_hittime_pi_p[i]){
  flag_time_pi_p=hodo_hittime_pi_p[i];
  flag_seg_pi_p=i;
  //std::cout<<"Pion hodo time: "<<flag_time_pi_p<<" hodo seg: "<<flag_seg_pi_p<<std::endl;
}
}
}
  /*
    if(count_pi_m>0 && count_pi_p>0){
    hist_deltatof[class_hodo] -> Fill(flag_time_pi_m-flag_time_pi_p);
    if(flag_time_pi_m < flag_time_pi_p) hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi_m, flag_seg_pi_p);
    else hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi_p, flag_seg_pi_m);
    }
  */
  if(count_pi_p>0 && count_pi_m==0){
  hist_xzpattern2D2[class_hodo] -> Fill(theta_xz_pi_p,theta_xz_pi_m);
  hist_yzpattern2D2[class_hodo] -> Fill(theta_yz_pi_p,theta_yz_pi_m);
}
  if(count_pi_p>0 && count_pi_m==0){
  hist_xzpattern2D2[class_hodo] -> Fill(theta_xz_pi_p,theta_xz_pi_m);
  hist_yzpattern2D2[class_hodo] -> Fill(theta_yz_pi_p,theta_yz_pi_m);

}
  if(count_pi_p>0 && count_pi_m>0){
  if(flag_seg_pi_p<32) flag_seg_pi_p=(flag_seg_pi_p+24)%32;
  if(flag_seg_pi_m<32) flag_seg_pi_m=(flag_seg_pi_m+24)%32;

  hist_deltatof[class_hodo] -> Fill(flag_time_pi_p-flag_time_pi_m);
  hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi_p, flag_seg_pi_m);
  hist_xzpattern2D[class_hodo] -> Fill(theta_xz_pi_p,theta_xz_pi_m);
  hist_yzpattern2D[class_hodo] -> Fill(theta_yz_pi_p,theta_yz_pi_m);
}

  for(int i=0;i<38;i++){
  if(hodo_seg[i]>0) multiplicity++;
  check_flag += hodo_seg[i];
  hodo_seg[i]=0;
  hodo_seg_pi_m[i]=0;
  hodo_seg_pi_p[i]=0;
  hodo_hittime_pi_m[i]=0;
  hodo_hittime_pi_p[i]=0;
}

  //counting check
  if(check_flag!=count_pi_m+count_pi_p+count_e+count_mu) std::cout<< "fatal error !!! "<<std::endl;

  //M-1
  if(multiplicity==1){
  if(count_pi_m==1 && count_pi_p==1 && flag_seg_pi_p==flag_seg_pi_m) count_m1_coin++;
  if(count_pi_m==1 && count_pi_p==0) count_m1_pi_m++;
  if(count_pi_p==1 && count_pi_m==0) count_m1_pi_p++;
  if(count_pi_p==0 && count_pi_m==0) count_m1_dummy++;
  count_m1++;
}

  //M-2
  else if(multiplicity>1){
  if(count_pi_m==1 && count_pi_p==1 && flag_seg_pi_p==flag_seg_pi_m) count_m2_coin++;
  if(count_pi_m>0 && count_pi_p>0) count_m2_multi++;
  if(count_pi_m>0 && count_pi_p==0){
  count_m2_pi_m++;
}
  if(count_pi_p>0 && count_pi_m==0){
  count_m2_pi_p++;
}
  if(count_pi_p==0 && count_pi_m==0) count_m2_dummy++;
  count_m2++;
}

  //M-0
  else count_m0++;
  hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
}

  //can_pi_pd[class_hodo] -> cd();
  //hist_pi_pd[class_hodo]->Draw();
  //hist_pi_pd[class_hodo]->GetXaxis()->SetTitle("PID(PDG encording) ");
  //hist_pi_pd[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_multi[class_hodo]->cd();
  hist_multi[class_hodo]->Draw();
  hist_multi[class_hodo]->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_segpattern2D[class_hodo]->cd();
  hist_segpattern2D[class_hodo]->Draw("colz");
  hist_segpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,250);

  can_xzpattern2D[class_hodo]->cd();
  hist_xzpattern2D[class_hodo]->Draw("colz");
  hist_xzpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_yzpattern2D[class_hodo]->cd();
  hist_yzpattern2D[class_hodo]->Draw("colz");
  hist_yzpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_xzpattern2D2[class_hodo]->cd();
  hist_xzpattern2D2[class_hodo]->Draw("colz");
  hist_xzpattern2D2[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_yzpattern2D2[class_hodo] -> cd();
  hist_yzpattern2D2[class_hodo] -> Draw("colz");
  hist_yzpattern2D2[class_hodo] -> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();

  can_deltatof[class_hodo]->cd();
  hist_deltatof[class_hodo]->Draw();

  if(printout==false){
  can_multi[class_hodo] -> Close();
}

  acceptance= (double) count_m2_multi/nevent;

  if(printout){
  std::cout<<""<<std::endl;
  std::cout<<"Event number   : "<<nevent<<std::endl;
  if(nevent!=count_m0+count_m1+count_m2){
  std::cout<<"event number error !!! "<<std::endl;
  std::cout<<"M-0,1,2 number : "<<count_m0+count_m1+count_m2<<std::endl;
}
  //std::cout<<"hit window     : "<<dummy<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"multi-0 number   : "<<(double) count_m0<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"acceptance     : "<<(double) acceptance*100<< " %"<<std::endl;
  std::cout<<"multi-1 number   : "<<count_m1<<std::endl;
  std::cout<<"multi-1 coin     : "<<count_m1_coin<<std::endl;
  std::cout<<"multi-1 Pi- hit    : "<<count_m1_pi_m<<std::endl;
  std::cout<<"multi-1 pi+ hit   : "<<count_m1_pi_p<<std::endl;
  std::cout<<"multi-1 e,mu hit : "<<count_m1_dummy<<std::endl;

  if(count_m1!=count_m1_coin+count_m1_pi_m+count_m1_pi_p+count_m1_dummy){
  std::cout<< "M-1 error !!! "<<std::endl;
  std::cout<<"multi-1 check    : "<<count_m1_coin+count_m1_pi_m+count_m1_pi_p+count_m1_dummy<<std::endl;
}

  std::cout<<""<<std::endl;
  std::cout<<"multi-2 number   : "<<count_m2<<std::endl;
  std::cout<<"multi-2 coin     : "<<count_m2_coin<<std::endl;
  std::cout<<"multi-2 Ppi hits : "<<count_m2_multi<<std::endl;
  std::cout<<"multi-2 pi- only   : "<<count_m2_pi_m<<std::endl;
  std::cout<<"multi-2 pi+ only  : "<<count_m2_pi_p<<std::endl;
  std::cout<<"multi-2 e,mu hits: "<<count_m2_dummy<<std::endl;
  if(count_m2!=count_m2_multi+count_m2_pi_m+count_m2_pi_p+count_m2_dummy){
  std::cout<<"M-2 error !!! "<<std::endl;
  std::cout<<"multi-2 check    : "<<count_m2_multi+count_m2_pi_m+count_m2_pi_p+count_m2_dummy<<std::endl;
}

  std::cout<<""<<std::endl;
  std::cout<<"multi-0 ratio    : "<<(double) count_m0/nevent*100<< " %"<<std::endl;
  std::cout<<"multi-1 ratio    : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
  std::cout<<"multi-2 ratio    : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

  std::cout<<"m1 coin/m1       : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
  std::cout<<"m2 coin/m2       : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"m-0,1,2 ratio    : "<<(double) (count_m0+count_m1+count_m2)/nevent*100<< " %"<<std::endl;
}

  class_hodo++;
}






  void hodo_3body::pipin2(TFile *file, double ecut, bool window, double &acceptance, bool printout=false){

  TTree *tree = (TTree*)file->Get("tree");

  //bool pdf=false;
  //bool pdf=true;

  gStyle->SetOptStat(0);

  tree->SetBranchAddress("event",&event);

  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("toftime",toftime);
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
  //tree->SetBranchAddress("targetpid",targetpid);

  can_multi[class_hodo] = new TCanvas(Form("can_multi_%d",class_hodo),"",800,800);
  can_deltatof[class_hodo] = new TCanvas(Form("can_deltatof_%d",class_hodo),"",800,800);
  can_segpattern2D[class_hodo] = new TCanvas(Form("can_segpattern2D_%d",class_hodo),"",800,800);
  can_xzpattern2D[class_hodo] = new TCanvas(Form("can_xzpattern2D_%d",class_hodo),"",800,800);
  can_yzpattern2D[class_hodo] = new TCanvas(Form("can_yzpattern2D_%d",class_hodo),"",800,800);
  can_xzpattern2D2[class_hodo] = new TCanvas(Form("can_xzpattern2D2_%d",class_hodo),"",800,800);
  can_yzpattern2D2[class_hodo] = new TCanvas(Form("can_yzpattern2D2_%d",class_hodo),"",800,800);

  //can_pi_pd[class_hodo] = new TCanvas(Form("can_pi_pd_%d",class_hodo),"",800,800);
  //hist_pi_pd[class_hodo] =new TH1D(Form("hist_pi_pd_%d",class_hodo),"Hodoscope PID",4600,-2300,2300);
  hist_multi[class_hodo] =new TH1D(Form("hist_multi_%d",class_hodo),"Hodoscope Multiplicity",7,0,7);
  hist_deltatof[class_hodo] =new TH1D(Form("hist_deltatof_%d",class_hodo),"Hodoscope #Delta tof",100,0,5);
  if(window==true) hist_segpattern2D[class_hodo] = new TH2D(Form("hist_segpattern2D_%d",class_hodo),"Hodoscope hitpattern;ID of #pi^{+} hit for TPC-Hodo;ID of #pi^{+} hit for TPC-Hodo;",38,0,38,38,0,38);
  if(window==false) hist_segpattern2D[class_hodo] = new TH2D(Form("hist_segpattern2D_%d",class_hodo),"Hodoscope hitpattern;ID of #pi^{+} hit for TPC-Hodo;ID of #pi^{+} hit for TPC-Hodo;",32,0,32,32,0,32);

  hist_xzpattern2D[class_hodo] = new TH2D(Form("hist_xzpattern2D_%d",class_hodo),"#theta_{xz} pattern;#theta_{xz} of #pi^{+} track for TPC-Hodo;#theta_{xz} of #pi^{+} track for TPC-Hodo",100,-180,180,100,-180,180);
  hist_yzpattern2D[class_hodo] = new TH2D(Form("hist_yzpattern2D_%d",class_hodo),"#theta_{yz} pattern;#theta_{yz} of #pi^{+} track for TPC-Hodo;#theta_{yz} of #pi^{+} track for TPC-Hodo",100,-90,90,100,-90,90);
  hist_xzpattern2D2[class_hodo] = new TH2D(Form("hist_xzpattern2D2_%d",class_hodo),"single hit, #theta_{xz} pattern;#theta_{xz} of #pi^{+} track for TPC-Hodo;#theta_{xz} of #pi^{+} track for TPC-Hodo",100,-180,180,100,-180,180);
  hist_yzpattern2D2[class_hodo] = new TH2D(Form("hist_yzpattern2D2_%d",class_hodo),"single hit, #theta_{yz} pattern;#theta_{yz} of #pi^{+} track for TPC-Hodo;#theta_{yz} of #pi^{+} track for TPC-Hodo",100,-90,90,100,-90,90);

  int nevent=tree->GetEntries();
  int count_pi_p1,count_pi_p2,count_mu,count_e; //flag for counting particle
  int multiplicity, check_flag, flag_seg_pi_p2, flag_seg_pi_p1;

  int hodo_seg[38]={0};
  int hodo_seg_pi_p1[38]={0};
  int hodo_seg_pi_p2[38]={0};

  double flag_time_pi_p2, flag_time_pi_p1;
  double hodo_hittime_pi_p1[38]={0};
  double hodo_hittime_pi_p2[38]={0};

  double m_pi = 0.139570;
  double p_beam = 2.0;
  double e_beam = sqrt(p_beam*p_beam + m_pi*m_pi);
  TLorentzVector lv_beam(0,0,p_beam,e_beam); //p1
  double theta_xz_pi_p2,theta_xz_pi_p1;
  double theta_yz_pi_p2,theta_yz_pi_p1;
  double E_pi_p2,E_pi_p1;

  int count_m0=0;
  int count_m1=0, count_m1_coin=0, count_m1_pi_p1=0, count_m1_pi_p2=0, count_m1_dummy=0;
  int count_m2=0, count_m2_coin=0, count_m2_multi=0, count_m2_pi_p1=0, count_m2_pi_p2=0, count_m2_dummy=0;
  int dummy=0;

  for(int i=0;i<nevent;i++){
  tree->GetEntry(i);

  count_pi_p1=0;
  count_pi_p2=0;
  count_e=0;
  count_mu=0;
  multiplicity=0;
  check_flag=0;
  flag_seg_pi_p2, flag_seg_pi_p1;
  flag_time_pi_p2=9999.;
  flag_time_pi_p1=9999.;
  theta_xz_pi_p2=9999;
  theta_xz_pi_p1=9999;
  theta_yz_pi_p2=9999;
  theta_yz_pi_p1=9999;

  for(int j=0;j<nEvt;j++){
  if(evttrid[j]==1&&evtpid[j]==211){
  E_pi_p2 = TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]+m_pi*m_pi);
  TLorentzVector lv_pi_p2(evtpx[j],evtpy[j],evtpz[j],E_pi_p2);
  theta_xz_pi_p2 = (evtpx[j]/TMath::Abs(evtpx[j]))*GetAng_xz(lv_beam,lv_pi_p2);
  //theta_yz_pi_p2 = (evtpy[j]/TMath::Abs(evtpy[j]))*GetAng_yz(lv_beam,lv_pi_p2);
  theta_yz_pi_p2 = GetAng_yz(lv_pi_p2);
}
  if(evttrid[j]==2&&evtpid[j]==211){
  E_pi_p1 = TMath::Sqrt(evtpx[j]*evtpx[j]+evtpy[j]*evtpy[j]+evtpz[j]*evtpz[j]+m_pi*m_pi);
  TLorentzVector lv_pi_p1(evtpx[j],evtpy[j],evtpz[j],E_pi_p1);

  theta_xz_pi_p1 =(evtpx[j]/TMath::Abs(evtpx[j]))*GetAng_xz(lv_beam,lv_pi_p1);
  //theta_yz_pi_p1 =(evtpy[j]/TMath::Abs(evtpy[j]))*GetAng_yz(lv_beam,lv_pi_p1);
  theta_yz_pi_p1 =GetAng_yz(lv_pi_p1);
}
}
  for(int j=0;j<nhTof;j++){
  //hist_pi_p2d[class_hodo] -> Fill(tofpid[j]); //PID
  if(window==true&&tofseg[j]==0&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==1&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==2&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==3&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==17&&TMath::Abs(tofposy[j])<window_size) dummy++;
  else if(window==true&&tofseg[j]==18&&TMath::Abs(tofposy[j])<window_size) dummy++;

  else if(toftrid[j]==1&&TMath::Abs(tofpid[j])==211&&tofedep[j]>ecut){ //for selecting events w/ pi
  count_pi_p2++;
  if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
  hodo_seg[32+0] += 1;
  hodo_seg_pi_p2[32+0] += 1;
  hodo_hittime_pi_p2[32+0]=toftime[j];
}
  else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
  hodo_seg[32+1] += 1;
  hodo_seg_pi_p2[32+1] += 1;
  hodo_hittime_pi_p2[32+1]=toftime[j];
}
  else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
  hodo_seg[32+2] += 1;
  hodo_seg_pi_p2[32+2] += 1;
  hodo_hittime_pi_p2[32+2]=toftime[j];
}
  else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
  hodo_seg[32+3] += 1;
  hodo_seg_pi_p2[32+3] += 1;
  hodo_hittime_pi_p2[32+3]=toftime[j];
}
  else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
  hodo_seg[32+4] += 1;
  hodo_seg_pi_p2[32+4] += 1;
  hodo_hittime_pi_p2[32+4]=toftime[j];
}
  else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
  hodo_seg[32+5] += 1;
  hodo_seg_pi_p2[32+5] += 1;
  hodo_hittime_pi_p2[32+5]=toftime[j];
}
  else{
  hodo_seg[tofseg[j]] += 1;
  hodo_seg_pi_p2[tofseg[j]] += 1;
  hodo_hittime_pi_p2[tofseg[j]]=toftime[j];
}
}
  else if(toftrid[j]==2&&tofpid[j]==211&&tofedep[j]>ecut){ //for selecting events w/ p
  count_pi_p1++;
  if(window==true&&tofseg[j]==0&&tofposy[j]<-window_size){
  hodo_seg[32+0] += 1;
  hodo_seg_pi_p1[32+0] += 1;
  hodo_hittime_pi_p1[32+0]=toftime[j];
}
  else if(window==true&&tofseg[j]==1&&tofposy[j]<-window_size){
  hodo_seg[32+1] += 1;
  hodo_seg_pi_p1[32+1] += 1;
  hodo_hittime_pi_p1[32+1]=toftime[j];
}
  else if(window==true&&tofseg[j]==2&&tofposy[j]<-window_size){
  hodo_seg[32+2] += 1;
  hodo_seg_pi_p1[32+2] += 1;
  hodo_hittime_pi_p1[32+2]=toftime[j];
}
  else if(window==true&&tofseg[j]==3&&tofposy[j]<-window_size){
  hodo_seg[32+3] += 1;
  hodo_seg_pi_p1[32+3] += 1;
  hodo_hittime_pi_p1[32+3]=toftime[j];
}
  else if(window==true&&tofseg[j]==17&&tofposy[j]<-window_size){
  hodo_seg[32+4] += 1;
  hodo_seg_pi_p1[32+4] += 1;
  hodo_hittime_pi_p1[32+4]=toftime[j];
}
  else if(window==true&&tofseg[j]==18&&tofposy[j]<-window_size){
  hodo_seg[32+5] += 1;
  hodo_seg_pi_p1[32+5] += 1;
  hodo_hittime_pi_p1[32+5]=toftime[j];
}
  else{
  hodo_seg[tofseg[j]] += 1;
  hodo_seg_pi_p1[tofseg[j]] += 1;
  hodo_hittime_pi_p1[tofseg[j]]=toftime[j];
}
}
  else if(TMath::Abs(tofpid[j])==11&&tofedep[j]>ecut){
  count_e++;
  hodo_seg[tofseg[j]] += 1;
}
  else if(TMath::Abs(tofpid[j])==13&&tofedep[j]>ecut){
  count_mu++;
  hodo_seg[tofseg[j]] += 1;
}
  //std::cout<<"event: "<<i<<std::endl;
}

  if(count_pi_p1>0){
  for(int i=0;i<38;i++){
  if(hodo_seg_pi_p1[i]>0&&flag_time_pi_p1>hodo_hittime_pi_p1[i]){
  flag_time_pi_p1=hodo_hittime_pi_p1[i];
  flag_seg_pi_p1=i;
  //std::cout<<"Proton hodo time: "<<flag_time_pi_p1<<" hodo seg: "<<flag_seg_pi_p1<<std::endl;
}
}
}
  if(count_pi_p2>0){
  for(int i=0;i<38;i++){
  if(hodo_seg_pi_p2[i]>0&&flag_time_pi_p2>hodo_hittime_pi_p2[i]){
  flag_time_pi_p2=hodo_hittime_pi_p2[i];
  flag_seg_pi_p2=i;
  //std::cout<<"Pion hodo time: "<<flag_time_pi_p2<<" hodo seg: "<<flag_seg_pi_p2<<std::endl;
}
}
}
  /*
    if(count_pi_p1>0 && count_pi_p2>0){
    hist_deltatof[class_hodo] -> Fill(flag_time_pi_p1-flag_time_pi_p2);
    if(flag_time_pi_p1 < flag_time_pi_p2) hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi_p1, flag_seg_pi_p2);
    else hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi_p2, flag_seg_pi_p1);
    }
  */
  if(count_pi_p2>0 && count_pi_p1==0){
  hist_xzpattern2D2[class_hodo] -> Fill(theta_xz_pi_p2,theta_xz_pi_p1);
  hist_yzpattern2D2[class_hodo] -> Fill(theta_yz_pi_p2,theta_yz_pi_p1);
}
  if(count_pi_p2>0 && count_pi_p1==0){
  hist_xzpattern2D2[class_hodo] -> Fill(theta_xz_pi_p2,theta_xz_pi_p1);
  hist_yzpattern2D2[class_hodo] -> Fill(theta_yz_pi_p2,theta_yz_pi_p1);

}
  if(count_pi_p2>0 && count_pi_p1>0){
  if(flag_seg_pi_p2<32) flag_seg_pi_p2=(flag_seg_pi_p2+24)%32;
  if(flag_seg_pi_p1<32) flag_seg_pi_p1=(flag_seg_pi_p1+24)%32;

  hist_deltatof[class_hodo] -> Fill(flag_time_pi_p2-flag_time_pi_p1);
  hist_segpattern2D[class_hodo] -> Fill(flag_seg_pi_p2, flag_seg_pi_p1);
  hist_xzpattern2D[class_hodo] -> Fill(theta_xz_pi_p2,theta_xz_pi_p1);
  hist_yzpattern2D[class_hodo] -> Fill(theta_yz_pi_p2,theta_yz_pi_p1);
}
  for(int i=0;i<38;i++){
  if(hodo_seg[i]>0) multiplicity++;
  check_flag += hodo_seg[i];
  hodo_seg[i]=0;
  hodo_seg_pi_p1[i]=0;
  hodo_seg_pi_p2[i]=0;
  hodo_hittime_pi_p1[i]=0;
  hodo_hittime_pi_p2[i]=0;
}

  //counting check
  if(check_flag!=count_pi_p1+count_pi_p2+count_e+count_mu) std::cout<< "fatal error !!! "<<std::endl;

  //M-1
  if(multiplicity==1){
  if(count_pi_p1==1 && count_pi_p2==1 && flag_seg_pi_p2==flag_seg_pi_p1) count_m1_coin++;
  if(count_pi_p1==1 && count_pi_p2==0) count_m1_pi_p1++;
  if(count_pi_p2==1 && count_pi_p1==0) count_m1_pi_p2++;
  if(count_pi_p2==0 && count_pi_p1==0) count_m1_dummy++;
  count_m1++;
}

  //M-2
  else if(multiplicity>1){
  if(count_pi_p1==1 && count_pi_p2==1 && flag_seg_pi_p2==flag_seg_pi_p1) count_m2_coin++;
  if(count_pi_p1>0 && count_pi_p2>0) count_m2_multi++;
  if(count_pi_p1>0 && count_pi_p2==0){
  count_m2_pi_p1++;
}
  if(count_pi_p2>0 && count_pi_p1==0){
  count_m2_pi_p2++;
}
  if(count_pi_p2==0 && count_pi_p1==0) count_m2_dummy++;
  count_m2++;
}

  //M-0
  else count_m0++;
  hist_multi[class_hodo] -> Fill(multiplicity); //Hodoscope multiplicity
}

  //can_pi_p2d[class_hodo] -> cd();
  //hist_pi_p2d[class_hodo]->Draw();
  //hist_pi_p2d[class_hodo]->GetXaxis()->SetTitle("PID(PDG encording) ");
  //hist_pi_p2d[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_multi[class_hodo]->cd();
  hist_multi[class_hodo]->Draw();
  hist_multi[class_hodo]->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi[class_hodo]->GetYaxis()->SetTitle("Counts");

  can_segpattern2D[class_hodo]->cd();
  hist_segpattern2D[class_hodo]->Draw("colz");
  hist_segpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,250);

  can_xzpattern2D[class_hodo]->cd();
  hist_xzpattern2D[class_hodo]->Draw("colz");
  hist_xzpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_yzpattern2D[class_hodo]->cd();
  hist_yzpattern2D[class_hodo]->Draw("colz");
  hist_yzpattern2D[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_xzpattern2D2[class_hodo]->cd();
  hist_xzpattern2D2[class_hodo]->Draw("colz");
  hist_xzpattern2D2[class_hodo]-> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_yzpattern2D2[class_hodo] -> cd();
  hist_yzpattern2D2[class_hodo] -> Draw("colz");
  hist_yzpattern2D2[class_hodo] -> GetZaxis() -> SetRangeUser(0,100);
  gPad -> SetGrid();
  can_deltatof[class_hodo]->cd();
  hist_deltatof[class_hodo]->Draw();

  if(printout==false){
  can_multi[class_hodo] -> Close();
}

  acceptance= (double) count_m2_multi/nevent;

  if(printout){
  std::cout<<""<<std::endl;
  std::cout<<"Event number   : "<<nevent<<std::endl;
  if(nevent!=count_m0+count_m1+count_m2){
  std::cout<<"event number error !!! "<<std::endl;
  std::cout<<"M-0,1,2 number : "<<count_m0+count_m1+count_m2<<std::endl;
}
  //std::cout<<"hit window     : "<<dummy<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"multi-0 number   : "<<(double) count_m0<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"multi-1 number   : "<<count_m1<<std::endl;
  std::cout<<"multi-1 coin     : "<<count_m1_coin<<std::endl;
  std::cout<<"multi-1 one pi hit   : "<<count_m1_pi_p1+count_m1_pi_p2<<std::endl;
  std::cout<<"multi-1 e,mu hit : "<<count_m1_dummy<<std::endl;

  if(count_m1!=count_m1_coin+count_m1_pi_p1+count_m1_pi_p2+count_m1_dummy){
  std::cout<< "M-1 error !!! "<<std::endl;
  std::cout<<"multi-1 check    : "<<count_m1_coin+count_m1_pi_p1+count_m1_pi_p2+count_m1_dummy<<std::endl;
}

  std::cout<<""<<std::endl;
  std::cout<<"multi-2 number   : "<<count_m2<<std::endl;
  std::cout<<"multi-2 coin     : "<<count_m2_coin<<std::endl;
  std::cout<<"multi-2 pi,pi hits : "<<count_m2_multi<<std::endl;
  std::cout<<"multi-2 pi only  : "<<count_m2_pi_p1+count_m2_pi_p2<<std::endl;
  std::cout<<"multi-2 e,mu hits: "<<count_m2_dummy<<std::endl;

  if(count_m2!=count_m2_multi+count_m2_pi_p1+count_m2_pi_p2+count_m2_dummy){
  std::cout<<"M-2 error !!! "<<std::endl;
  std::cout<<"multi-2 check    : "<<count_m2_multi+count_m2_pi_p1+count_m2_pi_p2+count_m2_dummy<<std::endl;
}

  std::cout<<""<<std::endl;
  std::cout<<"acceptance     : "<<(double) acceptance*100<< " %"<<std::endl;
  std::cout<<"multi-0 ratio    : "<<(double) count_m0/nevent*100<< " %"<<std::endl;
  std::cout<<"multi-1 ratio    : "<<(double) count_m1/nevent*100<< " %"<<std::endl;
  std::cout<<"multi-2 ratio    : "<<(double) count_m2/nevent*100<< " %"<<std::endl;

  std::cout<<"m1 coin/m1       : "<<(double) count_m1_coin/count_m1*100<<" %"<<std::endl;
  std::cout<<"m2 coin/m2       : "<<(double) count_m2_coin/count_m2*100<<" %"<<std::endl;

  std::cout<<""<<std::endl;
  std::cout<<"m-0,1,2 ratio    : "<<(double) (count_m0+count_m1+count_m2)/nevent*100<< " %"<<std::endl;
}

  class_hodo++;
}
