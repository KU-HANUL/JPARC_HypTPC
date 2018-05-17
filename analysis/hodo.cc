#include <TFile.h>p
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <TPad.h>

using namespace std;

void hodo()
{
  gStyle->SetOptStat(0);
  //data load
  std::cout<<"data chain load"<<std::endl;
  std::cout<<std::endl;

  TFile *file = new TFile("../rootfile/e45/p1035_t1.root","READ");
  TTree *tree = (TTree*)file->Get("tree");
  if(tree == 0) std::cout << "data open error : " << std::endl;

  int event;
  int nEvt;
  int nhTof;
  int tofseg[nhTof];
  double tofedep[nhTof];
  double toftime[nhTof];
  double tofposx[nhTof];
  double tofposy[nhTof];
  double tofposz[nhTof];
  double tofmomx[nhTof];
  double tofmomy[nhTof];
  double tofmomz[nhTof];
  int toftrid[nhTof];
  double tofpath[nhTof];
  int tofpid[nhTof];
  int tofparentpid1[nhTof];
  int tofparentpid2[nhTof];
  int tofparentpid3[nhTof];

  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nEvt",&nEvt);
  tree->SetBranchAddress("nhTof",&nhTof);
  tree->SetBranchAddress("tofseg",tofseg);
  tree->SetBranchAddress("tofedep",tofedep);
  tree->SetBranchAddress("toftime",toftime);
  tree->SetBranchAddress("tofposx",tofposx);
  tree->SetBranchAddress("tofposy",tofposy);
  tree->SetBranchAddress("tofposz",tofposz);
  tree->SetBranchAddress("tofmomx",tofmomx);
  tree->SetBranchAddress("tofmomy",tofmomy);
  tree->SetBranchAddress("tofmomz",tofmomz);
  tree->SetBranchAddress("toftrid",toftrid);
  tree->SetBranchAddress("tofpath",tofpath);
  tree->SetBranchAddress("tofpid",tofpid);
  tree->SetBranchAddress("tofparentpid1",tofparentpid1);
  tree->SetBranchAddress("tofparentpid2",tofparentpid2);
  tree->SetBranchAddress("tofparentpid3",tofparentpid3);

  TCanvas *can_multi = new TCanvas("can_multi","",1200,1200);
  TH1D *hist_multi =new TH1D("hist_multi","Hodoscope Multiplicity",10,0,10);
  TH1D *hist_multi_pi=new TH1D("hist_multi_pi","Hodoscope Multiplicity",10,0,10);

  int dummy;
  int nevent=tree->GetEntries();
  std::cout<<"nevent : "<<std::endl;
  //for(int i=0;i<nevent;i++){
  for(int i=0;i<90000;i++){
    tree -> GetEntry(i);
    dummy=0;
    hist_multi -> Fill(nhTof);
    /*
    for(int j=0;j<nhTof;j++){
      if(tofpid[j]==211||tofpid[j]==-211) dummy++;
    }
    hist_multi_pi -> Fill(dummy);
    */
  }
  can_multi->cd();
  hist_multi->Draw();
  hist_multi_pi->Draw("same");
  hist_multi_pi->SetLineColor(kRed);
  hist_multi->GetXaxis()->SetTitle("Multiplicity of TPC hodo ");
  hist_multi->GetYaxis()->SetTitle("Counts");









  /*

    TCanvas *can_scint = new TCanvas("can_scint","",1200,1200);
  can_scint -> Divide(2,1);

  TH3D *hist_scint_time = new TH3D("hist_scint_time","Scintillator hitpattern",800,-40,40,160,-8,8,1000,-500,500);
  hist_scint_time->GetXaxis()->SetTitle("Hit PositionX(mm) ");
  hist_scint_time->GetYaxis()->SetTitle("Hit PositionY(mm) ");
  hist_scint_time->GetZaxis()->SetTitle("Hit PositionZ(mm) ");

  TH2D *hist_scint_pos = new TH2D("hist_scint_pos","Scintillator hitpattern",800,-40,40,160,-8,8);
  hist_scint_pos->GetXaxis()->SetTitle("Hit PositionX(mm) ");
  hist_scint_pos->GetYaxis()->SetTitle("Hit PositionY(mm) ");

  can_scint->cd(1);
  tree->Draw("scintposz:scintposy:scintposx>>hist_scint_time","scinttime>2");

  can_scint->cd(2);
  tree->Draw("scintposy:scintposx>>hist_scint_pos","scinttime>2","colz");

  TCanvas *can_mppc = new TCanvas("can_mppc","",1200,1200);
  can_mppc -> Divide(2,1);

  TH3D *hist_mppc_time = new TH3D("hist_mppc_time","MPPC hitpattern vs time",800,-40,40,160,-8,8,50,3,8);
  hist_mppc_time->GetXaxis()->SetTitle("Hit PositionX(mm) ");
  hist_mppc_time->GetYaxis()->SetTitle("Hit PositionY(mm) ");
  hist_mppc_time->GetZaxis()->SetTitle("time(ns) ");

  TH2D *hist_mppc_pos = new TH2D("hist_mppc_pos","MPPC hitpattern",700,-35,35,100,-5,5);
  hist_mppc_pos->GetXaxis()->SetTitle("Hit PositionX(mm) ");
  hist_mppc_pos->GetYaxis()->SetTitle("Hit PositionY(mm) ");

  can_mppc->cd(1);
  //tree->Draw("MPPCtime:MPPCposy:MPPCposx>>hist_mppc_time","MPPCtime>2");
  tree->Draw("MPPCtime:MPPCposy:MPPCposx>>hist_mppc_time");

  can_mppc->cd(2);
  //tree->Draw("MPPCposy:MPPCposx>>hist_mppc_pos","MPPCtime>2","colz");
  tree->Draw("MPPCposy:MPPCposx>>hist_mppc_pos","","colz");

  TCanvas *can_slice = new TCanvas("can_slice","",1200,1200);
  can_slice -> Divide(2,2);

  TString cut[5]={"3.0","3.5","4.0","4.5","5.0"};
  TH2D *hist_slice[4];
  for(int i=0;i<4;i++){
    hist_slice[i] = new TH2D(Form("hist_slice%d",i),Form("MPPC hitpattern, %s ns < propagation time < %s ns",cut[i].Data(),cut[i+1].Data()),700,-35,35,100,-5,5);
    can_slice->cd(i+1);
    tree->Draw(Form("MPPCposy:MPPCposx>>hist_slice%d",i),Form("MPPCseg==0 && MPPCtime>%f && MPPCtime<%f",3+0.5*i,3.5+0.5*i),"colz");

    hist_slice[i]->GetXaxis()->SetTitle("Hit PositionX(mm) ");
    hist_slice[i]->GetYaxis()->SetTitle("Hit PositionY(mm) ");
  }

  TCanvas *can_hit = new TCanvas("can_hit","",1200,1200);


  TCanvas *can_time = new TCanvas("can_time","",1200,1200);
  TH2D *hist_time =new TH2D("hist_time","Propagation time:Path length",100,3,5,100,400,600);
  can_time->cd();
  tree->Draw("MPPCpath:MPPCtime>>hist_time","MPPCseg==1","colz");

  hist_time->GetXaxis()->SetTitle("Propagation time (ns) ");
  hist_time->GetYaxis()->SetTitle("Path length (mm) ");

  TCanvas *can_posz = new TCanvas("can_posz","",1200,1200);
  TH1D *hist_posz =new TH1D("hist_posz","Z Position",1000,-500,500);
  can_posz->cd();
  tree->Draw("MPPCposz>>hist_posz","MPPCseg==1");

  hist_posz->GetXaxis()->SetTitle("Z Position (mm) ");
  hist_posz->GetYaxis()->SetTitle("Counts (mm) ");
  */
}//end
