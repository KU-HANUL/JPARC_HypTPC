//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 18 18:30:16 2018 by ROOT version 6.10/04
// from TTree tree/EvtGen tree
// found on file: p1035_pipip.root
//////////////////////////////////////////////////////////

#ifndef e45_h
#define e45_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class e45 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           nEvt;
   Int_t           evtid[1];   //[nEvt]
   Int_t           evtpid[1];   //[nEvt]
   Int_t           evttrid[1];   //[nEvt]
   Int_t           evtfm[1];   //[nEvt]
   Int_t           evtlm[1];   //[nEvt]
   Int_t           evtfd[1];   //[nEvt]
   Int_t           evtld[1];   //[nEvt]
   Double_t        evttime[1];   //[nEvt]
   Double_t        evtpx[1];   //[nEvt]
   Double_t        evtpy[1];   //[nEvt]
   Double_t        evtpz[1];   //[nEvt]
   Double_t        evtvx[1];   //[nEvt]
   Double_t        evtvy[1];   //[nEvt]
   Double_t        evtvz[1];   //[nEvt]
   Double_t        evtm[1];   //[nEvt]
   Int_t           nBeam;
   Double_t        bpx[1];   //[nBeam]
   Double_t        bpy[1];   //[nBeam]
   Double_t        bpz[1];   //[nBeam]
   Double_t        bvx[1];   //[nBeam]
   Double_t        bvy[1];   //[nBeam]
   Double_t        bvz[1];   //[nBeam]
   Int_t           nhTpcPad;
   Int_t           tpcpad[182];   //[nhTpcPad]
   Double_t        tpcpadedep[182];   //[nhTpcPad]
   Double_t        tpcpadtime[182];   //[nhTpcPad]
   Double_t        tpcpadposx[182];   //[nhTpcPad]
   Double_t        tpcpadposy[182];   //[nhTpcPad]
   Double_t        tpcpadposz[182];   //[nhTpcPad]
   Double_t        tpcpadmomx[182];   //[nhTpcPad]
   Double_t        tpcpadmomy[182];   //[nhTpcPad]
   Double_t        tpcpadmomz[182];   //[nhTpcPad]
   Int_t           tpcpadtrid[182];   //[nhTpcPad]
   Double_t        tpcpadpath[182];   //[nhTpcPad]
   Int_t           tpcpadpid[182];   //[nhTpcPad]
   Int_t           nhTof;
   Int_t           tofseg[10];   //[nhTof]
   Double_t        tofedep[10];   //[nhTof]
   Double_t        toftime[10];   //[nhTof]
   Double_t        tofposx[10];   //[nhTof]
   Double_t        tofposy[10];   //[nhTof]
   Double_t        tofposz[10];   //[nhTof]
   Double_t        tofmomx[10];   //[nhTof]
   Double_t        tofmomy[10];   //[nhTof]
   Double_t        tofmomz[10];   //[nhTof]
   Int_t           toftrid[10];   //[nhTof]
   Double_t        tofpath[10];   //[nhTof]
   Int_t           tofpid[10];   //[nhTof]
   Int_t           tofparentid1[10];   //[nhTof]
   Int_t           tofparentid2[10];   //[nhTof]
   Int_t           tofparentid3[10];   //[nhTof]
   Int_t           tofparentpid1[10];   //[nhTof]
   Int_t           tofparentpid2[10];   //[nhTof]
   Int_t           tofparentpid3[10];   //[nhTof]
   Int_t           nhTarget;
   Int_t           targetseg[10];   //[nhTarget]
   Double_t        targetedep[10];   //[nhTarget]
   Double_t        targettime[10];   //[nhTarget]
   Double_t        targetposx[10];   //[nhTarget]
   Double_t        targetposy[10];   //[nhTarget]
   Double_t        targetposz[10];   //[nhTarget]
   Double_t        targetmomx[10];   //[nhTarget]
   Double_t        targetmomy[10];   //[nhTarget]
   Double_t        targetmomz[10];   //[nhTarget]
   Int_t           targettrid[10];   //[nhTarget]
   Double_t        targetpath[10];   //[nhTarget]
   Int_t           targetpid[10];   //[nhTarget]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_nEvt;   //!
   TBranch        *b_evtid;   //!
   TBranch        *b_evtpid;   //!
   TBranch        *b_evttrid;   //!
   TBranch        *b_evtfm;   //!
   TBranch        *b_evtlm;   //!
   TBranch        *b_evtfd;   //!
   TBranch        *b_evtld;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_evtpx;   //!
   TBranch        *b_evtpy;   //!
   TBranch        *b_evtpz;   //!
   TBranch        *b_evtvx;   //!
   TBranch        *b_evtvy;   //!
   TBranch        *b_evtvz;   //!
   TBranch        *b_evtm;   //!
   TBranch        *b_nBeam;   //!
   TBranch        *b_bpx;   //!
   TBranch        *b_bpy;   //!
   TBranch        *b_bpz;   //!
   TBranch        *b_bvx;   //!
   TBranch        *b_bvy;   //!
   TBranch        *b_bvz;   //!
   TBranch        *b_nhTpcPad;   //!
   TBranch        *b_tpcpad;   //!
   TBranch        *b_tpcpadedep;   //!
   TBranch        *b_tpcpadtime;   //!
   TBranch        *b_tpcpadposx;   //!
   TBranch        *b_tpcpadposy;   //!
   TBranch        *b_tpcpadposz;   //!
   TBranch        *b_tpcpadmomx;   //!
   TBranch        *b_tpcpadmomy;   //!
   TBranch        *b_tpcpadmomz;   //!
   TBranch        *b_tpcpadtrid;   //!
   TBranch        *b_tpcpadpath;   //!
   TBranch        *b_tpcpadpid;   //!
   TBranch        *b_nhTof;   //!
   TBranch        *b_tofseg;   //!
   TBranch        *b_tofedep;   //!
   TBranch        *b_toftime;   //!
   TBranch        *b_tofposx;   //!
   TBranch        *b_tofposy;   //!
   TBranch        *b_tofposz;   //!
   TBranch        *b_tofmomx;   //!
   TBranch        *b_tofmomy;   //!
   TBranch        *b_tofmomz;   //!
   TBranch        *b_toftrid;   //!
   TBranch        *b_tofpath;   //!
   TBranch        *b_tofpid;   //!
   TBranch        *b_tofparentid1;   //!
   TBranch        *b_tofparentid2;   //!
   TBranch        *b_tofparentid3;   //!
   TBranch        *b_tofparentpid1;   //!
   TBranch        *b_tofparentpid2;   //!
   TBranch        *b_tofparentpid3;   //!
   TBranch        *b_nhTarget;   //!
   TBranch        *b_targetseg;   //!
   TBranch        *b_targetedep;   //!
   TBranch        *b_targettime;   //!
   TBranch        *b_targetposx;   //!
   TBranch        *b_targetposy;   //!
   TBranch        *b_targetposz;   //!
   TBranch        *b_targetmomx;   //!
   TBranch        *b_targetmomy;   //!
   TBranch        *b_targetmomz;   //!
   TBranch        *b_targettrid;   //!
   TBranch        *b_targetpath;   //!
   TBranch        *b_targetpid;   //!

   e45(TTree *tree=0);
   virtual ~e45();

   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif

#ifdef e45_cxx
e45::e45(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../rootfiles/e45/p1035_pipip.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("p1035_pipip.root");
      }
      f->GetObject("tree",tree);
   }
   Init(tree);
}

e45::~e45()
{
   if (!fChain) return 0;
   delete fChain->GetCurrentFile();
}

Int_t e45::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void e45::Init(TTree *tree)
{

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nEvt", &nEvt, &b_nEvt);
   fChain->SetBranchAddress("evtid", &evtid, &b_evtid);
   fChain->SetBranchAddress("evtpid", &evtpid, &b_evtpid);
   fChain->SetBranchAddress("evttrid", &evttrid, &b_evttrid);
   fChain->SetBranchAddress("evtfm", &evtfm, &b_evtfm);
   fChain->SetBranchAddress("evtlm", &evtlm, &b_evtlm);
   fChain->SetBranchAddress("evtfd", &evtfd, &b_evtfd);
   fChain->SetBranchAddress("evtld", &evtld, &b_evtld);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("evtpx", &evtpx, &b_evtpx);
   fChain->SetBranchAddress("evtpy", &evtpy, &b_evtpy);
   fChain->SetBranchAddress("evtpz", &evtpz, &b_evtpz);
   fChain->SetBranchAddress("evtvx", &evtvx, &b_evtvx);
   fChain->SetBranchAddress("evtvy", &evtvy, &b_evtvy);
   fChain->SetBranchAddress("evtvz", &evtvz, &b_evtvz);
   fChain->SetBranchAddress("evtm", &evtm, &b_evtm);
   fChain->SetBranchAddress("nBeam", &nBeam, &b_nBeam);
   fChain->SetBranchAddress("bpx", bpx, &b_bpx);
   fChain->SetBranchAddress("bpy", bpy, &b_bpy);
   fChain->SetBranchAddress("bpz", bpz, &b_bpz);
   fChain->SetBranchAddress("bvx", bvx, &b_bvx);
   fChain->SetBranchAddress("bvy", bvy, &b_bvy);
   fChain->SetBranchAddress("bvz", bvz, &b_bvz);
   fChain->SetBranchAddress("nhTpcPad", &nhTpcPad, &b_nhTpcPad);
   fChain->SetBranchAddress("tpcpad", tpcpad, &b_tpcpad);
   fChain->SetBranchAddress("tpcpadedep", tpcpadedep, &b_tpcpadedep);
   fChain->SetBranchAddress("tpcpadtime", tpcpadtime, &b_tpcpadtime);
   fChain->SetBranchAddress("tpcpadposx", tpcpadposx, &b_tpcpadposx);
   fChain->SetBranchAddress("tpcpadposy", tpcpadposy, &b_tpcpadposy);
   fChain->SetBranchAddress("tpcpadposz", tpcpadposz, &b_tpcpadposz);
   fChain->SetBranchAddress("tpcpadmomx", tpcpadmomx, &b_tpcpadmomx);
   fChain->SetBranchAddress("tpcpadmomy", tpcpadmomy, &b_tpcpadmomy);
   fChain->SetBranchAddress("tpcpadmomz", tpcpadmomz, &b_tpcpadmomz);
   fChain->SetBranchAddress("tpcpadtrid", tpcpadtrid, &b_tpcpadtrid);
   fChain->SetBranchAddress("tpcpadpath", tpcpadpath, &b_tpcpadpath);
   fChain->SetBranchAddress("tpcpadpid", tpcpadpid, &b_tpcpadpid);
   fChain->SetBranchAddress("nhTof", &nhTof, &b_nhTof);
   fChain->SetBranchAddress("tofseg", tofseg, &b_tofseg);
   fChain->SetBranchAddress("tofedep", tofedep, &b_tofedep);
   fChain->SetBranchAddress("toftime", toftime, &b_toftime);
   fChain->SetBranchAddress("tofposx", tofposx, &b_tofposx);
   fChain->SetBranchAddress("tofposy", tofposy, &b_tofposy);
   fChain->SetBranchAddress("tofposz", tofposz, &b_tofposz);
   fChain->SetBranchAddress("tofmomx", tofmomx, &b_tofmomx);
   fChain->SetBranchAddress("tofmomy", tofmomy, &b_tofmomy);
   fChain->SetBranchAddress("tofmomz", tofmomz, &b_tofmomz);
   fChain->SetBranchAddress("toftrid", toftrid, &b_toftrid);
   fChain->SetBranchAddress("tofpath", tofpath, &b_tofpath);
   fChain->SetBranchAddress("tofpid", tofpid, &b_tofpid);
   fChain->SetBranchAddress("tofparentid1", tofparentid1, &b_tofparentid1);
   fChain->SetBranchAddress("tofparentid2", tofparentid2, &b_tofparentid2);
   fChain->SetBranchAddress("tofparentid3", tofparentid3, &b_tofparentid3);
   fChain->SetBranchAddress("tofparentpid1", tofparentpid1, &b_tofparentpid1);
   fChain->SetBranchAddress("tofparentpid2", tofparentpid2, &b_tofparentpid2);
   fChain->SetBranchAddress("tofparentpid3", tofparentpid3, &b_tofparentpid3);
   fChain->SetBranchAddress("nhTarget", &nhTarget, &b_nhTarget);
   fChain->SetBranchAddress("targetseg", targetseg, &b_targetseg);
   fChain->SetBranchAddress("targetedep", targetedep, &b_targetedep);
   fChain->SetBranchAddress("targettime", targettime, &b_targettime);
   fChain->SetBranchAddress("targetposx", targetposx, &b_targetposx);
   fChain->SetBranchAddress("targetposy", targetposy, &b_targetposy);
   fChain->SetBranchAddress("targetposz", targetposz, &b_targetposz);
   fChain->SetBranchAddress("targetmomx", targetmomx, &b_targetmomx);
   fChain->SetBranchAddress("targetmomy", targetmomy, &b_targetmomy);
   fChain->SetBranchAddress("targetmomz", targetmomz, &b_targetmomz);
   fChain->SetBranchAddress("targettrid", targettrid, &b_targettrid);
   fChain->SetBranchAddress("targetpath", targetpath, &b_targetpath);
   fChain->SetBranchAddress("targetpid", targetpid, &b_targetpid);
}


#endif // #ifdef e45_cxx
