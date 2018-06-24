#include "AnalysisManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "TpcHit.hh"
#include "TofHit.hh"
#include "TargetHit.hh"
//#include "TrackerHit.hh"
#include "Randomize.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TString.h"

#include "EvtGenBase/EvtVector4R.hh"

//#include "RootHelper.hh"
#include <string>
#include <sstream>

#include "TF1.h"
#include "TMath.h"

AnalysisManager::AnalysisManager( const G4String & histname )
  :outfile(histname), fActive_(true)
{
}

AnalysisManager::~AnalysisManager()
{
  SaveFile();
}

void AnalysisManager::SaveFile( void ) const
{
  if( fActive_ )
    hfile->Write();
}

void AnalysisManager::Terminate ( void ) const
{
  G4cout << "[Analysis] Terminate()"<< G4endl;
  if ( fActive_ )
    {
      hfile->Write();
      hfile->Close();
    }
}

void AnalysisManager::BeginOfRun(const G4Run*)
{
  G4cout<<"[AnalysisManager] Begin of Run: "<<G4endl;
  G4SDManager* SDManager = G4SDManager::GetSDMpointer();

  hfile = new TFile(outfile, "RECREATE");
  tree = new TTree("tree","EvtGen tree");

  G4cout<<"Output file made: "<< outfile << G4endl;

  tree->Branch("event", &event, "event/I");
  tree->Branch("nEvt", &nEvt, "nEvt/I");
  tree->Branch("evtid", evtid, "evtid[nEvt]/I");
  tree->Branch("evtpid", evtpid, "evtpid[nEvt]/I");
  tree->Branch("evttrid", evttrid, "evttrid[nEvt]/I");
  tree->Branch("evtfm", evtfm, "evtfm[nEvt]/I");
  tree->Branch("evtlm", evtlm, "evtlm[nEvt]/I");
  tree->Branch("evtfd", evtfd, "evtfd[nEvt]/I");
  tree->Branch("evtld", evtld, "evtld[nEvt]/I");
  tree->Branch("evttime", evttime, "evttime[nEvt]/D");
  tree->Branch("evtpx", evtpx, "evtpx[nEvt]/D");
  tree->Branch("evtpy", evtpy, "evtpy[nEvt]/D");
  tree->Branch("evtpz", evtpz, "evtpz[nEvt]/D");
  tree->Branch("evtvx", evtvx, "evtvx[nEvt]/D");
  tree->Branch("evtvy", evtvy, "evtvy[nEvt]/D");
  tree->Branch("evtvz", evtvz, "evtvz[nEvt]/D");
  tree->Branch("evtm", evtm, "evtm[nEvt]/D");

  tree->Branch("nBeam", &nBeam, "nBeam/I");
  tree->Branch("bpx", bpx, "bpx[nBeam]/D");
  tree->Branch("bpy", bpy, "bpy[nBeam]/D");
  tree->Branch("bpz", bpz, "bpz[nBeam]/D");
  tree->Branch("bvx", bvx, "bvx[nBeam]/D");
  tree->Branch("bvy", bvy, "bvy[nBeam]/D");
  tree->Branch("bvz", bvz, "bvz[nBeam]/D");

  tree->Branch("nhTpcPad", &nhTpcPad, "nhTpcPad/I");
  tree->Branch("tpcpad", tpcpad, "tpcpad[nhTpcPad]/I");
  tree->Branch("tpcpadedep", tpcpadedep, "tpcpadedep[nhTpcPad]/D");
  tree->Branch("tpcpadtime", tpcpadtime, "tpcpadtime[nhTpcPad]/D");
  tree->Branch("tpcpadposx", tpcpadposx, "tpcpadposx[nhTpcPad]/D");
  tree->Branch("tpcpadposy", tpcpadposy, "tpcpadposy[nhTpcPad]/D");
  tree->Branch("tpcpadposz", tpcpadposz, "tpcpadposz[nhTpcPad]/D");
  tree->Branch("tpcpadmomx", tpcpadmomx, "tpcpadmomx[nhTpcPad]/D");
  tree->Branch("tpcpadmomy", tpcpadmomy, "tpcpadmomy[nhTpcPad]/D");
  tree->Branch("tpcpadmomz", tpcpadmomz, "tpcpadmomz[nhTpcPad]/D");
  tree->Branch("tpcpadtrid", tpcpadtrid, "tpcpadtrid[nhTpcPad]/I");
  tree->Branch("tpcpadpath", tpcpadpath, "tpcpadpath[nhTpcPad]/D");
  tree->Branch("tpcpadpid", tpcpadpid, "tpcpadpid[nhTpcPad]/I");

  tree->Branch("nhTof", &nhTof, "nhTof/I");
  tree->Branch("tofseg", tofseg, "tofseg[nhTof]/I");
  tree->Branch("tofedep", tofedep, "tofedep[nhTof]/D");
  tree->Branch("toftime", toftime, "toftime[nhTof]/D");
  tree->Branch("tofposx", tofposx, "tofposx[nhTof]/D");
  tree->Branch("tofposy", tofposy, "tofposy[nhTof]/D");
  tree->Branch("tofposz", tofposz, "tofposz[nhTof]/D");
  tree->Branch("tofmomx", tofmomx, "tofmomx[nhTof]/D");
  tree->Branch("tofmomy", tofmomy, "tofmomy[nhTof]/D");
  tree->Branch("tofmomz", tofmomz, "tofmomz[nhTof]/D");
  tree->Branch("toftrid", toftrid, "toftrid[nhTof]/I");
  tree->Branch("tofpath", tofpath, "tofpath[nhTof]/D");
  tree->Branch("tofvtxx", tofvtxx, "tofvtxx[nhTof]/D");
  tree->Branch("tofvtxy", tofvtxy, "tofvtxy[nhTof]/D");
  tree->Branch("tofvtxz", tofvtxz, "tofvtxz[nhTof]/D");
  tree->Branch("tofpid", tofpid, "tofpid[nhTof]/I");
  tree->Branch("tofparentid1", tofparentid1, "tofparentid1[nhTof]/I");
  tree->Branch("tofparentid2", tofparentid2, "tofparentid2[nhTof]/I");
  tree->Branch("tofparentid3", tofparentid3, "tofparentid3[nhTof]/I");
  tree->Branch("tofparentpid1", tofparentpid1, "tofparentpid1[nhTof]/I");
  tree->Branch("tofparentpid2", tofparentpid2, "tofparentpid2[nhTof]/I");
  tree->Branch("tofparentpid3", tofparentpid3, "tofparentpid3[nhTof]/I");

  tree->Branch("nhTarget", &nhTarget, "nhTarget/I");
  tree->Branch("targetseg", targetseg, "targetseg[nhTarget]/I");
  tree->Branch("targetedep", targetedep, "targetedep[nhTarget]/D");
  tree->Branch("targettime", targettime, "targettime[nhTarget]/D");
  tree->Branch("targetposx", targetposx, "targetposx[nhTarget]/D");
  tree->Branch("targetposy", targetposy, "targetposy[nhTarget]/D");
  tree->Branch("targetposz", targetposz, "targetposz[nhTarget]/D");
  tree->Branch("targetmomx", targetmomx, "targetmomx[nhTarget]/D");
  tree->Branch("targetmomy", targetmomy, "targetmomy[nhTarget]/D");
  tree->Branch("targetmomz", targetmomz, "targetmomz[nhTarget]/D");
  tree->Branch("targettrid", targettrid, "targettrid[nhTarget]/I");
  tree->Branch("targetpath", targetpath, "targetpath[nhTarget]/D");
  tree->Branch("targetpid", targetpid, "targetpid[nhTarget]/I");

  event = 0;


}


void AnalysisManager::EndOfRun(const G4Run*)
{
  G4cout<<"[AnalysisManager] End of Event: "<<G4endl;
  tree->Write();
  hfile->Write();
  hfile->Close();
}

void AnalysisManager::BeginOfEvent(const G4Event*)
{
  //nEvt = 0;
  G4cout<<"[AnalysisManager] Begin of Event: "<<nEvt<<G4endl;
}

void AnalysisManager::EndOfEvent(const G4Event* anEvent)
{

  G4HCofThisEvent* HCTE = anEvent-> GetHCofThisEvent();
  if(!HCTE) return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhtpcpad=0, nhtof=0, nhtarget=0;
  TpcHitsCollection *TpcHC=0;
  TofHitsCollection *TofHC=0;
  TargetHitsCollection *TargetHC=0;

  G4int colIdTpc = SDMan->GetCollectionID( "TpcCollection" );
  if(colIdTpc>=0)
    {
      TpcHC=dynamic_cast<TpcHitsCollection *>( HCTE->GetHC( colIdTpc ) );
      if(TpcHC)
	{
	  nhtpcpad=TpcHC->entries();
	}
    }

  G4int colIdTof = SDMan->GetCollectionID( "TofCollection" );
  if(colIdTof>=0)
    {
      TofHC=dynamic_cast<TofHitsCollection *>( HCTE->GetHC( colIdTof ) );
      if(TofHC)
	{
	  nhtof=TofHC->entries();
	}
    }

  G4int colIdTarget = SDMan->GetCollectionID( "TargetCollection" );
  if(colIdTarget>=0)
    {
      TargetHC=dynamic_cast<TargetHitsCollection *>( HCTE->GetHC( colIdTarget ) );
      if(TargetHC)
	{
	  nhtarget=TargetHC->entries();
	}
    }



  if(nhtpcpad > 200 -1)
    {
      G4cout<<"[AnalysisManager] Number of TPC Hit > 200"<<G4endl;
      return;
    }
  if(nhtof > 100 -1)
    {
      G4cout<<"[AnalysisManager] Number of Tof Hit > 100"<<G4endl;
      return;
    }
  if(nhtarget > 100 -1)
    {
      G4cout<<"[AnalysisManager] Number of Target Hit > 100"<<G4endl;
      return;
    }

  for( int i=0; i<nhtpcpad; ++i )
    {
      TpcHit *aHit=(*TpcHC)[i];
      tpcpad[i] = aHit->GetPadID();
      tpcpadedep[i] = aHit->GetEdep();
      tpcpadtime[i] = aHit->GetTime();
      tpcpadposx[i] = aHit->GetPos().x();
      tpcpadposy[i] = aHit->GetPos().y();
      tpcpadposz[i] = aHit->GetPos().z();
      tpcpadmomx[i] = aHit->GetMom().x()*0.001;
      tpcpadmomy[i] = aHit->GetMom().y()*0.001;
      tpcpadmomz[i] = aHit->GetMom().z()*0.001;
      tpcpadtrid[i] = aHit->GetTrackNo();
      tpcpadpath[i] = aHit->GetPathLength();
      tpcpadpid[i] = aHit->GetHitParticleID();
    }
  nhTpcPad = nhtpcpad;

  for( int i=0; i<nhtof; ++i )
    {
      TofHit *aHit=(*TofHC)[i];
      tofseg[i] = aHit->GetSegmentID();
      tofedep[i] = aHit->GetEdep();
      toftime[i] = aHit->GetTime();
      tofposx[i] = aHit->GetPos().x();
      tofposy[i] = aHit->GetPos().y();
      tofposz[i] = aHit->GetPos().z();
      tofmomx[i] = aHit->GetMom().x()*0.001;
      tofmomy[i] = aHit->GetMom().y()*0.001;
      tofmomz[i] = aHit->GetMom().z()*0.001;
      toftrid[i] = aHit->GetTrackNo();
      tofpath[i] = aHit->GetPathLength();
      tofvtxx[i] = aHit->GetVtxPos().x();
      tofvtxy[i] = aHit->GetVtxPos().y();
      tofvtxz[i] = aHit->GetVtxPos().z();
      tofpid[i] = aHit->GetHitParticleID();
      tofparentid1[i] = aHit->GetParentID1();
      tofparentid2[i] = aHit->GetParentID2();
      tofparentid3[i] = aHit->GetParentID3();
      tofparentpid1[i] = aHit->GetParentPID1();
      tofparentpid2[i] = aHit->GetParentPID2();
      tofparentpid3[i] = aHit->GetParentPID3();

    }
  nhTof = nhtof;

  for( int i=0; i<nhtarget; ++i )
    {
      TargetHit *aHit=(*TargetHC)[i];
      targetseg[i] = aHit->GetLayerID();
      targetedep[i] = aHit->GetEdep();
      targettime[i] = aHit->GetTime();
      targetposx[i] = aHit->GetPos().x();
      targetposy[i] = aHit->GetPos().y();
      targetposz[i] = aHit->GetPos().z();
      targetmomx[i] = aHit->GetMom().x()*0.001;
      targetmomy[i] = aHit->GetMom().y()*0.001;
      targetmomz[i] = aHit->GetMom().z()*0.001;
      targettrid[i] = aHit->GetTrackNo();
      targetpath[i] = aHit->GetPathLength();
      targetpid[i] = aHit->GetHitParticleID();
    }
  nhTarget = nhtarget;

  tree->Fill();
  event++;

  G4cout<<"[AnalysisManager]Event: "<<event<<" nEvt: "<<nEvt << " nhTarget: " << nhTarget<<G4endl;
  nEvt = 0;
  nhTpcPad = 0;
  nhTof = 0;
  nhTarget = 0;
  nBeam = 0;
}

void AnalysisManager::SetEvtGen(int j, int partnum,
				int jmotherfirst,
				int jmotherlast,
				int jdaugfirst,
				int jdauglast,
				int tr_id,
				EvtVector4R x,
				EvtVector4R p
				)
{
  G4cout<<"[AnalysisManager] setEvtGen : "<<nEvt<<G4endl;
  if(nEvt < j) nEvt =j;
  evtid[j-1]=j;
  evtpid[j-1] = partnum;
  evttrid[j-1] = tr_id;
  evtfm[j-1] = jmotherfirst;
  evtlm[j-1] = jmotherlast;
  evtfd[j-1] = jdaugfirst;
  evtld[j-1] = jdauglast;
  evttime[j-1] = x.get(0);
  evtpx[j-1] = p.get(1);
  evtpy[j-1] = p.get(2);
  evtpz[j-1] = p.get(3);
  evtvx[j-1] = x.get(1);
  evtvy[j-1] = x.get(2);
  evtvz[j-1] = x.get(3);
  evtm[j-1] = p.mass();
  //G4cout<<"[Analysis Manager] j: "<<j << " ParticleNum: " << partnum << " nEvt: "<<nEvt<<G4endl;
}

void AnalysisManager::SetBeam(int j,
			      G4ThreeVector D,
			      G4ThreeVector P
			      )
{
  //xG4cout<<"[AnalysisManager] SetBeam : "<<nEvt<<G4endl;
  nBeam = j;
  bpx[j-1] = P.x();
  bpy[j-1] = P.y();
  bpz[j-1] = P.z();
  bvx[j-1] = D.x();
  bvy[j-1] = D.y();
  bvz[j-1] = D.z();

}
