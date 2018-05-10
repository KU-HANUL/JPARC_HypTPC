/*
  TofSD.cc from SKS
  2017/10  Yang
*/

#include "TofSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "T01Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"


#include "ConfMan.hh"

const double PositionSeparationThreshold = 2.0*cm;
const double TimeSeparationThreshold     = 5.0*ns;

TofSD::TofSD( G4String name )
  : G4VSensitiveDetector(name)
{
  collectionName.insert( G4String( "TofCollection" ) );
}

TofSD::~TofSD()
{
}

void TofSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  TofCollection = new TofHitsCollection( SensitiveDetectorName,
					 collectionName[0] );
  if( HCID<0 )
    {
      HCID = GetCollectionID(0);
    }
  HCE->AddHitsCollection( HCID, TofCollection );
}

G4bool TofSD::ProcessHits( G4Step *aStep,
			   G4TouchableHistory *ROhist )
{
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.) return false;
  
  G4Track *aTrack = aStep->GetTrack();
  const G4VTouchable *theTouchable = aStep->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume *vol=theTouchable->GetVolume();
  G4String hitName = vol->GetName();
  G4int hitSegment = vol->GetCopyNo();
  
  //G4cout << "[TofSD]TofSD: " << hitName <<  " " << hitSegment << G4endl;

  
  G4int nHits = TofCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime = aTrack->GetGlobalTime();
  G4int trackNo = aTrack->GetTrackID();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4double path = aTrack->GetTrackLength();
  //G4cout<< "pathlength=" << path<<G4endl;
  
  //Decay Particle Tracking
  G4int PartID= aTrack->GetDefinition()->GetPDGEncoding();

  //G4int MID1 = aTrack->GetParentID();
  G4int MID1 = aTrack->GetParentID();
  G4int MID[3];
  G4int MPID[3];

  /*
  if(MID1 != 0 )
    {
      T01Trajectory *parentTrajectory = GetParentTrajectory(trackNo);
      if(parentTrajectory==0)
	{
	  G4cout << "[TofSD] Trajectory trace back failed - no parent found." << G4endl;
	  break;
	}
      else
      {
      MID[0] = MID1;
      MPID[0] = parentTrajectory -> GetPDGEncoding(); 
      }
	}
  */
  for(int i=0 ; i < 3; i++)
    {
      if(MID1 != 0)
	{
	  T01Trajectory* parentTrajectory = GetParentTrajectory(MID1);
	  if(parentTrajectory==0)
	    {
	      G4cout << "[TofSD] Trajectory trace back failed - no parent found." << G4endl;
	      break;
	    }
	  MID[i] = MID1;
	  MPID[i] = parentTrajectory->GetPDGEncoding();
	  MID1 = parentTrajectory->GetParentID();
	  //MID[i+1] = MID1;
	  //std::cout<<"PartID: "<< PartID << " TrackID: " << trackNo << " i: "<< i << " MID: " << MID[i] << " MPID: "<< MPID[i] << std::endl;
	  
	}
      else
	{
	  MID[i] = 0;
	  MPID[i] = 0; 
	}
    }
  
  
  /*
    while(MID!=0)
    {
    
    if(parentTrajectory==0)
    {
    G4cout << "[TofSD] Trajectory trace back failed - no parent found." << G4endl;
    break;
    }
    MID = parentTrajectory->GetParentID();
    //parentTrajectory->ShowTrajectory();
    }
  */


  G4ThreeVector hitposl = theTouchable->GetHistory()->GetTopTransform().TransformPoint( hitpos );
  G4ThreeVector hitmoml = theTouchable->GetHistory()->GetTopTransform().TransformAxis( hitmom );
  //G4double t=hitmoml.z();
  //hitmoml.setZ(-hitmoml.y()); hitmoml.setY(t);

  //G4cout<<"[Tof]PartID: "<< PartID << " MID: "<<MID << " p: "<< hitmom <<G4endl;

  for( G4int i=0; i<nHits; ++i )
    {
      TofHit *aHit = (*TofCollection)[i];
      if( hitSegment==aHit->GetSegmentID() && trackNo == aHit->GetTrackNo())
	{
	  G4double time = aHit->GetTime();
	  if( fabs(hittime-time)<=TimeSeparationThreshold )
	    {
	      aHit->AddEdep( edep );
	      if( hittime<time )
		{
		  aHit->SetTime( hittime );
		  aHit->SetTrackNo( trackNo );
		  aHit->SetPos( hitpos );
		  aHit->SetMom( hitmom );
		  aHit->SetLMom( hitmoml );
		  aHit->SetLocalPos( hitposl.x(), hitposl.z() );
		  aHit->SetPathLength( path );
		  //aHit->SetHitParticleName( PartName );
		}
	      return true;
	    }
	}
    }
  
  TofHit *aHit = new TofHit();
  aHit->SetPass();
  aHit->SetSegmentID( hitSegment );
  aHit->SetEdep( edep );
  aHit->SetTime( hittime );
  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom );
  aHit->SetLMom( hitmoml );
  aHit->SetTrackNo( trackNo );
  aHit->SetLocalPos( hitposl.x(), hitposl.z() );
  aHit->SetPathLength( path );
  aHit->SetHitParticleID( PartID );
  aHit->SetParentID1( MID[0] );
  aHit->SetParentID2( MID[1] );
  aHit->SetParentID3( MID[2] );
  aHit->SetParentPID1( MPID[0] );
  aHit->SetParentPID2( MPID[1] );
  aHit->SetParentPID3( MPID[2] );

  TofCollection->insert( aHit );

#if 0
  G4cout << "[TofSD] " << "Layer=" << hitLayer 
	 << " Seg=" << hitSegment 
	 << " edep=" << edep/keV << "keV"  
	 << " G: " << hitpos << "  L: " << hitposl 
	 <<" P: "<< hitmom << G4endl;
#endif

  return true;
}

void TofSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void TofSD::clear()
{
  G4int nHits = TofCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    {
      delete (*TofCollection)[i];
    }
}

void TofSD::DrawAll() const
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

  if( pVisManager )
    {
      G4int nHits = TofCollection->entries();
      for( G4int i=0; i<nHits; ++i )
	{
	  (*TofCollection)[i]->Draw();
	}
    }
}

void TofSD::PrintAll() const
{
  G4int nHits = TofCollection->entries();
  for( G4int i=0; i<nHits; ++i)
    {
      (*TofCollection)[i]->Print();
    }
}

T01Trajectory* TofSD::GetParentTrajectory(G4int parentID)
{
  G4TrajectoryContainer* container = 
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(container==0) return 0;
  size_t nTraj = container->size();
  for(size_t i=0;i<nTraj;i++)
    {
      T01Trajectory* tr1 = (T01Trajectory*)((*container)[i]);
      if(tr1->GetTrackID()==parentID) return tr1;
    }
  return 0;
}
