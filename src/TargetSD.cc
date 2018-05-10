/*
  TargetSD.cc
  2017/10  Yang
*/

#include "TargetSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"


#include "ConfMan.hh"

const double PositionSeparationThreshold = 2.0*cm;
const double TimeSeparationThreshold     = 5.0*ns;

TargetSD::TargetSD( G4String name )
  : G4VSensitiveDetector(name)
{
  collectionName.insert( G4String( "TargetCollection" ) );
}

TargetSD::~TargetSD()
{
}

void TargetSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  TargetCollection = new TargetHitsCollection( SensitiveDetectorName,
					       collectionName[0] );
  if( HCID<0 )
    {
      HCID = GetCollectionID(0);
    }
  HCE->AddHitsCollection( HCID, TargetCollection );
}

G4bool TargetSD::ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist )
{

  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;
  
  G4Track *aTrack = aStep->GetTrack();
  const G4VTouchable *theTouchable = aStep->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume *vol=theTouchable->GetVolume();
  G4String hitName = vol->GetName();
  G4int hitLayer=0;
  if( hitName=="TPC_Target" )            hitLayer=0;
  else if( hitName=="TPC_TargetHolder_SideWall1" )    hitLayer=1;
  else if( hitName=="TPC_TargetHolder_SideWall2" )    hitLayer=2;
  else if( hitName=="TPC_TargetHolder_SideWall3" )    hitLayer=3;
  else if( hitName=="TPC_TargetHolder_BottomPlate1" )    hitLayer=11;
  else if( hitName=="TPC_TargetHolder_BottomPlate2" )    hitLayer=12;
  else if( hitName=="TPC_TargetHolder_BottomPlate3" )    hitLayer=13;
  else hitLayer=20;

  G4int hitSegment = vol->GetCopyNo();


  
  G4int nHits = TargetCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime = aTrack->GetGlobalTime();
  G4int trackNo = aTrack->GetTrackID();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4double path = aTrack->GetTrackLength();

  G4String PartName= aTrack->GetDefinition()->GetParticleName();
  G4int PartID= aTrack->GetDefinition()->GetPDGEncoding();
  G4int MID= aTrack->GetParentID();
  G4ThreeVector hitposl = theTouchable->GetHistory()->GetTopTransform().TransformPoint( hitpos );
  G4ThreeVector hitmoml = theTouchable->GetHistory()->GetTopTransform().TransformAxis( hitmom );

  //G4cout<<"[TargetSD]PartID: "<< PartID << " MID: "<<MID<<G4endl;
  for( G4int i=0; i<nHits; ++i )
    {
      TargetHit *aHit = (*TargetCollection)[i];
      if( hitLayer==aHit->GetLayerID() && trackNo == aHit->GetTrackNo())
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

  
  TargetHit *aHit = new TargetHit();
  aHit->SetPass();
  aHit->SetLayerID( hitLayer );
  aHit->SetSegmentID( hitSegment );
  aHit->SetEdep( edep );
  aHit->SetTime( hittime );
  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom );
  aHit->SetLMom( hitmoml );
  aHit->SetTrackNo( trackNo );
  aHit->SetPathLength( path );
  aHit->SetLocalPos( hitposl.x(), hitposl.y() );
  aHit->SetHitParticleID( PartID );
  TargetCollection->insert( aHit );

#if 0
  G4cout << "[TargetSD] " << "Layer=" << hitLayer 
	 << " edep=" << edep/keV << "keV"  
	 << " G: " << hitpos << "  L: " << hitposl 
	 <<" P: "<< hitmom << G4endl;
#endif

  return true;
}

void TargetSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void TargetSD::clear()
{
  G4int nHits = TargetCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    {
      delete (*TargetCollection)[i];
    }
}

void TargetSD::DrawAll()
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

  if( pVisManager )
    {
      G4int nHits = TargetCollection->entries();
      for( G4int i=0; i<nHits; ++i )
	{
	  (*TargetCollection)[i]->Draw();
	}
    }
}

void TargetSD::PrintAll()
{
  G4int nHits = TargetCollection->entries();
  for( G4int i=0; i<nHits; ++i)
    {
      (*TargetCollection)[i]->Print();
    }
}

