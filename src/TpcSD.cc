/*
  TpcSD.cc from SKS
  2017/10  Yang
*/

#include "TpcSD.hh"

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

TpcSD::TpcSD( G4String name )
  : G4VSensitiveDetector(name)
{
  collectionName.insert( G4String( "TpcCollection" ) );
}

TpcSD::~TpcSD()
{
}

void TpcSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  TpcCollection = new TpcHitsCollection( SensitiveDetectorName,
					 collectionName[0] );
  if( HCID<0 )
    {
      HCID = GetCollectionID(0);
    }
  HCE->AddHitsCollection( HCID, TpcCollection );
}

G4bool TpcSD::ProcessHits( G4Step *aStep,
			     G4TouchableHistory *ROhist )
{
  //Select Charged particles
  if(aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0) return false;

  G4double edep = aStep->GetTotalEnergyDeposit();

  G4Track *aTrack = aStep->GetTrack();
  const G4VTouchable *theTouchable = aStep->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume *vol=theTouchable->GetVolume();
  G4String hitName = vol->GetName();
  G4int hitLayer=0;
  G4int hitPad = vol->GetCopyNo();
  
  //G4cout << "TpcSD: " << hitName <<  " " << hitPad << G4endl;
  
  G4int nHits = TpcCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime = aTrack->GetGlobalTime();
  G4int trackNo = aTrack->GetTrackID();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4double path = aTrack->GetTrackLength();
  //G4cout<< "Edep=" << edep<<G4endl;
  
  G4String PartName= aTrack->GetDefinition()->GetParticleName();
  G4int PartID= aTrack->GetDefinition()->GetPDGEncoding();
  G4ThreeVector hitposl = theTouchable->GetHistory()->GetTopTransform().TransformPoint( hitpos );
  G4ThreeVector hitmoml = theTouchable->GetHistory()->GetTopTransform().TransformAxis( hitmom );
  //G4double t=hitmoml.z();
  //hitmoml.setZ(-hitmoml.y()); 
  //hitmoml.setY(t);

  //G4cout<<"particle id: "<< PartName <<" "<< aTrack->GetDefinition()->GetPDGEncoding() << G4endl;
  for( G4int i=0; i<nHits; ++i )
    {
      TpcHit *aHit = (*TpcCollection)[i];
      if( hitPad == aHit->GetPadID() && trackNo == aHit->GetTrackNo())
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

  TpcHit *aHit = new TpcHit();
  aHit->SetPass();
  aHit->SetPadID( hitPad );
  aHit->SetEdep( edep );
  aHit->SetTime( hittime );
  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom );
  aHit->SetLMom( hitmoml );
  aHit->SetTrackNo( trackNo );
  aHit->SetLocalPos( hitposl.x(), hitposl.z() );
  aHit->SetPathLength( path );
  aHit->SetHitParticleID( PartID );
  TpcCollection->insert( aHit );

#if 0
  G4cout << "[TpcSD] " << "Layer=" << hitLayer 
	 << " Pad=" << hitPad 
	 << " edep=" << edep/keV << "keV"  
	 << " G: " << hitpos << "  L: " << hitposl 
	 <<" P: "<< hitmom
         <<" PartID: " << PartID
	 <<" trackNoL: "<< trackNo
	 << G4endl;
#endif

  return true;
}

void TpcSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void TpcSD::clear()
{
  G4int nHits = TpcCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    {
      delete (*TpcCollection)[i];
    }
}

void TpcSD::DrawAll() const
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

  if( pVisManager )
    {
      G4int nHits = TpcCollection->entries();
      for( G4int i=0; i<nHits; ++i )
	{
	  (*TpcCollection)[i]->Draw();
	}
    }
}

void TpcSD::PrintAll() const
{
  G4int nHits = TpcCollection->entries();
  for( G4int i=0; i<nHits; ++i)
    {
      (*TpcCollection)[i]->Print();
    }
}

