/*
  SteppingAction.cc

  2017/8  Yang
*/

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "ConfMan.hh"

SteppingAction::SteppingAction( DetectorConstruction* det,
				EventAction* evt )
  : detector(det), eventaction(evt)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* theStep )
{ 
  ConfMan *confMan = ConfMan::GetConfManager();
  int StepFlag = confMan->StepFlag();

  //if(!StepFlag) return;

  G4Track * theTrack = theStep->GetTrack();

  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(theStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 

  G4ThreeVector pos = theTrack->GetPosition();
  G4StepPoint *preStepPoint = theStep->GetPreStepPoint();
  G4TouchableHandle Touchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector posl = Touchable->GetHistory()->GetTopTransform().TransformPoint( pos );

  G4String particle = theTrack->GetDefinition()->GetParticleName();
  /* 
  G4cout<<"[SteppingAction] particle: "<<particle <<G4endl;
  G4cout<<"[SteppingAction] physVol: "<<physVol->GetName()<<G4endl;
  G4cout<<"[SteppingAction] pos: "<<pos<<G4endl;
  G4cout<<"[SteppingAction] Track ID: "<<theTrack->GetTrackID()<<G4endl;
  */
  // All hits with SC are killed.                                                              
  if( GetFStopSC(StepFlag) && detector->IsVolumeStopper(physVol) )
    {
      theTrack->SetTrackStatus( fStopAndKill );
      //G4cout<<"SC hit killed"<<G4endl;
      return;
    }
  
  // All neutrino are killed.                                                                         
  if( GetFStopNu(StepFlag) &&( particle =="anti_nu_e" || particle =="anti_nu_mu"
			       || particle =="nu_e" || particle =="nu_mu" ) )
    {
      theTrack->SetTrackStatus(fStopAndKill);
      //G4cout<<"neutrino killed"<<G4endl;
      return;
    }
  // Gamma is killed.                                                                                 
  if( GetFStopGam(StepFlag) && particle =="gamma" )
    {
      theTrack->SetTrackStatus(fStopAndKill);
      //G4cout<<"gamma killed"<<G4endl;
      return;
    }
  // e+/e- are killed.                                                                                
  if( GetFStopE(StepFlag) && ( particle =="e-" || particle =="e+" ) )
    {
      theTrack->SetTrackStatus(fStopAndKill);
      //G4cout<<"electron killed"<<G4endl;
      return;
    }


}

