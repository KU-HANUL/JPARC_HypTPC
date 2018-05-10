
#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "T01Trajectory.hh"

TrackingAction::TrackingAction()
{
  ;
}

TrackingAction::~TrackingAction()
{
  ;
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  fpTrackingManager->SetStoreTrajectory(true);
  //if(aTrack->GetParentID() == 0)
  //  {
  //    fpTrackingManager->SetStoreTrajectory(true);
  //  }
  //else
  //  {
  //    fpTrackingManager->SetStoreTrajectory(false); }
  //  }
  //Trajectory * trajectory = new Trajectory(aTrack);
  fpTrackingManager->SetTrajectory(new T01Trajectory(aTrack));
  //fpTrackingManager->SetTrajectory( trajectory);
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  ;
}
