/*
  EventAction.cc

  2017/8  Yang
*/

#include "EventAction.hh"
#include "AnalysisManager.hh"

EventAction::EventAction( AnalysisManager *analysisManager )
  : G4UserEventAction(), anaMan(analysisManager)
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction( const G4Event *anEvent )
{
  //G4int eventID = anEvent->GetEventID();
  if (anaMan) anaMan->BeginOfEvent( anEvent );
}

void EventAction::EndOfEventAction( const G4Event *anEvent )
{
  if (anaMan) anaMan->EndOfEvent(anEvent);
}
