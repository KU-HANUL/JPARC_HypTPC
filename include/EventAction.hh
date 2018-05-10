/*
  EventAction.hh

  2017/8  Yang
*/

#ifndef EventAction_h 
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class AnalysisManager;

class EventAction : public G4UserEventAction
{
public:
  EventAction( AnalysisManager *analysisManager=0 );
  ~EventAction();

private:
  EventAction( const EventAction & );
  EventAction & operator = ( const EventAction & );
  
public:
  void BeginOfEventAction( const G4Event *anEvent );
  void EndOfEventAction( const G4Event *anEvent );
  
protected:
  AnalysisManager *anaMan;

};

#endif
