/*
  SteppingAction.hh

  2017/9  Yang
*/

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class DetectorConstruction;
class EventAction;

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction( DetectorConstruction*, EventAction* );
  ~SteppingAction();

  void UserSteppingAction(const G4Step*);

private:
  DetectorConstruction* detector;
  EventAction*          eventaction;  

};

#endif
