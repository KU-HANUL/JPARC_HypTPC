//---------------------------------------------------------------
//
// TrackingAction.hh
//
// Description:
//   This class represents actions taken place by the user at each
//   end of stepping. 
//
//---------------------------------------------------------------

#ifndef TrackingAction_h
#define TrackingAction_h 1


#include "G4UserTrackingAction.hh"
//#include "G4TrackingManager.hh"
///////////////////////////

class G4TrackingManager;

class TrackingAction : public G4UserTrackingAction
///////////////////////////
{
  //--------
public:
  //--------
  // Constructor & Destructor
  TrackingAction();
   ~TrackingAction();
  // Member functions
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);
  //----------- 
protected:
  //----------- 
  // Member data
  //G4TrackingManager* fpTrackingManager;
};

#endif
