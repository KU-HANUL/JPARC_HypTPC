/*
  Transportation.hh
*/

#ifndef Transportation_h
#define Transportation_h 1

#include "G4Transportation.hh"

class Transportation : public G4Transportation
{
public:
  Transportation() : G4Transportation() {}
  ~Transportation() {}
  
  G4double AlongStepGetPhysicalInteractionLength( const G4Track &track,
                                                  G4double previousStepSize, 
                                                  G4double currentMinimumStep,
                                                  G4double &currerntSafety, 
                                                  G4GPILSelection *selection );
};

#endif
