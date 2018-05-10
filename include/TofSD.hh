/*
  TofSD.hh from SKS
  2017/10  Yang
*/

#ifndef TofSD_h
#define TofSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TofHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class T01Trajectory;

class TofSD : public G4VSensitiveDetector
{
public:
  TofSD( G4String name );
  ~TofSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();



private:
  TofHitsCollection *TofCollection;
  T01Trajectory* GetParentTrajectory(G4int parentID);


};

#endif
