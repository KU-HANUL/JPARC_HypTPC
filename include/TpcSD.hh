/*
  TpcSD.hh from SKS
  2017/10  Yang
*/

#ifndef TpcSD_h
#define TpcSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TpcHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class TpcSD : public G4VSensitiveDetector
{
public:
  TpcSD( G4String name );
  ~TpcSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  TpcHitsCollection *TpcCollection;
};

#endif
