/*
  RunAction.hh

  2017/8  Yang
*/

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

class G4Run;
class AnalysisManager;

class RunAction : public G4UserRunAction
{
public:
  RunAction( AnalysisManager *analysisManager=0 );
  ~RunAction();

public:
  void BeginOfRunAction( const G4Run *aRun );
  void EndOfRunAction( const G4Run *aRun );

private:
  AnalysisManager *anaMan;
};

#endif
