/*
  RunAction.cc

  2017/8  Yang
*/

#include "RunAction.hh"
#include "AnalysisManager.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

RunAction::RunAction( AnalysisManager *analysisManager )
  : G4UserRunAction(), anaMan(analysisManager)
{}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction( const G4Run *aRun )
{
  if( G4VVisManager::GetConcreteInstance() ){
    G4UImanager *UI = G4UImanager::GetUIpointer();
    //    UI->ApplyCommand( "/vis/scene/notifyHandlers" );
  }
  if( anaMan ) anaMan->BeginOfRun( aRun );
}

void RunAction::EndOfRunAction( const G4Run *aRun )
{
  if( G4VVisManager::GetConcreteInstance() ){
    G4UImanager *UI = G4UImanager::GetUIpointer();
    //    UI->ApplyCommand( "/vis/viewer/update" );
  }
  
  if( anaMan ) anaMan->EndOfRun( aRun );
} 

