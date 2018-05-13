###########

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "ConfMan.hh"
#include "TrackingAction.hh"


#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtStdlibRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenExternal/EvtExternalGenFactory.hh"


//////////////////////////
//<mac> <conf> <outfile>//
//////////////////////////

int main(int argc, char** argv)
{
  CLHEP::HepRandom::setTheSeed((unsigned)time(NULL));

  //Conf file//
  G4String confFile;
  if( argc>=3 )
    {
      confFile = argv[2];
    }
  else
    {
      //confFile = "conf/accept/p760_t0.6_electron.conf";
      //confFile = "conf/trigger/p730_zm650_f1.0.conf";
      //confFile = "conf/trigger/p730_zm1000_f0.6.conf";
      confFile = "conf/p730_t1.conf";
      //confFile = "conf/beam.conf";
    }
  ConfMan *confManager = new ConfMan( confFile );
  confManager->Initialize();

  //Root file//
  G4String histname;
  if( argc >= 4 )
    {
      histname = argv[3];
    }
  else
    {
      histname = "geant4.root";
    }

  //EvtGen//
  G4cout<< "************" <<G4endl;
  EvtStdlibRandomEngine eng;
  EvtRandom::setRandomEngine((EvtRandomEngine*)&eng);

  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;

  //#ifdef EVTGEN_EXTERNAL
  EvtExternalGenList genList;
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
  //#endif

  EvtExternalGenFactory* externalGenerators = EvtExternalGenFactory::getInstance();
  std::string xmlDir("./param/xmldoc");

  bool convertPhysCode(true);
  externalGenerators->definePythiaGenerator(xmlDir, convertPhysCode);

  EvtGen *evtgen(0);
  ConfMan *confMan = ConfMan::GetConfManager();
  if( confMan->ExistEvtData1() && confMan->ExistEvtData1()){
    //Initialize the generator - read in the decay table and particle properties
    evtgen = new EvtGen( confMan->EvtGenDecayName().c_str(),
                         confMan->EvtGenPDLName().c_str(),
                         (EvtRandomEngine*)&eng,
                         radCorrEngine, &extraModels);
  }
  G4cout<< "************" <<G4endl;


  G4RunManager * runManager = new G4RunManager;

  DetectorConstruction *detector = new DetectorConstruction();
  runManager->SetUserInitialization( detector );
  //runManager->SetUserInitialization(new DetectorConstruction);

  G4cout<<"Test1"<<G4endl;
  PhysicsList *physList = new PhysicsList();
  runManager->SetUserInitialization( physList );
  //runManager->SetUserInitialization(new PhysicsList);
  runManager->Initialize();


  AnalysisManager *anaMan = new AnalysisManager( histname );

  PrimaryGeneratorAction *priGen = new PrimaryGeneratorAction( detector, anaMan, evtgen );
  //PrimaryGeneratorAction *priGen = new PrimaryGeneratorAction( detector, anaMan);
  RunAction *runAction = new RunAction( anaMan );
  EventAction *eventAction = new EventAction( anaMan );
  SteppingAction* stepAction = new SteppingAction( detector, eventAction );


  runManager->SetUserAction(priGen);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(stepAction);

  TrackingAction* tracking_action = new TrackingAction;
  runManager->SetUserAction(tracking_action);


  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive *ui = 0;
  if (argc == 1)
    {
      ui = new G4UIExecutive(argc, argv);
    }

  if (!ui)   // batch mode
    {
      //G4String fileName1 = argv[3];
      //AnalysisManager::GetInstance()->SetRootFileName(fileName1);

      G4String command = "/control/execute ";
      G4String fileName2 = argv[1];
      UImanager->ApplyCommand(command+fileName2);
    }
  else  // interactive mode : define UI session
    {
      //AnalysisManager::GetInstance()->SetRootFileName("default.root");
      //G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      G4cout<<"Test2"<<G4endl;
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    }

  delete visManager;
  delete runManager;

  return 0;
}
