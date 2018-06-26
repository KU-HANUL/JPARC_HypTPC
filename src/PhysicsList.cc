#include "PhysicsList.hh"

#include "G4VModularPhysicsList.hh"
//#include "G4Scintillation.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

#include "G4SystemOfUnits.hh"

#include "globals.hh"
#include "PhysicsList.hh"
#include "Transportation.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4Ions.hh"
#include "G4ios.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"
#include "G4Threading.hh"

#include <iomanip>

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4FastSimulationManagerProcess.hh"

//#include "ConfMan.hh"


PhysicsList::PhysicsList() : G4VModularPhysicsList(){
  defaultCutValue = 0.1*cm;
}
////////////////////////////////////////////////////////

PhysicsList::~PhysicsList(){}

void PhysicsList::ConstructParticle(){
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  //ConstructHeavyIon();
}

void PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  // gamma
  G4Gamma::GammaDefinition();
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void PhysicsList::ConstructLeptons()
{
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void PhysicsList::ConstructMesons()
{
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();

  //Particles
  G4DecayTable *decayTable;
  G4VDecayChannel *mode;
  G4ParticleDefinition *particle;

  //test
  particle
    = new G4ParticleDefinition("kaonT1",    0.493677*GeV,  5.352e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,        -321,
			       false,         1.238e-8*s,          NULL,
			       false,       "kaon");

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaonT1", 1.0, 2, "mu-", "anti_nu_mu" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );



  //K- -> pi- pi0
  #if 0
  particle
    = new G4ParticleDefinition("kaonM1",    0.493677*GeV,  5.352e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,        -321,
			       true,         -1.0,          NULL,
			       false,       "kaon");

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaonM1", 1.0, 2, "pi-", "pi0" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //K- -> mu- an-nu-mu
  particle
    = new G4ParticleDefinition("kaonM2",    0.493677*GeV,  5.352e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,         -321,
			       true,         -1.0,          NULL,
			       false,       "kaon");

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaonM2", 1.0, 2, "mu-", "anti_nu_mu" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //pi- no decay
  particle
    = new G4ParticleDefinition("pionM1",    0.1395700*GeV,  2.5452e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       2,              -2,             -1,
			       "meson",         0,             0,         -211,
			       true,         -1.0,          NULL,
			       false,       "pi");

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "pionM1", 1.0, 2, "mu-", "anti_nu_mu" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //pi+ no decay
  particle
    = new G4ParticleDefinition("pionP1",    0.1395700*GeV,  2.5452e-14*MeV,    +1.*eplus,
			       0,              -1,             0,
			       2,              +2,             -1,
			       "meson",         0,             0,          211,
			       true,         -1.0,          NULL,
			       false,       "pi");

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "pionP1", 1.0, 2, "mu+", "nu_mu" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //K->Pi+ Pi- Pi-
  particle
    = new G4ParticleDefinition(
			       "kaon1-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,            -321,
			       false,          0.0,           NULL,
			       false,         "kaon" );

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaon1-", 1.0, 3, "pi+", "pi-", "pi-" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //K->e- Pi0 nue
  particle
    = new G4ParticleDefinition (
				"kaon2-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus,
				0,              -1,             0,
				1,              -1,             0,
				"meson",         0,             0,           -321,
				false,         0.0,            NULL,
				false,        "kaon" );

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("kaon2-",1.0, 3,
				       "e-","pi0", "anti_nu_e");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //K->Mu- Pi0 nuMu
  particle
    = new G4ParticleDefinition (
				"kaon3-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus,
				0,              -1,             0,
				1,              -1,             0,
				"meson",         0,             0,            -321,
				false,         0.0,            NULL,
				false,        "kaon" );

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("kaon3-",1.0, 3,
				       "mu-","pi0", "anti_nu_mu");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //K->Pi- Pi0 Pi0
  particle
    = new G4ParticleDefinition (
				"kaon4-",    0.493677*GeV, 5.315e-14*MeV,    -1.*eplus,
				0,              -1,             0,
				1,              -1,             0,
				"meson",         0,             0,            -321,
				false,         0.0,            NULL,
				false,        "kaon" );

  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel("kaon4-",1.0, 3,
				       "pi-","pi0", "pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  #endif


}

void PhysicsList::ConstructBaryons(){
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

void PhysicsList::ConstructHeavyIon(){
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


void PhysicsList::ConstructProcess(){
  //ConfMan *confMan = ConfMan::GetConfManager();
  //int flag = confMan->PhysFlag();

  AddTransportationSks();
  ConstructEM();
  //ConstructOp();
  ConstructHadronic();
  ConstructDecay();


  G4StepLimiter *stepLimiter = new G4StepLimiter();
  G4UserSpecialCuts *userCuts = new G4UserSpecialCuts();
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager ->AddDiscreteProcess(stepLimiter);
    pmanager ->AddDiscreteProcess(userCuts);
  }
}

///////Transportation
void PhysicsList::AddTransportationSks(){
  Transportation* theTransportationProcess= new Transportation();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if(!particle->IsShortLived()){
      pmanager->AddProcess(theTransportationProcess);
      pmanager->SetProcessOrderingToFirst(theTransportationProcess,idxAlongStep);
      pmanager ->SetProcessOrderingToFirst(theTransportationProcess,idxPostStep);
    }
  }
}

///////Cut
void PhysicsList::SetCuts()
{
  // Suppress error message int case e/gamma/proton do not exist
  G4int temp = GetVerboseLevel();
  // Retrive verbose level
  SetVerboseLevel(temp);

  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  SetCutsWithDefault();
  //  theParticleIterator->reset();
  //  while( (*theParticleIterator)() ){
  //    G4ParticleDefinition *particle=theParticleIterator->value();
  //    particle->SetApplyCutsFlag( true );
    //////////////////////////////////////////////////////////////////////
    //    G4cout << particle->GetParticleName() << " ==> ApplyCutFlag = "
    //     << particle->GetApplyCutsFlag() << G4endl;
    //////////////////////////////////////////////////////////////////////
  //  }
}

///////EM
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"

//#include "G4VEnergyLossProcess.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4UserSpecialCuts.hh"

void PhysicsList::ConstructEM(){
  /*
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  RegisterPhysics( opticalPhysics );
  opticalPhysics->SetScintillationYieldFactor(1.0);
  opticalPhysics->SetScintillationExcitationRatio(0.0);
  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  */
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if( particleName == "gamma" ){
      pManager->AddProcess( new G4UserSpecialCuts(), -1, -1, 1 );
      pManager->AddDiscreteProcess( new G4PhotoElectricEffect() );
      pManager->AddDiscreteProcess( new G4ComptonScattering() );
      pManager->AddDiscreteProcess( new G4GammaConversion() );
    }
    else if( particleName == "e-" ){
      pManager->AddProcess( new G4UserSpecialCuts(),     -1, -1, 1 );
      pManager->AddProcess( new G4eMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4eIonisation(),         -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung(),     -1, -1, 4 );
    }
    else if( particleName == "e+" ){
      pManager->AddProcess( new G4UserSpecialCuts(),     -1, -1, 1 );
      pManager->AddProcess( new G4eMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4eIonisation(),         -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung(),     -1, -1, 4 );
      pManager->AddProcess( new G4eplusAnnihilation(),    1, -1, 5 );
    }
    else if( particleName == "mu+" || particleName == "mu-" ){
      pManager->AddProcess( new G4MuMultipleScattering, -1,  1, 1 );
      pManager->AddProcess( new G4MuIonisation,       -1,  2, 2 );
      pManager->AddProcess( new G4MuBremsstrahlung,   -1,  3, 3 );
      pManager->AddProcess( new G4MuPairProduction,   -1,  4, 4 );
      /*
      pManager->AddProcess( new G4UserSpecialCuts(),      -1, -1, 1 );
      pManager->AddProcess( new G4MuMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4MuIonisation(),         -1,  2, 3 );
      pManager->AddProcess( new G4MuBremsstrahlung(),     -1,  3, 4 );
      pManager->AddProcess( new G4MuPairProduction(),     -1,  4, 5 );
      */
    }
    else if( !(particle->IsShortLived()) && particle->GetPDGCharge()!=0 &&
	     !( particleName=="chargedgeantino"
		|| particleName=="antichargedgeantino") ){
      pManager->AddProcess(new G4UserSpecialCuts(),     -1, -1, 1 );
      pManager->AddProcess(new G4hMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess(new G4hIonisation(),         -1,  2, 3 );
    }
  }
}
void PhysicsList::ConstructOp(){


  //G4Cerenkov* cerenkovProcess = new G4Cerenkov("Cerenkov");
  G4Scintillation* scintillationProcess = new G4Scintillation("Scintillation");
  G4OpAbsorption* absorptionProcess = new G4OpAbsorption();
  //G4OpRayleigh* rayleighScatteringProcess = new G4OpRayleigh();
  //G4OpMieHG* mieHGScatteringProcess = new G4OpMieHG();
  G4OpBoundaryProcess* boundaryProcess = new G4OpBoundaryProcess();
  /*
  int fVerboseLebel = 1;
  cerenkovProcess->SetVerboseLevel(fVerboseLebel);
  scintillationProcess->SetVerboseLevel(fVerboseLebel);
  absorptionProcess->SetVerboseLevel(fVerboseLebel);
  rayleighScatteringProcess->SetVerboseLevel(fVerboseLebel);
  mieHGScatteringProcess->SetVerboseLevel(fVerboseLebel);
  boundaryProcess->SetVerboseLevel(fVerboseLebel);
  */
  // Use Birks Correction in the Scintillation process
  if(!G4Threading::IsWorkerThread()){
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    //G4Scintillation::AddSaturation(emSaturation);
  }


  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    /*
    if (cerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(cerenkovProcess);
      pmanager->SetProcessOrdering(cerenkovProcess,idxPostStep);
    }
    */
    if (scintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(scintillationProcess);
      pmanager->SetProcessOrderingToLast(scintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(scintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(absorptionProcess);
      //pmanager->AddDiscreteProcess(rayleighScatteringProcess);
      //pmanager->AddDiscreteProcess(mieHGScatteringProcess);
      pmanager->AddDiscreteProcess(boundaryProcess);
    }
  }
}


///////Hadron
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
//#include "G4HadronInelasticProcess.hh"
//#include "G4LElastic.hh"

////Pion
//#include "G4PionPlusInelasticProcess.hh"
//#include "G4LEPionPlusInelastic.hh"
//#include "G4PionMinusInelasticProcess.hh"
//#include "G4LEPionMinusInelastic.hh"
//
////Kaon
//#include "G4KaonPlusInelasticProcess.hh"
//#include "G4LEKaonPlusInelastic.hh"
//#include "G4KaonMinusInelasticProcess.hh"
//#include "G4LEKaonMinusInelastic.hh"
//#include "G4KaonZeroSInelasticProcess.hh"
//#include "G4LEKaonZeroSInelastic.hh"
//#include "G4KaonZeroLInelasticProcess.hh"
//#include "G4LEKaonZeroLInelastic.hh"
//
////Nucleon
//#include "G4ProtonInelasticProcess.hh"
//#include "G4LEProtonInelastic.hh"
//#include "G4NeutronInelasticProcess.hh"
//#include "G4LENeutronInelastic.hh"
// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
//#include "G4GGNuclNuclCrossSection.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"

void PhysicsList:: ConstructHadronic()
{
  /*
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4LElastic* theElasticModel = new G4LElastic;
  theElasticProcess->RegisterMe(theElasticModel);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    //Pion
    if( particleName == "pi+" ) {
      pmanager->AddProcess( theElasticProcess );
      G4PionPlusInelasticProcess* theInelasticProcess
	= new G4PionPlusInelasticProcess( "inelastic" );
      G4LEPionPlusInelastic* theLEInelasticModel
	= new G4LEPionPlusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "pi-" ) {
      pmanager->AddProcess( theElasticProcess );
      G4PionMinusInelasticProcess* theInelasticProcess
	= new G4PionMinusInelasticProcess( "inelastic" );
      G4LEPionMinusInelastic* theLEInelasticModel
	= new G4LEPionMinusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    //Kaon
    else if( particleName == "kaon+" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonPlusInelasticProcess* theInelasticProcess
	= new G4KaonPlusInelasticProcess( "inelastic" );
      G4LEKaonPlusInelastic* theLEInelasticModel
	= new G4LEKaonPlusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon-" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonMinusInelasticProcess* theInelasticProcess
	= new G4KaonMinusInelasticProcess( "inelastic" );
      G4LEKaonMinusInelastic* theLEInelasticModel
	= new G4LEKaonMinusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon0S" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonZeroSInelasticProcess* theInelasticProcess
	= new G4KaonZeroSInelasticProcess( "inelastic" );
      G4LEKaonZeroSInelastic* theLEInelasticModel
	= new G4LEKaonZeroSInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon0L" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonZeroLInelasticProcess* theInelasticProcess
	= new G4KaonZeroLInelasticProcess( "inelastic" );
      G4LEKaonZeroLInelastic* theLEInelasticModel
	= new G4LEKaonZeroLInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    //Nucleon
    else if( particleName == "proton" ) {
      pmanager->AddProcess( theElasticProcess );
      G4ProtonInelasticProcess* theInelasticProcess
	= new G4ProtonInelasticProcess( "inelastic" );
      G4LEProtonInelastic* theLEInelasticModel
	= new G4LEProtonInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "neutron" ) {
      pmanager->AddProcess( theElasticProcess );
      G4NeutronInelasticProcess* theInelasticProcess
	= new G4NeutronInelasticProcess( "inelastic" );
      G4LENeutronInelastic* theLEInelasticModel
	= new G4LENeutronInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
  */
}

///////Decay
#include "G4Decay.hh"

void PhysicsList::ConstructDecay()
{
  G4Decay *theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    if( theDecayProcess->IsApplicable(*particle) ){
      pManager->AddProcess( theDecayProcess );
      pManager->SetProcessOrdering( theDecayProcess, idxPostStep );
      pManager->SetProcessOrdering( theDecayProcess, idxAtRest );
    }
  }
}
