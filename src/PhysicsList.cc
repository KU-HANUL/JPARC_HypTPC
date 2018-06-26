/*
  PhysicsList.cc

  2017/9  Yang
*/

#include "globals.hh"
#include "PhysicsList.hh"
#include "Transportation.hh"

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
#include <iomanip>

#include "G4SystemOfUnits.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4FastSimulationManagerProcess.hh"

#include "ConfMan.hh"

PhysicsList::PhysicsList()
  : G4VUserPhysicsList()
{
  defaultCutValue = 0.1*cm;
}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructHeavyIon();
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

void PhysicsList::ConstructBaryons()
{
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

void PhysicsList::ConstructHeavyIon()
{
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


void PhysicsList::ConstructProcess()
{
  ConfMan *confMan = ConfMan::GetConfManager();
  //int flag = confMan->PhysFlag();

  AddTransportationSks();
  ConstructEM();
  ConstructHadronic();
  ConstructDecay();

  //if( GetFPhysProcEM(flag) ){
  //  ConstructEM();
  //  if( GetFPhysProcHD(flag) ) ConstructHadronic();
  // }
  //if( GetFPhysProcDCY(flag)  ) ConstructDecay();

  // G4StepLimiter *stepLimiter = new G4StepLimiter();
  // G4UserSpecialCuts *userCuts = new G4UserSpecialCuts();
  // theParticleIterator->reset();
  // while ((*theParticleIterator)()){
  //   G4ParticleDefinition* particle = theParticleIterator->value();
  //   G4ProcessManager* pmanager = particle->GetProcessManager();
  //   G4String particleName = particle->GetParticleName();
  //   pmanager ->AddDiscreteProcess(stepLimiter);
  //   pmanager ->AddDiscreteProcess(userCuts);
  // }
}

/*
void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructDecay();
}
*/

///////Transportation
void PhysicsList::AddTransportationSks()
{
  Transportation* theTransportationProcess= new Transportation();

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if ( !particle->IsShortLived() ) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        //G4Exception("PhysicsList::AddTransportation : no process manager!");
      }
      else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager->AddProcess(theTransportationProcess);
        pmanager->SetProcessOrderingToFirst(theTransportationProcess,
                                            idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess,
                                             idxPostStep);
      }
    }
    else {
      // shortlived particle case
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

//#include "G4MultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4UserSpecialCuts.hh"

void PhysicsList::ConstructEM()
{
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
      pManager->AddProcess( new G4UserSpecialCuts(),      -1, -1, 1 );
      pManager->AddProcess( new G4MuMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4MuIonisation(),         -1,  2, 3 );
      pManager->AddProcess( new G4MuBremsstrahlung(),     -1, -1, 4 );
      pManager->AddProcess( new G4MuPairProduction(),     -1, -1, 5 );
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

///////Hadron
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4HadronElastic.hh"

//#include "G4LElastic.hh"

//Pion
//#include "G4LEPionPlusInelastic.hh"
//#include "G4LEPionMinusInelastic.hh"

//Kaon
//#include "G4LEKaonPlusInelastic.hh"
//#include "G4LEKaonMinusInelastic.hh"
//#include "G4LEKaonZeroSInelastic.hh"
//#include "G4LEKaonZeroLInelastic.hh"

//Nucleon
//#include "G4LEProtonInelastic.hh"
//#include "G4LENeutronInelastic.hh"

//Inelastic
#include "G4HadronInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
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

//Cross sections
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
#include "G4ComponentGGNuclNuclXsc.hh"

//#include "G4GGNuclNuclCrossSection.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

//Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

//Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"

//version 4.10.4
void PhysicsList:: ConstructHadronic()
{
  //Elastic models
  //const G4double elastic_elimitPi = 1.0*GeV;
  const G4double elastic_elimitPi = 2.5*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE();
  elastic_he->SetMinEnergy( elastic_elimitPi );


  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4ComponentGGNuclNuclXsc * ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet * theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);

  //auto particleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while ((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionPlusInelasticProcess* theInelasticProcess =
	    new G4PionPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "pi-")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionMinusInelasticProcess* theInelasticProcess =
	    new G4PionMinusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  //Absorption
	  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
	}
      /*
      else if (particleName == "kaon+")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4KaonPlusInelasticProcess* theInelasticProcess =
	    new G4KaonPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon0S")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4KaonZeroSInelasticProcess* theInelasticProcess =
	    new G4KaonZeroSInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon0L")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4KaonZeroLInelasticProcess* theInelasticProcess =
	    new G4KaonZeroLInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon-")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4KaonMinusInelasticProcess* theInelasticProcess =
	    new G4KaonMinusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
	}
      */
      else if (particleName == "proton")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
					GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4ProtonInelasticProcess* theInelasticProcess =
	    new G4ProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "anti_proton")
	{
	  // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
          G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4AntiProtonInelasticProcess* theInelasticProcess =
	    new G4AntiProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  // Absorption
	  pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
	}
      /*
      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy( 19.0*MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
	G4ParticleHPElastic * theElasticNeutronHP = new G4ParticleHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
	theElasticProcess->RegisterMe( theElasticNeutronHP );
	theElasticProcess->AddDataSet( new G4ParticleHPElasticData );
	pmanager->AddDiscreteProcess( theElasticProcess );
	// inelastic scattering
	G4NeutronInelasticProcess* theInelasticProcess =
	  new G4NeutronInelasticProcess("inelastic");
	theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
	theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
	G4ParticleHPInelastic * theNeutronInelasticHPModel = new G4ParticleHPInelastic;
        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
	theInelasticProcess->AddDataSet( new G4ParticleHPInelasticData );
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4HadronCaptureProcess* theCaptureProcess =
	  new G4HadronCaptureProcess;
	G4ParticleHPCapture * theLENeutronCaptureModel = new G4ParticleHPCapture;
	theLENeutronCaptureModel->SetMinEnergy(theHPMin);
	theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
	theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
	theCaptureProcess->AddDataSet( new G4ParticleHPCaptureData);
	pmanager->AddDiscreteProcess(theCaptureProcess);

      }
      else if (particleName == "anti_neutron")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering (include annihilation on-fly)
	  G4AntiNeutronInelasticProcess* theInelasticProcess =
	    new G4AntiNeutronInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "deuteron")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4DeuteronInelasticProcess* theInelasticProcess =
	    new G4DeuteronInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "triton")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4TritonInelasticProcess* theInelasticProcess =
	    new G4TritonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "alpha")
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4AlphaInelasticProcess* theInelasticProcess =
	    new G4AlphaInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      */
    }
}




//version 4.9.6
/*
void PhysicsList:: ConstructHadronic()
{
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  //G4LElastic* theElasticModel = new G4LElastic;
  G4HadronElastic* theElasticModel = new G4HadronElastic;
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
}
*/

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
