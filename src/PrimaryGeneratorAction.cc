/*
  PrimaryGeneratorAction.cc

  2017/10  Yang
*/

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"
#include "ConfMan.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4LorentzVector.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"


#include "Randomize.hh"

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

#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtSecondary.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include "EvtGenBase/EvtStatus.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtRadCorr.hh"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <zlib.h>

#include <TRandom.h>

//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/gzip.hpp>

PrimaryGeneratorAction::PrimaryGeneratorAction( DetectorConstruction *det,
						AnalysisManager *analysisManager,
						EvtGen *evtgen)
  : G4VUserPrimaryGeneratorAction(), det_(det),
    anaMan_(analysisManager),
    evtgen_(evtgen)
{

  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleTable = G4ParticleTable::GetParticleTable();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  //if(BP_) delete BP_;
  DeleteGuns();
}

void PrimaryGeneratorAction::GeneratePrimaries( G4Event *anEvent )
{
  //std::cout<<"Beam generated "<< std::endl;
  G4int event_num = anEvent->GetEventID()+1;
  //G4cout<<"Events: "<< event_num <<G4endl;
  if(event_num% 1000 == 0)
    {
      G4cout<<"Event# "<<event_num<<G4endl;
    }

  ConfMan *confMan = ConfMan::GetConfManager();
  int reactionMode = confMan->ReactionMode();
  int beamMode = confMan->BeamMode();

  int momMode = confMan->BeamMomentumMode();
  int C = confMan->PionCharge();

  //vertex point and beam momentum
  double bvx, bvy, bvz;
  bvx = confMan->GetBeamVX();
  bvy = confMan->GetBeamVY();
  bvz = confMan->GetBeamVZ();
  if(beamMode==1){ // if put beam profile
    /*
      bvx = G4RandGauss::shoot(bvx_,sigmabvx_);
      bvy = G4RandGauss::shoot(bvy_,sigmabvy_);
      bvz = G4RandGauss::shoot(bvz_,sigmabvz_);
    */
    TRandom *eventgen = new TRandom();
    bvx = eventgen->Gaus(bvx_,sigmabvx_);
    bvy = eventgen->Gaus(bvy_,sigmabvy_);
    bvz = eventgen->Gaus(bvz_,sigmabvz_);
  }
  G4ThreeVector D(bvx, bvy, bvz);
  bpx_ = confMan->GetBeamPX();
  bpy_ = confMan->GetBeamPY();
  bpz_ = confMan->GetBeamPZ();
  G4ThreeVector P(bpx_, bpy_, bpz_);
  anaMan_->SetBeam(1, D, P);

  //G4cout<<"Reaction Mode: "<< reactionMode << G4endl;
  switch(reactionMode)
    {
    case 1: GenerateTestKaon(anEvent, D, P); break;
    case 2: GenerateTestPion(anEvent, D, P, C); break;
    case 3: GenerateTest72(anEvent, evtgen_, D, P); break;
    case 4: GenerateTest45(anEvent, evtgen_, D, P, C, momMode ); break;
    }
}

void PrimaryGeneratorAction::DeleteGuns()
{
}

G4ParticleGun * PrimaryGeneratorAction::chooseGun( G4int Pid )
{
  //G4String name=PIDParticleName(Pid);
  //return chooseGun( name );
}

G4ParticleGun * PrimaryGeneratorAction::chooseGun( const G4String & name  )
{
}

void PrimaryGeneratorAction::GenerateTestKaon(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P)
{
  double mass_km = 0.493677;
  //particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaonT1"));
  particleGun -> SetParticleDefinition (particleTable -> FindParticle("kaon-"));
  //std::cout<<"Beam energy GenerateTest 2" << std::endl;
  //G4ThreeVector dir = G4RandomDirection();
  G4ThreeVector beamx ( D.x(), D.y(), D.z());
  G4ThreeVector beamp ( P.x(), P.y(), P.z());
  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_km*mass_km+beamp.mag2()) - mass_km )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
}

void PrimaryGeneratorAction::GenerateTestPion(G4Event* anEvent, G4ThreeVector D, G4ThreeVector P, int C)
{
  double mass_pi = 0.13957061;
  if(C==-1) particleGun -> SetParticleDefinition (particleTable -> FindParticle("pi-"));
  else if(C==1) particleGun -> SetParticleDefinition (particleTable -> FindParticle("pi+"));
  else{
    G4cout<<"### Put a pion charge information ###"<<G4endl;
    return;
  }
  //G4ThreeVector dir = G4RandomDirection();
  G4ThreeVector beamx ( D.x(), D.y(), D.z());
  G4ThreeVector beamp ( P.x(), P.y(), P.z());
  G4ThreeVector beampu =  beamp/beamp.mag();
  G4double energy = (sqrt(mass_pi*mass_pi+beamp.mag2()) - mass_pi )*GeV;
  //std::cout<<"Beam energy: "<< energy /GeV << std::endl;

  particleGun->SetParticleMomentumDirection ( beampu );
  particleGun->SetParticleTime ( 0.0 );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);
}

void PrimaryGeneratorAction::GenerateTest72(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P)
{
  //G4cout<<"Start Generate Test72"<<G4endl;
  /// mother particle momentum //
  double mass_km = 0.493677;
  double mass_proton = 0.938272081;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = P.mag();

  lv_beam.setX(P.x());
  lv_beam.setY(P.y());
  lv_beam.setZ(P.z());
  lv_beam.setE(sqrt(mass_km*mass_km + pbeam*pbeam));

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_proton*mass_proton));

  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* lam1663(0);
  static EvtId LAM1663 = EvtPDL::getId(std::string("Lambda(1663)0"));
  G4LorentzVector LvLam1663;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());
  G4ThreeVector TVx (D.x(), D.y(), D.z());
  //double mass_lam1663 = EvtPDL::getMass(LAM1663);
  double mass_lam1663 = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_lam1663 < 1.115683 + 0.547862 )
    {
      G4cout<<"### Beam momentum is not enough to generate Lam1663 ###"<<G4endl;
      return;
    }
  /*
  G4cout<<"########################### Test  ##############################"<<G4endl;
  G4cout<<"Momentum of K-: "<<pbeam << " GeV/c" <<G4endl;
  G4cout<<"Invariant mass of K + p: "<<lv_particle.m() << " GeV/c2" <<G4endl;
  G4cout<<"Momentum of Lam1663: "<< TVp.mag() << " GeV/c"<<G4endl;
  G4cout<<"Mass of Lam1663: "<< mass_lam1663 << " GeV/c2" <<G4endl;
  G4cout<<"################################################################" << G4endl;
  */
  LvLam1663.setVect(TVp);
  LvLam1663.setE(sqrt(mass_lam1663*mass_lam1663+TVp.mag2()));
  //LvLam1663.setE(sqrt(1.664*1.664+TVp.mag2()));

  EvtVector4R pInit_lam1663( LvLam1663.e(), LvLam1663.vect().x(), LvLam1663.vect().y(), LvLam1663.vect().z() );
  lam1663 = EvtParticleFactory::particleFactory(LAM1663, pInit_lam1663);
  GenerateDecay(anEvent, evtGenerator, lam1663, D);
}

void PrimaryGeneratorAction::GenerateTest45(G4Event* anEvent, EvtGen *evtGenerator, G4ThreeVector D, G4ThreeVector P, int C, int momMode)
{
  //G4cout<<"Start Generate Test45"<<G4endl;
  /// mother particle momentum //
  double mass_pi = 0.13957061;
  double mass_proton = 0.938272081;
  G4LorentzVector lv_beam;
  G4LorentzVector lv_target;
  G4LorentzVector lv_particle;
  double pbeam = P.mag();

  lv_beam.setX(P.x());
  lv_beam.setY(P.y());
  lv_beam.setZ(P.z());
  lv_beam.setE(sqrt(mass_pi*mass_pi + pbeam*pbeam));

  lv_target.setX(0.0);
  lv_target.setY(0.0);
  lv_target.setZ(0.0);
  lv_target.setE(sqrt(mass_proton*mass_proton));

  lv_particle = lv_beam + lv_target;

  // make mother particle //
  EvtParticle* Nstar(0);
  EvtId evtid_N;
  if(C==-1){
    if(momMode==0)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1460)0")); //p=0.635 GeV/c
    else if(momMode==1)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1580)0")); //p=0.835 GeV/c
    else if(momMode==2)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1690)0")); //p=1.035 GeV/c
    else if(momMode==3)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1800)0")); //p=1.235 GeV/c
    else if(momMode==4)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1900)0")); //p=1.435 GeV/c
    else if(momMode==5)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1990)0")); //p=1.635 GeV/c
    else if(momMode==6)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(2090)0")); //p=1.835 GeV/c
    else if(momMode==7)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(2160)0")); //p=2.000 GeV/c

    else if(momMode==8)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1460)0")); //p=0.635 GeV/c
    else if(momMode==9)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1580)0")); //p=0.835 GeV/c
    else if(momMode==10)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1690)0")); //p=1.035 GeV/c
    else if(momMode==11)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1800)0")); //p=1.235 GeV/c
    else if(momMode==12)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1900)0")); //p=1.435 GeV/c
    else if(momMode==13)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1990)0")); //p=1.635 GeV/c
    else if(momMode==14)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(2090)0")); //p=1.835 GeV/c
    else if(momMode==15)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(2160)0")); //p=2.000 GeV/c

    else{
      G4cout<<"### No Particle data in param/EVT ###"<<G4endl;
      return;
    }
  }
  else if(C==1){
    if(momMode==0)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1460)++")); //p=0.635 GeV/c
    else if(momMode==1)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1580)++")); //p=0.835 GeV/c
    else if(momMode==2)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1690)++")); //p=1.035 GeV/c
    else if(momMode==3)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1800)++")); //p=1.235 GeV/c
    else if(momMode==4)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1900)++")); //p=1.435 GeV/c
    else if(momMode==5)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(1990)++")); //p=1.635 GeV/c
    else if(momMode==6)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(2090)++")); //p=1.835 GeV/c
    else if(momMode==7)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipip(2160)++")); //p=2.000 GeV/c

    else if(momMode==8)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1460)++")); //p=0.635 GeV/c
    else if(momMode==9)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1580)++")); //p=0.835 GeV/c
    else if(momMode==10)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1690)++")); //p=1.035 GeV/c
    else if(momMode==11)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1800)++")); //p=1.235 GeV/c
    else if(momMode==12)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1900)++")); //p=1.435 GeV/c
    else if(momMode==13)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(1990)++")); //p=1.635 GeV/c
    else if(momMode==14)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(2090)++")); //p=1.835 GeV/c
    else if(momMode==15)  evtid_N = EvtPDL::getId(std::string("PhaseSpace_pipin(2160)++")); //p=2.000 GeV/c

    else{
      G4cout<<"### No Particle data in param/EVT ###"<<G4endl;
      return;
    }
  }
  else{
    G4cout<<"### Put a pion charge information ###"<<G4endl;
    return;
  }

  G4LorentzVector Lv_N;
  G4ThreeVector TVp (lv_particle.x(), lv_particle.y(), lv_particle.z());
  G4ThreeVector TVx (D.x(), D.y(), D.z());

  double mass_N = sqrt((lv_beam.e()+lv_target.e())*(lv_beam.e()+lv_target.e()) - pbeam*pbeam);

  // check total energy //
  if(mass_N < 0.938272081 + 2*0.13957061 ) //p, pi, pi mass sum
    {
      G4cout<<"### Beam momentum is not enough to generate Particle ###"<<G4endl;
      return;
    }
  /*
  G4cout<<"########################### Test  ##############################"<<G4endl;
  G4cout<<"Momentum of K-: "<<pbeam << " GeV/c" <<G4endl;
  G4cout<<"Invariant mass of K + p: "<<lv_particle.m() << " GeV/c2" <<G4endl;
  G4cout<<"Momentum of generated particle: "<< TVp.mag() << " GeV/c"<<G4endl;
  G4cout<<"Mass of generated particle: "<< mass_N << " GeV/c2" <<G4endl;
  G4cout<<"################################################################" << G4endl;
  */
  Lv_N.setVect(TVp);
  Lv_N.setE(sqrt(mass_N*mass_N+TVp.mag2()));

  EvtVector4R pInit_N( Lv_N.e(), Lv_N.vect().x(), Lv_N.vect().y(), Lv_N.vect().z() );
  Nstar = EvtParticleFactory::particleFactory(evtid_N, pInit_N);
  GenerateDecay(anEvent, evtGenerator, Nstar, D);
}

void PrimaryGeneratorAction::GenerateDecay(G4Event* anEvent, EvtGen *evtGenerator, EvtParticle* particle, G4ThreeVector D)
{


  static EvtStdHep evtstdhep;
  static EvtSecondary evtsecondary;

  EvtId        list_of_stable[10];
  EvtParticle* stable_parent[10];

  list_of_stable[0]=EvtId(-1,-1);
  stable_parent[0]=0;

  evtsecondary.init();
  evtstdhep.init();
  evtGenerator -> generateDecay(particle);
  particle->makeStdHep(evtstdhep,evtsecondary,list_of_stable);

  int npart = evtstdhep.getNPart();

  bool generate_flag = false;
  if(npart<100) generate_flag = true;
  while(!generate_flag)
    {
      G4cout<<"!!!Particles in EvtGen is more than 99."<<G4endl;
      evtsecondary.init();
      evtstdhep.init();
      evtGenerator -> generateDecay(particle);
      particle->makeStdHep(evtstdhep,evtsecondary,list_of_stable);
      if(evtstdhep.getNPart() < 100) generate_flag = true;
    }

  int j;
  int istat;
  int partnum;
  double px,py,pz,e,m;
  double x,y,z,t;

  EvtVector4R p4,x4;
  int n_beam = 0; // number of beams

  for(int i=0;i<evtstdhep.getNPart();i++)
    {
      j=i+1;
      int jmotherfirst=evtstdhep.getFirstMother(i)+1;
      int jmotherlast=evtstdhep.getLastMother(i)+1;
      int jdaugfirst=evtstdhep.getFirstDaughter(i)+1;
      int jdauglast=evtstdhep.getLastDaughter(i)+1;

      partnum=evtstdhep.getStdHepID(i);

      istat=evtstdhep.getIStat(i);

      p4=evtstdhep.getP4(i);
      x4=evtstdhep.getX4(i);

      px=p4.get(1);
      py=p4.get(2);
      pz=p4.get(3);
      e=p4.get(0);

      x=x4.get(1)+D.x();
      y=x4.get(2)+D.y();
      z=x4.get(3)+D.z();
      //t=x4.get(0)+D.t();
      t=x4.get(0)*0.001/2.999792458 *1.e1; //mm/c --> ns
      m=p4.mass();

      EvtVector4R evx4, evp4;
      evx4.set(1, x);
      evx4.set(2, y);
      evx4.set(3, z);
      evx4.set(0, t);

      evp4.set(1, px);
      evp4.set(2, py);
      evp4.set(3, pz);
      evp4.set(0, e);

      bool beam_flag = false;
      int tr_id = -1;
      //if(jdaugfirst ==0 && jdauglast ==0 && (partnum == 2212 || partnum == -211) && jmotherfirst == 2)
      if(jdaugfirst ==0 && jdauglast ==0)
	{
	  makeGun(anEvent, partnum, evx4, evp4);
	  beam_flag = true;
	  n_beam++;
	}
      if(beam_flag == true)
	{
	  tr_id = n_beam;
	}
      anaMan_->SetEvtGen(j, partnum, jmotherfirst, jmotherlast, jdaugfirst, jdauglast, tr_id, evx4, evp4);
    }

#if 0
  G4cout<<"############# Particle decay table  ##############"<<G4endl;
  G4cout<<"Npart: "<<npart<<G4endl;
  for(int i=0;i<evtstdhep.getNPart();i++)
    {
      j=i+1;
      int jmotherfirst=evtstdhep.getFirstMother(i)+1;
      int jmotherlast=evtstdhep.getLastMother(i)+1;
      int jdaugfirst=evtstdhep.getFirstDaughter(i)+1;
      int jdauglast=evtstdhep.getLastDaughter(i)+1;

      p4=evtstdhep.getP4(i);
      x4=evtstdhep.getX4(i);

      px=p4.get(1);
      py=p4.get(2);
      pz=p4.get(3);
      e=p4.get(0);

      x=x4.get(1)+D.x();
      y=x4.get(2)+D.y();
      z=x4.get(3)+D.z();
      t=x4.get(0);
      m=p4.mass();
      G4cout<<"x : "<<x<<" y " <<y<<" z: "<<z<<G4endl;

      partnum=evtstdhep.getStdHepID(i);
      G4cout<<"ID: " << j<< "  Particle Num: "<<partnum<<"  mf: "<<jmotherfirst<< "  ml: "<<jmotherlast << "  df: "<<jdaugfirst << "  dl: "<<jdauglast<<G4endl;
      G4cout<< "   p: "<<(float)sqrt(px*px+py*py+pz*pz) << " e: " << (float)e << " t: "<< (float)t<< " m: "<< (float)m <<G4endl;
    }

  G4cout<<"##################################################"<<G4endl;
#endif
  particle->deleteTree();

}

void PrimaryGeneratorAction::makeGun(G4Event* anEvent, int partnum, EvtVector4R x4, EvtVector4R p4)
{
  particleGun -> SetParticleDefinition (particleTable -> FindParticle(partnum));
  G4ThreeVector beamx ( x4.get(1), x4.get(2), x4.get(3));
  G4ThreeVector beamp ( p4.get(1), p4.get(2), p4.get(3));
  G4ThreeVector beampu = beamp/beamp.mag();
  G4double energy = (p4.get(0) - p4.mass())*GeV;
  //G4double energy = beamp.mag()*GeV;
  double time = x4.get(0);
  G4double g4time =  time * ns;

  //std::cout<<"Beam energy: "<< energy << std::endl;
  //G4cout<<"###############  beam generated  #################" << G4endl;
  //G4cout<<"particle: "<<partnum<< " vertex(x,y,z) ("<< beamx.x() << ", " << beamx.y() << ", " << beamx.z() << ") "
  //	<< "energy: " <<energy/GeV<< " momentum(x,y,z) ("<< p4.get(1) << ", " << p4.get(2) << ", " <<p4.get(3) << ") "<< beamp.mag()<< " time: "<<time
  //	<< G4endl;
  //G4cout<<"##################################################" << G4endl;
  particleGun->SetParticleMomentumDirection ( beampu );
  //particleGun->SetParticleMomentum ( beamp );
  particleGun->SetParticleTime ( g4time );
  particleGun->SetParticlePosition( beamx );
  particleGun->SetParticleEnergy( energy );
  particleGun->GeneratePrimaryVertex( anEvent);

}
