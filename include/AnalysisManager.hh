#ifndef AnalysisManager_h
#define AnalysisManager_h

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <vector>

class G4Run;
class G4Event;
class TFile;
class TTree;
class EvtVector4R;

const int num_beam = 10;
const int num_evtgen = 100;
const int num_tpchit = 300;
const int num_tofhit = 100;
const int num_targethit = 100;


class AnalysisManager
{
public:
  void BeginOfRun(const G4Run*);
  void EndOfRun(const G4Run*); 
  void BeginOfEvent(const G4Event*); 
  void EndOfEvent(const G4Event*);
  void SetRootFileName(G4String file){ outfile=file; }

  
private:
  G4String outfile;
  TFile *hfile;
  TTree *tree;

public:
  AnalysisManager( const G4String & histname );
  virtual ~AnalysisManager();

  //static AnalysisManager* GetInstance()
  //{
  //  static AnalysisManager instance;
  //  return &instance;
  //}
  //AnalysisManager():outfile("default.root"){}
  //~AnalysisManager(){}
  AnalysisManager(const AnalysisManager&);
  AnalysisManager& operator=(const AnalysisManager&);

  void SetEvtGen(int j, int partnum, int jmotherfirst, int jmotherlast, int jdaugfirst, int jdauglast, int tr_id, EvtVector4R x, EvtVector4R p);
  void SetBeam(int j, G4ThreeVector D, G4ThreeVector P);
  void Terminate( void ) const;
  void SaveFile ( void ) const;

private:
  G4bool fActive_;

  int event;
  //EvtGen
  int nEvt;
  int evtid[num_evtgen];
  int evtpid[num_evtgen];
  int evttrid[num_evtgen];
  int evtfm[num_evtgen];
  int evtlm[num_evtgen];
  int evtfd[num_evtgen];
  int evtld[num_evtgen];
  double evttime[num_evtgen];
  double evtpx[num_evtgen];
  double evtpy[num_evtgen];
  double evtpz[num_evtgen];
  double evtvx[num_evtgen];
  double evtvy[num_evtgen];
  double evtvz[num_evtgen];
  double evtm[num_evtgen];

  //Beam
  int nBeam;
  double bpx[num_beam];
  double bpy[num_beam];
  double bpz[num_beam];
  double bvx[num_beam];
  double bvy[num_beam];
  double bvz[num_beam];

  //Tpc
  int nhTpcPad;
  int tpcpad[num_tpchit];
  double tpcpadedep[num_tpchit];
  double tpcpadtime[num_tpchit];
  double tpcpadposx[num_tpchit];
  double tpcpadposy[num_tpchit];
  double tpcpadposz[num_tpchit];
  double tpcpadmomx[num_tpchit];
  double tpcpadmomy[num_tpchit];
  double tpcpadmomz[num_tpchit];
  int tpcpadtrid[num_tpchit];
  double tpcpadpath[num_tpchit];
  int tpcpadpid[num_tpchit];

  //Tof
  int nhTof;
  int tofseg[num_tofhit];
  double tofedep[num_tofhit];
  double toftime[num_tofhit];
  double tofposx[num_tofhit];
  double tofposy[num_tofhit];
  double tofposz[num_tofhit];
  double tofmomx[num_tofhit];
  double tofmomy[num_tofhit];
  double tofmomz[num_tofhit];
  int toftrid[num_tofhit];
  double tofpath[num_tofhit];
  int tofpid[num_tofhit];
  int tofparentid1[num_tofhit];
  int tofparentid2[num_tofhit];
  int tofparentid3[num_tofhit];
  int tofparentpid1[num_tofhit];
  int tofparentpid2[num_tofhit];
  int tofparentpid3[num_tofhit];

  //Target
  int nhTarget;
  int targetseg[num_targethit];
  double targetedep[num_targethit];
  double targettime[num_targethit];
  double targetposx[num_targethit];
  double targetposy[num_targethit];
  double targetposz[num_targethit];
  double targetmomx[num_targethit];
  double targetmomy[num_targethit];
  double targetmomz[num_targethit];
  int targettrid[num_targethit];
  double targetpath[num_targethit];
  int targetpid[num_targethit];


};

#endif
