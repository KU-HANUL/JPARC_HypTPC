/*
  TpcHit.hh from SKS
  2017/10  Yang
*/

#ifndef TpcHit_h
#define TpcHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

class TpcHit : public G4VHit
{
public:
  TpcHit();
  ~TpcHit() { }
private:
  TpcHit( const TpcHit & );
  TpcHit & operator = ( const TpcHit & );
public:
  int operator == ( const TpcHit & ) const { return 0; }

  inline void * operator new ( size_t size );
  inline void operator delete( void * aHit );

  void Draw();
  void Print();

private:
  G4int pass_;
  G4int padID_;
  G4double time_;
  G4double edep_;
  G4ThreeVector pos_;
  G4ThreeVector mom_;
  G4ThreeVector lmom_;
  G4int trackNo_;
  G4bool fSignal_;
  G4double xl_, yl_;
  G4double path_;
  G4int particleID_;

public:
  void SetPass( void ) { pass_=1; }
  void SetPadID( G4int id ) { padID_ = id; }
  void SetTime( G4double time ) { time_=time; }
  void SetEdep( G4double edep ) { edep_ = edep; }
  void AddEdep( G4double edep ) { edep_ += edep; }
  void SetPos( const G4ThreeVector &pos ) { pos_=pos; }
  void SetMom( const G4ThreeVector &mom ) { mom_=mom; }
  void SetLMom( const G4ThreeVector &mom ) { lmom_=mom; }
  void SetPathLength( G4double path ) { path_ = path; }
  void SetTrackNo( G4int no ) { trackNo_=no; }
  void SetTrueSignal() { fSignal_=true; }
  void SetFalseSignal() { fSignal_=false; }
  void SetLocalPos( G4double x, G4double y ) { xl_=x; yl_=y; }
  void SetHitParticleID( G4int id ) 
  { particleID_=id; }

  G4int GetPass() const { return pass_; }
  G4int GetPadID() const { return padID_; }
  G4double GetTime( void ) const { return time_; }
  G4double GetEdep() const { return edep_; }
  G4ThreeVector GetPos( void ) const { return pos_; }
  G4ThreeVector GetMom( void ) const { return mom_; }
  G4ThreeVector GetLMom( void ) const { return lmom_; }
  G4double GetPathLength( void ) const { return path_; }
  G4int GetTrackNo( void ) const { return trackNo_; }
  G4bool IsTrueSignal() const { return fSignal_; }
  G4double GetXLocal( void ) const { return xl_; }
  G4double GetYLocal( void ) const { return yl_; }
  G4int GetHitParticleID( void ) const 
  { 
    return particleID_;
  }
};

typedef G4THitsCollection<TpcHit> TpcHitsCollection;
extern G4Allocator<TpcHit> TpcHitAllocator;

inline void * TpcHit::operator new( size_t )
{
  return static_cast<void *>( TpcHitAllocator.MallocSingle() );
}

inline void TpcHit::operator delete( void *aHit )
{
  TpcHitAllocator.FreeSingle( static_cast<TpcHit *>( aHit ) );
}
		       

#endif
