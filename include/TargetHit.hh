/*
  TargetHit.hh from SKS
  2017/10  Yang
*/

#ifndef TargetHit_h
#define TargetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

class TargetHit : public G4VHit
{
public:
  TargetHit();
  ~TargetHit() { }
private:
  TargetHit( const TargetHit & );
  TargetHit & operator = ( const TargetHit & );
public:
  int operator == ( const TargetHit & ) const { return 0; }
  
  inline void * operator new ( size_t size );
  inline void operator delete( void * aHit );

  void Draw() ;
  void Print() ;

private:
  G4int pass_;
  G4int layerID_; 
  G4int segID_;
  G4double time_;
  G4double edep_;
  G4double path_;
  G4ThreeVector pos_;
  G4ThreeVector mom_;
  G4ThreeVector lmom_;
  G4int trackNo_;
  G4bool fSignal_;
  G4double xl_, yl_;
  G4int HitPartID_;

public:
  void SetPass( void ) { pass_=1; }
  void SetLayerID( G4int id ) { layerID_=id; }
  void SetSegmentID( G4int id ) { segID_ = id; }
  void SetTime( G4double time ) { time_=time; }
  void SetEdep( G4double edep ) { edep_ = edep; }
  void AddEdep( G4double edep ) { edep_ += edep; }
  void SetPos( const G4ThreeVector &pos ) { pos_=pos; }
  void SetMom( const G4ThreeVector &mom ) { mom_=mom; }
  void SetLMom( const G4ThreeVector &mom ) { lmom_=mom; }
  void SetTrackNo( G4int no ) { trackNo_=no; }
  void SetTrueSignal() { fSignal_=true; }
  void SetFalseSignal() { fSignal_=false; }
  void SetLocalPos( G4double x, G4double y ) { xl_=x; yl_=y; }
  void SetHitParticleID( G4int id ) 
  { HitPartID_ = id; }
  void SetPathLength ( G4double path ) { path_ = path;}

  G4int GetPass() const { return pass_; }
  G4int GetLayerID( void ) const { return layerID_; }
  G4int GetSegmentID() const { return segID_; }
  G4double GetTime( void ) const { return time_; }
  G4double GetEdep() const { return edep_; }
  G4double GetPathLength( void ) const { return path_; }
  G4ThreeVector GetPos( void ) const { return pos_; }
  G4ThreeVector GetMom( void ) const { return mom_; }
  G4ThreeVector GetLMom( void ) const { return lmom_; }
  G4int GetTrackNo( void ) const { return trackNo_; }
  G4bool IsTrueSignal() const { return fSignal_; }
  G4double GetXLocal( void ) const { return xl_; }
  G4double GetYLocal( void ) const { return yl_; }
  G4int GetHitParticleID( void ) const 
  { 
    return HitPartID_;
  }
};

typedef G4THitsCollection<TargetHit> TargetHitsCollection;
extern G4Allocator<TargetHit> TargetHitAllocator;

inline void * TargetHit::operator new( size_t )
{
  return static_cast<void *>( TargetHitAllocator.MallocSingle() );
}

inline void TargetHit::operator delete( void *aHit )
{
  TargetHitAllocator.
    FreeSingle( static_cast<TargetHit *>( aHit ) );
}


#endif
