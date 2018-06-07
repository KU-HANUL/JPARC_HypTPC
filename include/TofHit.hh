/*
  TofHit.hh from SKS
  2017/10  Yang
*/

#ifndef TofHit_h
#define TofHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

class TofHit : public G4VHit
{
public:
  TofHit();
  ~TofHit() { }
private:
  TofHit( const TofHit & );
  TofHit & operator = ( const TofHit & );
public:
  int operator == ( const TofHit & ) const { return 0; }

  inline void * operator new ( size_t size );
  inline void operator delete( void * aHit );

  void Draw();
  void Print();

private:
  G4int pass_;
  G4int segID_;
  G4double time_;
  G4double edep_;
  G4ThreeVector pos_;
  G4ThreeVector mom_;
  G4ThreeVector lmom_;
  G4ThreeVector vtxpos_;
  G4int trackNo_;
  G4bool fSignal_;
  G4double xl_, yl_;
  G4double path_;
  G4int HitPartID_;
  G4int parentID1_;
  G4int parentID2_;
  G4int parentID3_;
  G4int parentPID1_;
  G4int parentPID2_;
  G4int parentPID3_;

public:
  void SetPass( void ) { pass_=1; }
  void SetSegmentID( G4int id ) { segID_ = id; }
  void SetParentID1( G4int id ) { parentID1_ = id; }
  void SetParentID2( G4int id ) { parentID2_ = id; }
  void SetParentID3( G4int id ) { parentID3_ = id; }
  void SetParentPID1( G4int id ) { parentPID1_ = id; }
  void SetParentPID2( G4int id ) { parentPID2_ = id; }
  void SetParentPID3( G4int id ) { parentPID3_ = id; }
  void SetTime( G4double time ) { time_=time; }
  void SetEdep( G4double edep ) { edep_ = edep; }
  void AddEdep( G4double edep ) { edep_ += edep; }
  void SetPos( const G4ThreeVector &pos ) { pos_=pos; }
  void SetMom( const G4ThreeVector &mom ) { mom_=mom; }
  void SetLMom( const G4ThreeVector &mom ) { lmom_=mom; }
  void SetPathLength( G4double path ) { path_ = path; }
  void SetVtxPos( const G4ThreeVector &pos ) { vtxpos_=pos; }
  void SetTrackNo( G4int no ) { trackNo_=no; }
  void SetTrueSignal() { fSignal_=true; }
  void SetFalseSignal() { fSignal_=false; }
  void SetLocalPos( G4double x, G4double y ) { xl_=x; yl_=y; }
  void SetHitParticleID( G4int id )
  { HitPartID_ = id; }

  G4int GetPass() const { return pass_; }
  G4int GetSegmentID() const { return segID_; }
  G4int GetParentID1() const { return parentID1_; }
  G4int GetParentID2() const { return parentID2_; }
  G4int GetParentID3() const { return parentID3_; }
  G4int GetParentPID1() const { return parentPID1_; }
  G4int GetParentPID2() const { return parentPID2_; }
  G4int GetParentPID3() const { return parentPID3_; }
  G4double GetTime( void ) const { return time_; }
  G4double GetEdep() const { return edep_; }
  G4ThreeVector GetPos( void ) const { return pos_; }
  G4ThreeVector GetMom( void ) const { return mom_; }
  G4ThreeVector GetLMom( void ) const { return lmom_; }
  G4double GetPathLength( void ) const { return path_; }
  G4ThreeVector GetVtxPos( void ) const { return vtxpos_; }
  G4int GetTrackNo( void ) const { return trackNo_; }
  G4bool IsTrueSignal() const { return fSignal_; }
  G4double GetXLocal( void ) const { return xl_; }
  G4double GetYLocal( void ) const { return yl_; }
  G4int GetHitParticleID( void ) const
  {
    return HitPartID_;
  }
};

typedef G4THitsCollection<TofHit> TofHitsCollection;
extern G4Allocator<TofHit> TofHitAllocator;

inline void * TofHit::operator new( size_t )
{
  return static_cast<void *>( TofHitAllocator.MallocSingle() );
}

inline void TofHit::operator delete( void *aHit )
{
  TofHitAllocator.FreeSingle( static_cast<TofHit *>( aHit ) );
}


#endif
