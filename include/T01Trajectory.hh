#ifndef T01Trajectory_h
#define T01Trajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"     
#include <vector>
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"   
#include "G4Track.hh"
#include "G4Step.hh"

class G4Polyline;

typedef std::vector<G4VTrajectoryPoint*> T01TrajectoryPointContainer;

class T01Trajectory : public G4VTrajectory
{
public:
  T01Trajectory();
  T01Trajectory(const G4Track* aTrack);
  T01Trajectory(T01Trajectory &);
  virtual ~T01Trajectory();

  inline void* operator new(size_t);
  inline void  operator delete(void*);
  inline int operator == (const T01Trajectory& right) const
  {return (this==&right);} 

  inline G4int GetTrackID() const
  { return fTrackID; }
  inline G4int GetParentID() const
  { return fParentID; }
  inline G4String GetParticleName() const
  { return ParticleName; }
  inline G4double GetCharge() const
  { return PDGCharge; }
  inline G4int GetPDGEncoding() const
  { return PDGEncoding; }
  inline const G4ThreeVector& GetMomentum() const
  { return momentum; }
  inline const G4ThreeVector& GetVertexPosition() const
  { return vertexPosition; }
  inline G4double GetGlobalTime() const
  { return globalTime; }
  virtual int GetPointEntries() const
  { return positionRecord->size(); }
  virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
  { return (*positionRecord)[i]; }

  virtual void ShowTrajectory() const;
  virtual void DrawTrajectory(G4int i_mode=0) const;
  virtual void AppendStep(const G4Step* aStep);
  virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);
  inline G4ThreeVector GetInitialMomentum() const
  { return momentum; }

  G4ParticleDefinition* GetParticleDefinition();

private:
  T01TrajectoryPointContainer* positionRecord;
  G4int                        fTrackID;
  G4int                        fParentID;
  G4ParticleDefinition*        fpParticleDefinition;
  G4String                     ParticleName;
  G4double                     PDGCharge;
  G4int                        PDGEncoding;
  G4ThreeVector                momentum;
  G4ThreeVector                vertexPosition;
  G4double                     globalTime;

};

extern G4Allocator<T01Trajectory> myTrajectoryAllocator;

inline void* T01Trajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void T01Trajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((T01Trajectory*)aTrajectory);
}

#endif
