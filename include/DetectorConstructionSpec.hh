/*
  DetectorConstructionSpec.hh

  2012/5  K.Shirotori
*/

#ifndef DetectorConstructionSpec_h
#define DetectorConstructionSpec_h 1

#include "DetectorConstruction.hh"
#include "G4ThreeVector.hh"

class G4Material;
class G4VPhysicalVolume;
class DetectorConstruction;
class G4Material;

class DetectorConstructionSpec : public DetectorConstruction
{
public:
  DetectorConstructionSpec();
  ~DetectorConstructionSpec();

private:
  DetectorConstructionSpec( const DetectorConstructionSpec & );
  DetectorConstructionSpec &
  operator = ( const DetectorConstructionSpec & );

public:
  virtual G4ThreeVector TargetPosition( void ) const;
  virtual G4double TargetAngle( void ) const;
  virtual G4double TargetSizeX( void ) const;
  virtual G4double TargetSizeY( void ) const;
  virtual G4double TargetSizeZ( void ) const;
  virtual G4bool IsVolumeStopper( G4VPhysicalVolume *physVol ) const;

private:
  virtual G4VPhysicalVolume *ConstructPayload( void );
  void MakeTarget( G4VPhysicalVolume *pMother );
  void MakeLiqTarget( G4VPhysicalVolume *pMother );
  void MakeTrackers( G4VPhysicalVolume *pMother );
  void MakeTofCounters( G4VPhysicalVolume *pMother );
  void MakeVDetector( G4VPhysicalVolume *pMother );

  void SetRealMaterials( void );
  void PrintRealMaterialName( void ) const; 
  
  G4Material *matWorld, *matSpecGap;
  G4Material *matBT1Base, *matBT1Layer, *matBT1Box;
  G4Material *matBT2Base, *matBT2Layer, *matBT2Box;

  G4Material *matST1Base, *matST1Layer, *matST1Box;
  G4Material *matST2Base, *matST2Layer, *matST2Box;
  G4Material *matST3Gas,  *matST3Frame, *matST3Box;
  G4Material *matST4Gas,  *matST4Frame, *matST4Box;

  G4Material *matT0Scin, *matT0LG, *matT0PMT, *matT0Box;
  G4Material *matTofScin, *matTofLG, *matTofPMT, *matTofPMTBox,
    *matTofFrame, *matTofBox;
  G4Material *matVD;
  G4Material *matTarget;

  G4double targSizX, targSizY, targSizZ, targAng2; 
  G4ThreeVector targPos;
};

#endif


















