/*
  DetectorConstruction.hh

  2017/10  Yang
*/

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include <string>

class G4Material;
class MaterialList;
class G4MagneticField;
class G4NistManager;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

private:
  G4VPhysicalVolume* Construct();
  G4NistManager* NistMan;


public:
  //virtual G4VPhysicalVolume *Construct();
  //virtual G4ThreeVector TargetPosition( void ) const = 0;
  //G4double TargetLength( void ) const { return TargetSizeZ(); }
  //virtual G4double TargetAngle( void ) const = 0;
  //virtual G4double TargetSizeX( void ) const = 0;
  //virtual G4double TargetSizeY( void ) const = 0;
  //virtual G4double TargetSizeZ( void ) const = 0;
  virtual G4bool IsVolumeStopper( G4VPhysicalVolume *physVol ) const;
  void MakeSCMagnet(G4VPhysicalVolume *physVol);
  void MakeHypTPC(G4VPhysicalVolume *physVol, G4ThreeVector &pos, G4RotationMatrix &rot);
  void MakeHypTPC2(G4VPhysicalVolume *physVol, G4ThreeVector &pos, G4RotationMatrix &rot);
  void MakeTarget(G4VPhysicalVolume *physVol, G4ThreeVector &pos, G4RotationMatrix &rot);
  void MakeTargetDummy(G4VPhysicalVolume *physVol, G4ThreeVector &pos, G4RotationMatrix &rot);
  void MakeTargetH(G4VPhysicalVolume *physVol, G4ThreeVector &pos, G4RotationMatrix &rot);
  void MakeTOF(G4VPhysicalVolume *physVol, G4ThreeVector &pos, G4RotationMatrix &rot);
  //G4bool IsVolumeStopper (G4VPhysicalVolume *physVol);

protected:
  MaterialList *DefineMaterials( void );
  //virtual G4VPhysicalVolume *ConstructPayload( void ) = 0;
  //virtual void MakeSpecMagnet( G4VPhysicalVolume *pMother,
  //			       G4Material *matGap );
  G4MagneticField * MakeUniformMagField( G4double Bz );
  G4MagneticField * MakeMagFieldFromMap( const std::string &filename,
					 G4double NormFac=1.0 );

  MaterialList *mList_;
  G4Material *mat_p10;
  G4Material *mat_air;
  G4Material *mat_C;
  G4Material *mat_LH2;
  G4Material *mat_Scin;
  G4Material *mat_G10;
  G4Material *mat_Mylar;

  G4MagneticField *field_;
};

#endif
