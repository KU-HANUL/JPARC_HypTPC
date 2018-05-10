/*
  SimpleFieldElement.hh
*/

#ifndef SimpleFieldElement_h
#define SimpleFieldElement_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4AffineTransform;
class G4VPhysicalVolume;

class SimpleFieldElement 
{
public:
  SimpleFieldElement( const G4String &name, 
		      const G4ThreeVector &pos,
		      const G4RotationMatrix *rotMtx=0 );
  virtual ~SimpleFieldElement();

  virtual bool ExistMagneticField( void ) const = 0;
  virtual bool ExistElectricField( void ) const = 0;

  virtual G4ThreeVector GetMagneticField( const G4ThreeVector &gPos ) const;
  virtual G4ThreeVector GetElectricField( const G4ThreeVector &gPos ) const;

  const G4String & GetName( void ) const { return elemName_; } 

  virtual G4VPhysicalVolume * Place( G4VPhysicalVolume *pMother ) = 0;
protected:
  G4String elemName_;
  G4ThreeVector gPos_;
  G4RotationMatrix rotMtx_;

  G4AffineTransform *pGtoL, *pLtoG;

public:
  G4ThreeVector GetLocalPosition( const G4ThreeVector &gPos ) const;
  G4ThreeVector GetGlobalPosition( const G4ThreeVector &lPos ) const;
  G4ThreeVector GetLocalDirection( const G4ThreeVector &gDir ) const;
  G4ThreeVector GetGlobalDirection( const G4ThreeVector &lDir ) const;
};

#endif
