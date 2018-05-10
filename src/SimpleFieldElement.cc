/*
  SimpleFieldElement.cc
*/

#include "SimpleFieldElement.hh"
#include "G4AffineTransform.hh"

SimpleFieldElement::SimpleFieldElement(  const G4String &name,
					 const G4ThreeVector &pos,
					 const G4RotationMatrix *rotMtx )
  :  elemName_(name), gPos_(pos),
     rotMtx_( rotMtx ? (*rotMtx) : G4RotationMatrix() ),
     pGtoL(new G4AffineTransform(rotMtx_,gPos_)),
     pLtoG(new G4AffineTransform(rotMtx_,gPos_))
{
  pGtoL->Invert();
}

SimpleFieldElement::~SimpleFieldElement()
{
  delete pGtoL;
  delete pLtoG;
}

G4ThreeVector SimpleFieldElement::
GetMagneticField( const G4ThreeVector & ) const
{
  return G4ThreeVector(0.,0.,0.);
}

G4ThreeVector SimpleFieldElement::
GetElectricField( const G4ThreeVector & ) const
{
  return G4ThreeVector(0.,0.,0.);
}

G4ThreeVector SimpleFieldElement::
GetLocalPosition( const G4ThreeVector &gPos ) const
{
  return pGtoL->TransformPoint(gPos);
}

G4ThreeVector SimpleFieldElement::
GetGlobalPosition( const G4ThreeVector &lPos ) const
{
  return pLtoG->TransformPoint(lPos);
}

G4ThreeVector SimpleFieldElement::
GetLocalDirection( const G4ThreeVector &gDir ) const
{
  return pGtoL->TransformAxis(gDir);
}

G4ThreeVector SimpleFieldElement::
GetGlobalDirection( const G4ThreeVector &lDir ) const
{
  return pLtoG->TransformAxis(lDir);
}
