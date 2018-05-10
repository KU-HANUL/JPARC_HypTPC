/*
  SCField.cc
*/

#include "SCField.hh"
#include "SimpleFieldElement.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

SCField::SCField( const std::string &FieldMapName, double scaleFactor )
  : fMap( FieldMapName.c_str(),scaleFactor )
{
  fMap.Initialize();
}

SCField::~SCField()
{
}

void SCField::GetFieldValue( const double Point[4], 
			      double *Bfield ) const
{
  double X[3];
  X[0]=Point[0]/mm; X[1]=Point[1]/mm; X[2]=Point[2]/mm;

  if( fMap.GetFieldValue( X, Bfield ) )
    {
      Bfield[0] *= tesla;
      Bfield[1] *= tesla;
      Bfield[2] *= tesla;
    }
  else
    {
      Bfield[0]=Bfield[1]=Bfield[2]=0.0;
    }

  G4ThreeVector gPos( Point[0], Point[1], Point[2] );
  G4ThreeVector B( 0., 0., 0. );
  FMIterator end=elemList_.end();
  for( FMIterator itr=elemList_.begin(); itr!=end; ++itr )
    {
      if( (*itr)->ExistMagneticField() )
	B += (*itr)->GetMagneticField( gPos );
      std::cout<<"FMIterator!!!"<<std::endl;
    }
  
  Bfield[0] += B.x(); Bfield[1] += B.y(); Bfield[2] += B.z();

#if 0
  if(Bfield[1]/tesla > -0.5 && Bfield[1]/tesla < -0.4)
    {
      G4cout << "X=(" << X[0] << "," << X[1] << "," << X[2] << ") "
	     << "B=(" << Bfield[0]/tesla << "," << Bfield[1]/tesla
	     << "," << Bfield[2]/tesla << ")" << G4endl;
    }
#endif
}


void SCField::cleanupSimpleElementList( void )
{
  elemList_.clear();
}

void SCField::AddSimpleElement( SimpleFieldElement *elem )
{
  elemList_.push_back( elem );
}
