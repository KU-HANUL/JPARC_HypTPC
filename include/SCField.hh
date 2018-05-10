/*
  SCField.hh
*/

#ifndef SCField_h 
#define SCField_h 1

#include "SCFieldMap.hh"
#include "G4MagneticField.hh"

#include <vector>

class SimpleFieldElement;

class SCField : public G4MagneticField
{
public:
  SCField();
  explicit SCField( const std::string &FieldMapName, 
		     double scaleFactor=1.0 );
  ~SCField();

private:
  SCField( const SCField & );
  SCField & operator = ( const SCField & );

public:
  void GetFieldValue( const double Point[4], double *Bfield ) const;
  void cleanupSimpleElementList( void );
  void AddSimpleElement( SimpleFieldElement *elem );

private:
  SCFieldMap fMap;

  typedef std::vector <SimpleFieldElement *> FMContainer;
  typedef std::vector <SimpleFieldElement *>
  ::const_iterator FMIterator;

  FMContainer elemList_;
};

#endif

