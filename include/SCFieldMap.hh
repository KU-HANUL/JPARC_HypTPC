/*
  SCFieldMap.hh
*/

#ifndef SCFieldMap_h
#define SCFieldMap_h 1

#include <string>
#include <vector>

class SCFieldMap
{
public:
  SCFieldMap( const char *filename=0, double ScaleFactor=1.0 );
  ~SCFieldMap();

private:
  SCFieldMap( const SCFieldMap & );
  SCFieldMap & operator = ( const SCFieldMap & );

private:
  std::string filename_;
  double ScaleFactor_;

public:
  bool Initialize( void );
  bool GetFieldValue( const double point[3], double *Bfield ) const;
private:
  struct FD {
    float x, y, z;
  };

  typedef std::vector < std::vector < std::vector < FD > > > FDContainer;
  FDContainer B;
  double X0, Y0, Z0, dX, dY, dZ;
  int Nx, Ny, Nz;

  void cleanupMap( void );
};

#endif
