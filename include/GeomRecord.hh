/*
  GeomRecord.hh

  2017/9  Yang
*/

#ifndef GeomRecord_h
#define GeomRecord_h 1

#include <string>
#include <functional>

#include "G4ThreeVector.hh"

class GeomRecord
{
public:
  GeomRecord( int id, const char *name,
                double x, double y, double z, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
  { calcVectors(); }

  GeomRecord( int id, const std::string &name,
                double x, double y, double z, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
  { calcVectors(); }

  GeomRecord( int id, const char *name,
                const G4ThreeVector pos, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
  { calcVectors(); }

  GeomRecord( int id, const std::string &name,
                const G4ThreeVector pos, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
  { calcVectors(); }

  ~GeomRecord() {}
  GeomRecord( const GeomRecord & );
  GeomRecord & operator=( const GeomRecord );
public:
  const G4ThreeVector & Position( void ) const { return pos_; }
  G4ThreeVector NormalVector( void ) const 
  { return G4ThreeVector( dxdu_, dydu_, dzdu_ ); }
  G4ThreeVector UnitVector( void ) const
  { return G4ThreeVector( dxds_, dyds_, dzds_ ); }

  double dsdx( void ) const { return dsdx_; }
  double dsdy( void ) const { return dsdy_; }
  double dsdz( void ) const { return dsdz_; }
  double dtdx( void ) const { return dtdx_; }
  double dtdy( void ) const { return dtdy_; }
  double dtdz( void ) const { return dtdz_; }
  double dudx( void ) const { return dudx_; }
  double dudy( void ) const { return dudy_; }
  double dudz( void ) const { return dudz_; }

  double dxds( void ) const { return dxds_; }
  double dxdt( void ) const { return dxdt_; }
  double dxdu( void ) const { return dxdu_; }
  double dyds( void ) const { return dyds_; }
  double dydt( void ) const { return dydt_; }
  double dydu( void ) const { return dydu_; }
  double dzds( void ) const { return dzds_; }
  double dzdt( void ) const { return dzdt_; }
  double dzdu( void ) const { return dzdu_; }

  double WirePos( int wire ) const { return dd_*(double(wire)-w0_)+ofs_; }
  int WireNumber( double pos ) const; 

private:
  void calcVectors( void );

private:
  int id_;
  std::string name_;
  G4ThreeVector pos_;
  double tiltAngle_, rotAngle1_, rotAngle2_;
  double length_;
  double resol_;
  double w0_, dd_, ofs_;

  double dxds_, dxdt_, dxdu_;
  double dyds_, dydt_, dydu_;
  double dzds_, dzdt_, dzdu_;

  double dsdx_, dsdy_, dsdz_;
  double dtdx_, dtdy_, dtdz_;
  double dudx_, dudy_, dudz_;
 
  friend class GeomMan;
  friend class GeomRecordComp;
};

struct GeomRecordComp 
  : public std::binary_function <GeomRecord *, GeomRecord *, bool> 
{
  bool operator()( const GeomRecord * const p1,
		   const GeomRecord * const p2 ) const
  { return p1->id_ < p2->id_; }
};

#endif
