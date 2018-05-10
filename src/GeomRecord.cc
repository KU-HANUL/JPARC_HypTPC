/*
  GeomRecord.cc

  2017/7  Yang
*/

#include "GeomRecord.hh"

#include <cmath>

const double Deg2Rad = acos(-1.)/180.;

void GeomRecord::calcVectors( void )
{
  double ct1=cos(rotAngle1_*Deg2Rad), st1=sin(rotAngle1_*Deg2Rad);
  double ct2=cos(rotAngle2_*Deg2Rad), st2=sin(rotAngle2_*Deg2Rad);
  double ct0=cos(tiltAngle_*Deg2Rad), st0=sin(tiltAngle_*Deg2Rad);

  dxds_ =  ct0*ct2-st0*ct1*st2;
  dxdt_ = -st0*ct2-ct0*ct1*st2;
  dxdu_ =  st1*st2;

  dyds_ =  ct0*st2+st0*ct1*ct2;
  dydt_ = -st0*st2+ct0*ct1*ct2;
  dydu_ = -st1*ct2;

  dzds_ =  st0*st1;
  dzdt_ =  ct0*st1;
  dzdu_ =  ct1;


  dsdx_ =  ct0*ct2-st0*ct1*st2;
  dsdy_ =  ct0*st2+st0*ct1*ct2;
  dsdz_ =  st0*st1;

  dtdx_ = -st0*ct2-ct0*ct1*st2;
  dtdy_ = -st0*st2+ct0*ct1*ct2;
  dtdz_ =  ct0*st1;

  dudx_ =  st1*st2;
  dudy_ = -st1*ct2;
  dudz_ =  ct1;
}

int GeomRecord::WireNumber( double pos ) const
{
  double dw=((pos-ofs_)/dd_)+w0_;
  int iw=int(dw);
  if((dw-double(iw))>0.5)
    return iw+1;
  else
    return iw;
}
