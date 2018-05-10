/*
  SCFieldMap.cc
  2007/4  K.Shirotori
*/

#include "SCFieldMap.hh"
#include "ConfMan.hh"

#include "globals.hh"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TString.h"


#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

SCFieldMap::SCFieldMap( const char *filename, double ScaleFactor )
  : filename_(filename), ScaleFactor_(ScaleFactor), Nx(0), Ny(0), Nz(0)
{
}

SCFieldMap::~SCFieldMap()
{
  cleanupMap();
}

bool SCFieldMap::Initialize( void )
{
  static const std::string funcname = "SCFieldMap::Initialize";
  //std::ifstream fsin( filename_.c_str() );
  std::string fieldFilename = filename_.c_str();
  G4cout << "Read field map : " << fieldFilename << G4endl;
  TFile* tf = new TFile( fieldFilename.c_str());
  if( tf == NULL )
    {
      G4cerr << "No Such File." << G4endl;
    }
  TTree* trPar = (TTree*)tf->Get("parameter");
  if( trPar->GetEntries() != 1 )
    {
      G4cerr << "Invailied Entries." << G4endl;
    }
  int  nDiv[3];
  float posMin[3];
  float deltaPos[3];
  trPar->SetBranchAddress("ndiv",nDiv);
  trPar->SetBranchAddress("posMin",posMin);
  trPar->GetEntry(0);
  Nx = (int)nDiv[0];
  Ny = (int)nDiv[1];
  Nz = (int)nDiv[2];


  X0 = posMin[0];
  Y0 = posMin[1];
  Z0 = posMin[2];


  
  dX = -1.0*(double)((int)X0*2 / (Nx-1));
  dY = -1.0*(double)((int)Y0*2 / (Ny-1));
  dZ = -1.0*(double)((int)Z0*2 / (Nz-1));

  //std::cout<<"x0: "<<X0<<" "<<dX<<" Ny: "<<Y0<<" Nz: "<<Z0<<std::endl;

  B.resize(Nx);
  for( int ix=0; ix<Nx; ++ix )
    {
      B[ix].resize(Ny);
      for( int iy=0; iy<Ny; ++iy )
	{
	  B[ix][iy].resize(Nz);
	}
    }

  double x,y,z,bx,by,bz;  
  int npoint=0;
  TTree* trMap = (TTree*)tf->Get("field");
  G4double pos[3];
  G4double field[3];
  trMap->SetBranchAddress("pos",pos);
  trMap->SetBranchAddress("field",field);

  std::cout << "Now reading Field Map " <<std::endl;
  if( trMap->GetEntries() != Nx*Ny*Nz )
    {
      G4cerr << "Invailied Entries : map " << G4endl;
    }

  int ievent = 0;
  for( int ix = 0; ix < Nx; ix++)
    {
      for( int iy = 0; iy < Ny; iy++)
	{
	  for( int iz = 0; iz < Nz; iz++)
	    {
	      trMap->GetEntry( ievent );
	      int ixx = int((pos[0]-X0+0.1*dX)/dX);
	      int iyy = int((pos[1]-Y0+0.1*dY)/dY);
	      int izz = int((pos[2]-Z0+0.1*dZ)/dZ);
	      if( ixx>=0 && ixx<Nx && iyy>=0 && iyy<Ny && izz>=0 && izz<Nz )
		{
		  B[ixx][iyy][izz].x = field[0]*ScaleFactor_;
		  B[ixx][iyy][izz].y = field[1]*ScaleFactor_;
		  B[ixx][iyy][izz].z = field[2]*ScaleFactor_;
		  if(ix%100 == 0 && iy%100 == 0 && iz%100 == 0)
		    {
		      //std::cout<<" Bx:"<<pos[0]<<" "<<ixx<<" "<<B[ixx][iyy][izz].x 
		      //       <<" By:"<<pos[1]<<" "<<iyy<<" "<<B[ixx][iyy][izz].y
		      //       <<" Bz:"<<pos[2]<<" "<<izz<<" "<<B[ixx][iyy][izz].z<<std::endl;
		    }
		}
	      if (ievent%50000 == 0)
		{
		  //std::cout << " (-_-)p[Wait]q "<<ievent<<" ixx: "<<ixx;
		  std::cout << " (-_-)p[Wait]q ";
		  fflush( stdout );
		}

	      ievent++;
	    }
	}
    }
  
  
  std::cout << std::endl << "Finished reading Field Map " << std::endl;;
  tf->Close();
  return true;
}

bool SCFieldMap::GetFieldValue( const double point[3], double *Bfield ) const
{
  static const std::string funcname = "SCFieldMap::GetFieldValue";
  double xt=point[0], yt=point[1], zt=point[2];
  
  int ix1, ix2, iy1, iy2, iz1, iz2;
  ix1=int( (xt-X0)/dX );
  iy1=int( (yt-Y0)/dY );
  iz1=int( (zt-Z0)/dZ );

  double wx1, wx2, wy1, wy2, wz1, wz2;
  if( ix1<0 ) { ix1=ix2=0; wx1=1.; wx2=0.; }
  else if( ix1>=Nx-1 ) { ix1=ix2=Nx-1; wx1=1.; wx2=0.; }
  else { ix2=ix1+1; wx1=(X0+dX*ix2-xt)/dX; wx2=1.-wx1; }

  if( iy1<0 ) { iy1=iy2=0; wy1=1.; wy2=0.; }
  else if( iy1>=Ny-1 ) { iy1=iy2=Ny-1; wy1=1.; wy2=0.; }
  else { iy2=iy1+1; wy1=(Y0+dY*iy2-yt)/dY; wy2=1.-wy1; }

  if( iz1<0 ) { iz1=iz2=0; wz1=1.; wz2=0.; }
  else if( iz1>=Nz-1 ) { iz1=iz2=Nz-1; wz1=1.; wz2=0.; }
  else { iz2=iz1+1; wz1=(Z0+dZ*iz2-zt)/dZ; wz2=1.-wz1; }

  double bx1=wx1*wy1*B[ix1][iy1][iz1].x+wx1*wy2*B[ix1][iy2][iz1].x
    +wx2*wy1*B[ix2][iy1][iz1].x+wx2*wy2*B[ix2][iy2][iz1].x;
  double bx2=wx1*wy1*B[ix1][iy1][iz2].x+wx1*wy2*B[ix1][iy2][iz2].x
    +wx2*wy1*B[ix2][iy1][iz2].x+wx2*wy2*B[ix2][iy2][iz2].x;
  double bx=wz1*bx1+wz2*bx2;
  double by1=wx1*wy1*B[ix1][iy1][iz1].y+wx1*wy2*B[ix1][iy2][iz1].y
    +wx2*wy1*B[ix2][iy1][iz1].y+wx2*wy2*B[ix2][iy2][iz1].y;
  double by2=wx1*wy1*B[ix1][iy1][iz2].y+wx1*wy2*B[ix1][iy2][iz2].y
    +wx2*wy1*B[ix2][iy1][iz2].y+wx2*wy2*B[ix2][iy2][iz2].y;
  double by=wz1*by1+wz2*by2;
  double bz1=wx1*wy1*B[ix1][iy1][iz1].z+wx1*wy2*B[ix1][iy2][iz1].z
    +wx2*wy1*B[ix2][iy1][iz1].z+wx2*wy2*B[ix2][iy2][iz1].z;
  double bz2=wx1*wy1*B[ix1][iy1][iz2].z+wx1*wy2*B[ix1][iy2][iz2].z
    +wx2*wy1*B[ix2][iy1][iz2].z+wx2*wy2*B[ix2][iy2][iz2].z;
  double bz=wz1*bz1+wz2*bz2;

  Bfield[0]=bx; Bfield[1]=by; Bfield[2]=bz; 
  return true;
}

void SCFieldMap::cleanupMap( void )
{
  for( int ix=0; ix<Nx; ++ix ){
    for( int iy=0; iy<Ny; ++iy ){
      B[ix][iy].clear();
    }
    B[ix].clear();
  }
  B.clear();
}
