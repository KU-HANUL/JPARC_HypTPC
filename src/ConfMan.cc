/*
  ConfMan.cc
  2017/8  Yang
*/

#include "ConfMan.hh"
#include "GeomMan.hh"
#include <iostream>
#include <iomanip>
#include <sstream>

ConfMan * ConfMan::confManager_ = 0;
const std::string defFieldMapFile="FieldMap.dat";
const std::string defEvtGenDecayFile="param/EvtGenDecay.dat";
const std::string defEvtGenPDLFile="param/EvtGenPDL.dat";

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),
    GeomManager_(0),
    fField_(false), FieldMapName_(defFieldMapFile), FieldScale_(0),
    fStepping(0),
    EvtGenDecayName_(defEvtGenDecayFile), EvtGenPDLName_(defEvtGenPDLFile),
    ReactionMode_(0),
    BeamMomentumMode_(0),
    bpx_(0), bpy_(0), bpz_(0), bvx_(0), bvy_(0), bvz_(0)
{
  static const std::string funcname="[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constring twice" << std::endl;
    exit(-1);
  }
  confManager_=this;
}

ConfMan::~ConfMan( )
{
  //if(DCGeomManager_) delete DCGeomManager_;
  confManager_=0;
}

const G4int BufSize=300;

bool ConfMan::Initialize( void )
{
  static const std::string funcname="[ConfMan::Initialize]";

  char buf[BufSize], buf1[BufSize];
  G4double val1, val2, val3;
  G4int id;

  FILE *fp;
  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while( fgets( buf, BufSize, fp ) != 0 ){
    if( buf[0]!='#' ){
      // Geometry
      //if( sscanf(buf,"GEO: %s",buf1)==1 ){
      //GeomFileName_=buf1;
      //}
      //Field
      if( sscanf(buf,"FLDMAP: %s",buf1)==1 ){
        FieldMapName_=buf1; fField_=true;
      }

      if( sscanf(buf,"MAPSCALE: %lf", &val1)==1 )
	{
	  FieldScale_=val1;
	}

      // Physics Process
      /*
      else if( sscanf(buf,"EM: %d", &id )==1 ){
        if(id) fPhysProc=FPhysProEMOn(fPhysProc);
        else   fPhysProc=FPhysProEMOff(fPhysProc);
      }
      else if( sscanf(buf,"DECAY: %d", &id )==1 ){
        if(id) fPhysProc=FPhysProDCYOn(fPhysProc);
        else   fPhysProc=FPhysProDCYOff(fPhysProc);
      }
      else if( sscanf(buf,"HADRON: %d", &id )==1 ){
        if(id) fPhysProc=FPhysProHDOn(fPhysProc);
        else   fPhysProc=FPhysProHDOff(fPhysProc);
      }
      */
      //
      // Stepping Action
      else if( sscanf(buf,"SCSTOP: %lf", &val1 )==1 )
        if(val1) fStepping=FStopSCOn(fStepping);
        else     fStepping=FStopSCOff(fStepping);
      else if( sscanf(buf,"SCGSTOP: %lf", &val1 )==1 )
        if(val1) fStepping=FStopSCGamOn(fStepping);
        else     fStepping=FStopSCGamOff(fStepping);
      else if( sscanf(buf,"NUSTOP: %lf", &val1 )==1 )
        if(val1) fStepping=FStopNuOn(fStepping);
        else     fStepping=FStopNuOff(fStepping);
      else if( sscanf(buf,"GSTOP: %lf", &val1 )==1 )
        if(val1) fStepping=FStopGamOn(fStepping);
        else     fStepping=FStopGamOff(fStepping);
      else if( sscanf(buf,"ESTOP: %lf", &val1 )==1 )
        if(val1) fStepping=FStopEOn(fStepping);
        else     fStepping=FStopEOff(fStepping);
      //
      // Evt Gen //
      else if( sscanf(buf,"EVTGENDATA1: %s",buf1)==1 )
	{
	  EvtGenDecayName_=buf1;
	  fEvtData1_=true;
	}
      else if( sscanf(buf,"EVTGENDATA2: %s",buf1)==1 )
	{
	  EvtGenPDLName_=buf1;
	  fEvtData2_=true;
	}
      //
      // Beam Momentum//
      else if( sscanf(buf,"BEAMP: %lf %lf %lf", &val1, &val2, &val3  )==3 )
	{
	  bpx_=val1;
	  bpy_=val2;
	  bpz_=val3;
	}
      // Vertex //
      else if( sscanf(buf,"VERTEX: %lf %lf %lf", &val1, &val2, &val3  )==3 )
	{
	  bvx_=val1;
	  bvy_=val2;
	  bvz_=val3;
	}

      //
      // Reaction //
      else if( sscanf(buf,"REACTION: %d", &id )==1 )
        ReactionMode_=id;
      else if( sscanf(buf,"BEAMMOMENTUMMODE: %d", &id )==1 )
        BeamMomentumMode_=id;
      //
      else {
	std::cout << funcname << ": un-recognized record\n"
		  << buf << std::endl;
      }
    }
  }

  fclose(fp);
  PrintParameters();
  InitializeParameterFiles();

  return true;
}


bool ConfMan::InitializeParameterFiles( void )
{
  //GeomManager_ = & GeomMan::GetInstance();
  //GeomManager_->Initialize(GeomFileName_);
  return true;
}

void ConfMan::PrintParameters( void )
{
  std::cout << "-----------------" << ConfFileName_ << "------" << std::endl;
  std::cout << "************** Geometry **************" << std::endl;
  //std::cout << "Geom. Param.: " << GeomFileName_ << std::endl;
  std::cout << "************ Physics Process *******" << std::endl;
  //std::cout << "EM:     " << GetFPhysProcEM(fPhysProc) << std::endl;
  //std::cout << "DECAY:  " << GetFPhysProcDCY(fPhysProc) << std::endl;
  //std::cout << "HADRON: " << GetFPhysProcHD(fPhysProc) << std::endl;
  std::cout << "************ Field Map *************" << std::endl;
  std::cout << "Map File: " << FieldMapName_ << std::endl;
  std::cout << "Field Scale: " << FieldScale_ << std::endl;
  std::cout << "************ Stepping Action *******" << std::endl;
  std::cout << "SC Stop;           " << GetFStopSC(fStepping) << std::endl;
  std::cout << "SC gamma-ray Stop; " << GetFStopSCGam(fStepping) << std::endl;
  std::cout << "Neutrino Stop;      " << GetFStopNu(fStepping) << std::endl;
  std::cout << "Gamma Stop;         " << GetFStopGam(fStepping) << std::endl;
  std::cout << "Electron Stop;      " << GetFStopE(fStepping) << std::endl;
  std::cout << "************ EventGen *************" << std::endl;
  std::cout << "Evt PDL: " << EvtGenPDLName_ << std::endl;
  std::cout << "Evt Decay: " << EvtGenDecayName_ << std::endl;

  std::cout << "************ Reaction Mode *************" << std::endl;
  std::cout << "Reaction Mode: " << ReactionMode_ << std::endl;

  std::cout << "************ Beam *************" << std::endl;
  std::cout << "Beam Momentum: (" << bpx_ <<","<<bpy_<<","<<bpz_<<")" << "GeV/c" << std::endl;
  std::cout << "Vertex: (" << bvx_ <<","<<bvy_<<","<<bvz_<<")" << "mm" << std::endl;

  std::cout << "--------------------------------------------" << std::endl;
}
