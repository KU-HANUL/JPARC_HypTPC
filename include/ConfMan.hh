/*
  ConfMan.hh

  2017/9  Yang
*/

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>

class GeomMan;

inline int FPhysProEMOn( int flag ){ return (flag | 0x01); }
inline int FPhysProDCYOn( int flag ){ return (flag | 0x02); }
inline int FPhysProHDOn( int flag ){ return (flag | 0x04); }

inline int FPhysProEMOff( int flag ){ return (flag & ~0x01); }
inline int FPhysProDCYOff( int flag ){ return (flag & ~0x02); }
inline int FPhysProHDOff( int flag ){ return (flag & ~0x04); }

inline bool GetFPhysProcEM( int flag ){ return (flag & 0x1); }
inline bool GetFPhysProcDCY( int flag ){ return (flag & 0x2); }
inline bool GetFPhysProcHD( int flag ){ return (flag & 0x4); }

inline bool GetFStopSC( int flag ){ return (flag & 0x01); }
inline bool GetFStopSCGam( int flag ){ return (flag & 0x02); }
inline bool GetFStopNu( int flag ){ return (flag & 0x04); }
inline bool GetFStopGam( int flag ){ return (flag & 0x08); }
inline bool GetFStopE( int flag ){ return (flag & 0x10); }

inline int FStopSCOn( int flag ){ return (flag | 0x01); }
inline int FStopSCGamOn( int flag ){ return (flag | 0x02); }
inline int FStopNuOn( int flag ){ return (flag | 0x04); }
inline int FStopGamOn( int flag ){ return (flag | 0x08); }
inline int FStopEOn( int flag ){ return (flag | 0x10); }

inline int FStopSCOff( int flag ){ return (flag & ~0x01); }
inline int FStopSCGamOff( int flag ){ return (flag & ~0x02); }
inline int FStopNuOff( int flag ){ return (flag & ~0x04); }
inline int FStopGamOff( int flag ){ return (flag & ~0x08); }
inline int FStopEOff( int flag ){ return (flag & ~0x10); }


class ConfMan
{
public:
  ConfMan( const std::string & filename );
  ~ConfMan();
private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );
public:
  static ConfMan *GetConfManager( void ) { return confManager_; }
  bool Initialize( void );
  void SetFileName( const std::string & filename ) { ConfFileName_=filename; }

  // Geometry
  GeomMan *GetGeomManager( void ) { return GeomManager_; }

  // Field
  bool ExistField( void ) const { return fField_; }
  const std::string & FieldMapName( void ) const { return FieldMapName_; }
  double GetFieldScale( void ) const { return FieldScale_; }

  // Stepping Action
  int StepFlag( void ) const { return fStepping; }
  bool DoesStopInSC( void ) const { return GetFStopSC(fStepping); }
  bool DoesGamStopInSC( void ) const { return GetFStopSCGam(fStepping); }
  bool DoesNuStop( void ) const { return GetFStopNu(fStepping); }
  bool DoesGamStop( void ) const { return GetFStopGam(fStepping); }
  bool DoesEStop( void ) const { return GetFStopE(fStepping); }

  // Physics Process
  /*
  int PhysFlag( void ) const { return fPhysProc; }
  bool ExistEMProc( void ) const { return GetFPhysProcEM(fPhysProc); }
  bool ExistDCYProc( void ) const { return GetFPhysProcDCY(fPhysProc); }
  bool ExistHDProc( void ) const { return GetFPhysProcHD(fPhysProc); }
  */

  // Evt Gen //
  bool ExistEvtData1( void ) const { return fEvtData1_; }
  bool ExistEvtData2( void ) const { return fEvtData2_; }
  const std::string & EvtGenDecayName( void ) const { return EvtGenDecayName_; }
  const std::string & EvtGenPDLName( void ) const { return EvtGenPDLName_; }

  // Reaction mode //
  int ReactionMode( void ) const { return ReactionMode_; }
  int BeamMomentumMode( void ) const { return BeamMomentumMode_; }
  int PionCharge( void ) const { return PionCharge_ }

  // Beam //
  double GetBeamPX( void ) const { return bpx_; }
  double GetBeamPY( void ) const { return bpy_; }
  double GetBeamPZ( void ) const { return bpz_; }
  double GetBeamVX( void ) const { return bvx_; }
  double GetBeamVY( void ) const { return bvy_; }
  double GetBeamVZ( void ) const { return bvz_; }

private:
  // ConfMan
  static ConfMan *confManager_;
  std::string ConfFileName_;

  // Field
  bool fField_;
  std::string FieldMapName_;
  double FieldScale_;

  // Geometry
  std::string GeomFileName_;
  GeomMan *GeomManager_;

  // Steping Action
  int fStepping;

  // EvtGen //
  bool fEvtData1_, fEvtData2_;
  std::string EvtGenDecayName_;
  std::string EvtGenPDLName_;

  // Reaction mode //
  int ReactionMode_;
  int BeamMomentumMode_;
  int PionCharge_;

  // Beam //
  double bpx_;
  double bpy_;
  double bpz_;
  double bvx_;
  double bvy_;
  double bvz_;

private:
  bool InitializeParameterFiles( void );
  void PrintParameters( void );
};


#endif
