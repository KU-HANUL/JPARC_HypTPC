/*
  MaterialList.hh
  2007/4  K.Shirotori
  2010/4/15 T.Takahashi

*/

#ifndef MaterialList_h
#define MaterialList_h 1

class G4Element;
class G4Material;

struct MaterialList
{
public:
  MaterialList();
  ~MaterialList();

  G4Material *chooseMaterial( int id );

  G4Element *elH;
  G4Element *elD;
  G4Element *elHe;
  G4Element *elLi;
  G4Element *elBe;
  G4Element *elB10;
  G4Element *elB11;
  G4Element *elC;
  G4Element *elN;
  G4Element *elO;
  G4Element *elNa;
  G4Element *elAl;
  G4Element *elSi;
  G4Element *elP;
  G4Element *elS;
  G4Element *elAr;
  G4Element *elTi;
  G4Element *elCr;
  G4Element *elMn;
  G4Element *elFe;
  G4Element *elNi;
  G4Element *elCu;
  G4Element *elZn;
  G4Element *elGe;
  G4Element *elMo;
  G4Element *elI;
  G4Element *elCs;
  G4Element *elW;
  G4Element *elPt;
  G4Element *elPb;
  G4Element *elBi;

  G4Material *HeGas;
  G4Material *HeLiq;
  G4Material *Li;
  G4Material *Be;
  G4Material *B10;
  G4Material *B11;
  G4Material *C;
  G4Material *Al;
  G4Material *ArGas;
  G4Material *Ti;
  G4Material *Fe;
  G4Material *Ni;
  G4Material *Cu;
  G4Material *Ge;
  G4Material *W;
  G4Material *Pt;
  G4Material *Pb;

  G4Material *Vacuum;
  G4Material *Air;
  G4Material *Water;
  G4Material *BGO;
  G4Material *PWO;
  G4Material *CsI;
  G4Material *NaI;
  G4Material *SUS316L;
  G4Material *Aerogel;
  G4Material *Scin;
  G4Material *Polyethylene;
  G4Material *LiO;
  G4Material *LiN;
  G4Material *Brass;

  G4Material *LiqH2;
  G4Material *LiqD2;

  G4Material *MethaneGas;
  G4Material *EthaneGas;
  G4Material *IsoButaneGas;

  G4Material *P10Gas;
  G4Material *Ar50Ethane50Gas;
  G4Material *Ar80IsoButane20Gas;

  G4Material *G10;
  G4Material *Mylar;
  G4Material *BC404;

private:
  MaterialList( const MaterialList & );
  MaterialList & operator=(const MaterialList & );

};

#endif
