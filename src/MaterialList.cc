/*
  MaterialList.cc
  2007/4  K.Shirotori
  2010/4/15 T.Takahashi
*/

#include "MaterialList.hh"

#include "G4String.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

MaterialList::MaterialList()
{
  // Elements
  elH  = new G4Element( "Hydrogen"  , "H"  ,  1.,   1.00794 *g/mole );
  elD  = new G4Element( "Deuteron"  , "D"  ,  1.,   2.01410 *g/mole );
  elHe = new G4Element( "Helium"    , "He" ,  2.,   4.002602*g/mole );
  elLi = new G4Element( "Lithium"   , "Li" ,  3.,   6.941   *g/mole );
  elBe = new G4Element( "Beryllium" , "Be" ,  4.,   9.012182*g/mole );
  elB10= new G4Element( "Boron10"   , "B10",  5.,  10.0     *g/mole );
  elB11= new G4Element( "Boron11"   , "B11",  5.,  11.0     *g/mole );
  elC  = new G4Element( "Carbon"    , "C"  ,  6.,  12.011   *g/mole );
  elN  = new G4Element( "Nitrogen"  , "N"  ,  7.,  14.00674 *g/mole );
  elO  = new G4Element( "Oxygen"    , "O"  ,  8.,  15.9994  *g/mole );
  elNa = new G4Element( "Na"        , "Na" , 11.,  23.0     *g/mole );
  elAl = new G4Element( "Aluminum"  , "Al" , 13.,  26.981539*g/mole );
  elSi = new G4Element( "Silicon"   , "Si" , 14.,  28.0855  *g/mole );
  elP  = new G4Element( "Phoshorus" , "P"  , 15.,  30.973762*g/mole );
  elS  = new G4Element( "Sulfur"    , "S"  , 16.,  32.066   *g/mole );
  elAr = new G4Element( "Argon"     , "Ar" , 18.,  39.948   *g/mole );
  elTi = new G4Element( "Titanium"  , "Ti" , 22.,  47.867   *g/mole );
  elCr = new G4Element( "Chrominum" , "Cr" , 24.,  51.9961  *g/mole );
  elMn = new G4Element( "Manganese" , "Mn" , 25.,  54.93805 *g/mole );
  elFe = new G4Element( "Iron"      , "Fe" , 26.,  55.847   *g/mole );
  elNi = new G4Element( "Nickel"    , "Ni" , 28.,  58.69    *g/mole );
  elCu = new G4Element( "Cupper"    , "Cu" , 29.,  63.546   *g/mole );
  elZn = new G4Element( "Zinc"      , "Zn" , 30.,  65.409   *g/mole );
  elGe = new G4Element( "Germanium" , "Ge" , 32.,  72.61    *g/mole );
  elMo = new G4Element( "Molybdenum", "Mo" , 42.,  95.94    *g/mole );
  elI  = new G4Element( "I"         , "I"  , 53., 127.0     *g/mole );
  elCs = new G4Element( "Cesium"    , "Cs" , 55., 132.9054  *g/mole) ;
  elW  = new G4Element( "Tungstem"  , "W"  , 74., 183.84    *g/mole );
  elPt = new G4Element( "Platinum"  , "Pt" , 78., 195.08    *g/mole );
  elPb = new G4Element( "Lead"      , "Pb" , 82., 207.2     *g/mole );
  elBi = new G4Element( "Bismuth"   , "Bi" , 83., 208.98    *g/mole );

  // Simple Materials, Compounds & Mixtures
  HeGas = new G4Material( "HeGas", 2.,  4.002602*g/mole, 0.1787*mg/cm3 );
  HeLiq = new G4Material( "HeLiq", 2.,  4.002602*g/mole, 0.1249* g/cm3 );
  Li    = new G4Material( "Li",    3.,  6.941   *g/mole, 0.534 * g/cm3 );
  Be    = new G4Material( "Be",    4.,  9.012182*g/mole, 1.85  * g/cm3 );
  B10   = new G4Material( "B10",   5., 10.0     *g/mole, 1.42  * g/cm3 );
  B11   = new G4Material( "B11",   5., 11.0     *g/mole, 2.38  * g/cm3 );
  C     = new G4Material( "C",     6., 12.0     *g/mole, 1.8   * g/cm3 );
  Al    = new G4Material( "Al",   13., 26.981539*g/mole, 2.70  * g/cm3 );
  ArGas = new G4Material( "ArGas",18., 39.948   *g/mole, 1.7834*mg/cm3 );
  Ti    = new G4Material( "Ti",   22., 47.867   *g/mole, 4.54  * g/cm3 );
  Fe    = new G4Material( "Fe",   26., 55.847   *g/mole, 7.87  * g/cm3 );
  Ni    = new G4Material( "Ni",   28., 58.69    *g/mole, 8.902 * g/cm3 );
  Cu    = new G4Material( "Cu",   29., 63.546   *g/mole, 8.96  * g/cm3 );
  Ge    = new G4Material( "Ge",   32., 72.61    *g/mole, 5.323 * g/cm3 );
  W     = new G4Material( "W",    74.,183.84    *g/mole,19.3   * g/cm3 );
  Pt    = new G4Material( "Pt",   78.,195.08    *g/mole,21.45  * g/cm3 );
  Pb    = new G4Material( "Pb",   82.,207.2     *g/mole,11.35  * g/cm3 );

  Vacuum= new G4Material( "Vacuum", universe_mean_density, 2 );
  Vacuum-> AddElement( elN, 70.*perCent );
  Vacuum-> AddElement( elO, 30.*perCent );

  Air = new G4Material( "Air", 1.290*mg/cm3, 2 );
  Air->AddElement( elN, 70.*perCent );
  Air->AddElement( elO, 30.*perCent );

  Water= new G4Material( "Water", 1.*g/cm3, 2 );
  Water-> AddElement( elH, 2 );
  Water-> AddElement( elO, 1 );

  BGO= new G4Material( "BGO", 7.23*g/cm3, 3 );
  BGO-> AddElement( elBi, 4 );
  BGO-> AddElement( elGe, 3 );
  BGO-> AddElement( elO, 12 );

  PWO= new G4Material( "PWO", 8.28*g/cm3, 3 );
  PWO-> AddElement(elPb, 1 );
  PWO-> AddElement(elW,  1 );
  PWO-> AddElement(elO,  4 );

  CsI= new G4Material( "CsI", 4.526*g/cm3, 2 );
  CsI-> AddElement(elCs, 1 );
  CsI-> AddElement(elI,  1 );

  NaI= new G4Material( "NaI", 3.67*g/cm3, 2 );
  NaI-> AddElement(elNa, 1 );
  NaI-> AddElement(elI,  1 );

  SUS316L = new G4Material( "SUS316L", 7.5*g/cm3, 10 );
  SUS316L->AddElement( elFe, 66.337*perCent );
  SUS316L->AddElement( elC,   0.028*perCent );
  SUS316L->AddElement( elO,   0.50 *perCent );
  SUS316L->AddElement( elNi, 12.50 *perCent );
  SUS316L->AddElement( elCr, 17.00 *perCent );
  SUS316L->AddElement( elMo,  2.43 *perCent );
  SUS316L->AddElement( elP,   0.018*perCent );
  SUS316L->AddElement( elS,   0.007*perCent );
  SUS316L->AddElement( elSi,  0.69 *perCent );
  SUS316L->AddElement( elMn,  0.49 *perCent );

  Aerogel = new G4Material( "Aerogel", 0.2000*g/cm3, 2 );
  Aerogel->AddElement( elSi, 1 );
  Aerogel->AddElement( elO,  2 );

  // PolyStylene
  Scin = new G4Material( "Scintillator", 1.032*g/cm3, 2 );
  Scin->AddElement( elC, 2 );
  Scin->AddElement( elH, 6 );

  // Polyethylene
  Polyethylene= new G4Material( "Polyethylene", 0.93*g/cm3, 2 );
  Polyethylene-> AddElement( elC, 1 );
  Polyethylene-> AddElement( elH, 2 );

  // LiO
  LiO = new G4Material( "LithiumOxide", 2.013*g/cm3, 2 );
  LiO->AddElement( elLi,2 );
  LiO->AddElement( elO, 1 );

  // LiN
  LiN = new G4Material( "LithiumNitride", 1.38*g/cm3, 2 );
  LiN->AddElement( elLi,3 );
  LiN->AddElement( elN, 1 );

  // Brass
  Brass = new G4Material("Brass", 8.4*g/cm3, 2);
  Brass->AddElement(elCu, 0.65);
  Brass->AddElement(elZn, 0.35);

  // Liq.H2
  LiqH2 = new G4Material( "Liq-H2", 1.,  1.00794*g/mole, 0.0709* g/cm3, kStateLiquid );

  // Liq.D2
  LiqD2 = new G4Material( "Liq-D2", 1.,  2.01410*g/mole, 0.1419* g/cm3, kStateLiquid );

  // MethaneGas
  MethaneGas = new G4Material("MethaneGas", 0.7162*mg/cm3, 2, kStateGas );
  MethaneGas->AddElement(elH,4);
  MethaneGas->AddElement(elC,1);

  // Etane Gas
  EthaneGas = new G4Material("EThanGas", 1.342*mg/cm3, 2, kStateGas  );
  EthaneGas->AddElement(elH,6);
  EthaneGas->AddElement(elC,2);

  // IsoButane Gas
  IsoButaneGas = new G4Material("IsobutaneGas",2.594*mg/cm3, 2, kStateGas );
  IsoButaneGas->AddElement(elH,10);
  IsoButaneGas->AddElement(elC,4);

  // P10 Gas
  // Ar (90) Methane (10) by volume
  P10Gas = new G4Material("P10Gas", 1.6767*mg/cm3, 2, kStateGas );
  P10Gas->AddMaterial(ArGas,      0.9573);
  P10Gas->AddMaterial(MethaneGas, 0.0427);

  // Ar:Ethane=50:50
  Ar50Ethane50Gas = new G4Material("Ar:Ethane=50:50", 1.5627*mg/cm3, 2, kStateGas );
  Ar50Ethane50Gas->AddMaterial(ArGas,     0.5706);
  Ar50Ethane50Gas->AddMaterial(EthaneGas, 0.4294);

  // Ar::IsoButhane=80:20
  Ar80IsoButane20Gas = new G4Material("Ar:IsoButane=80:20",
				      1.9455*mg/cm3, 2, kStateGas );
  Ar80IsoButane20Gas->AddMaterial(ArGas,        0.7333);
  Ar80IsoButane20Gas->AddMaterial(IsoButaneGas, 0.2667);

  G10 =  new G4Material("NemaG10",
			1.700*g/cm3,
			4
			);
  G10->AddElement(elSi, 1);
  G10->AddElement(elO, 2);
  G10->AddElement(elC, 3);
  G10->AddElement(elH, 3);

  Mylar = new G4Material("Mylar",
			 1.39*g/cm3,
			 3
			 );
  Mylar->AddElement(elO, 2);
  Mylar->AddElement(elC, 5);
  Mylar->AddElement(elH, 4);

  //Saint-Gobain BC404 Scintillator
  BC404 = new G4Material( "BC404", 1.032*g/cm3, 2 );
  BC404 -> AddElement( elC, 474 );
  BC404 -> AddElement( elH, 521 );

  //  G4cout << "MaterialList::Constructor" << G4endl;
}

MaterialList::~MaterialList()
{
  delete Ar80IsoButane20Gas;
  delete Ar50Ethane50Gas;
  delete P10Gas;
  delete IsoButaneGas;
  delete EthaneGas;
  delete MethaneGas;
  delete LiqD2;
  delete LiqH2;
  delete Brass;
  delete LiN;
  delete LiO;
  delete Polyethylene;
  delete Scin;
  delete Aerogel;
  delete SUS316L;
  delete NaI;
  delete CsI;
  delete PWO;
  delete BGO;
  delete Water;
  delete Air;
  delete Vacuum;

  delete Pb;
  delete Pt;
  delete W;
  delete Ge;
  delete Cu;
  delete Ni;
  delete Fe;
  delete Ti;
  delete ArGas;
  delete Al;
  delete C;
  delete B11;
  delete B10;
  delete Be;
  delete Li;
  delete HeLiq;
  delete HeGas;

  delete elBi;
  delete elPb;
  delete elPt;
  delete elW;
  delete elCs;
  delete elI;
  delete elMo;
  delete elZn;
  delete elGe;
  delete elCu;
  delete elNi;
  delete elFe;
  delete elMn;
  delete elCr;
  delete elTi;
  delete elAr;
  delete elS;
  delete elP;
  delete elSi;
  delete elAl;
  delete elNa;
  delete elO;
  delete elN;
  delete elC;
  delete elB11;
  delete elB10;
  delete elBe;
  delete elLi;
  delete elHe;
  delete elD;
  delete elH;
  delete G10;
  delete Mylar;
  delete BC404;
}

G4Material * MaterialList::chooseMaterial( int id )
{
  // see. ConfMan.hh for Sks Material ID

  G4Material *val=0;
  switch(id){
  case  1: val=Vacuum;             break;
  case  2: val=Air;                break;
  case  3: val=HeGas;              break;
  case  4: val=ArGas;              break;
  case  5: val=Ar50Ethane50Gas;    break;
  case  6: val=Ar80IsoButane20Gas; break;
  case  7: val=P10Gas;             break;

  case 11: val=Fe;                 break;
  case 12: val=SUS316L;            break;
  case 13: val=Al;                 break;
  case 14: val=Pb;                 break;
  case 15: val=Brass;              break;

  case 21: val=LiqH2;              break;
  case 22: val=LiqD2;              break;
  case 23: val=HeLiq;              break;
  case 24: val=Li;                 break;
  case 25: val=Be;                 break;
  case 26: val=B10;                break;
  case 27: val=B11;                break;
  case 28: val=C;                  break;
  case 29: val=LiO;                break;
  case 30: val=LiN;                break;
  case 31: val=Water;              break;
  case 32: val=Polyethylene;       break;
  }
  return val;
}
