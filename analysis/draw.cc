#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <TPad.h>
#include "e45.h"

using namespace std;

void draw()
{
  //gStyle -> SetOptFit(1000);

  bool pdf_saving = true;
  //bool pdf_saving = false;

  e45 e45_class;
  e45_class.beam_through();
  //e45_class.display_draw();
  //e45_class.draw_evtgen();
  //e45_class.draw_multiplicity();
  //e45_class.draw_edep();
  //e45_class.draw_elastic();

}//end
