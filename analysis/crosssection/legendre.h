#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TChain.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <iostream>
#include <TPad.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooFFTConvPdf.h>
#include <RooPlot.h>
#include <RooFormula.h>

using namespace RooFit;
using namespace std;

static int legendre_class = 0;

class legendre{

 private:

  TCanvas* can[100] = {0};
  TH1D *hist_pdf[100] = {0};
  TGraphErrors *gr[100] = {0};

 public:

  void make_pdf(double *amplitude);

  legendre() = default;
};

void legendre::make_pdf(double *amplitude)
{
  //gStyle -> SetOptFit(1000);

  can[legendre_class] = new TCanvas(Form("can%d",legendre_class),"",1200,1200);

  double min_hist = -1;
  double max_hist = 1;

  hist_pdf[legendre_class] = new TH1D(Form("hist_pdf[%d]",legendre_class),Form("hist_pdf[%d]",legendre_class),100,-1,1);

  // --- Observable --
  RooRealVar COS("COS","Cos(#theta_{C.M.})",min_hist,max_hist);

  // --- Build Gaussian signal PDF ---
  RooLegendre L0("L0","Legendre l=0 PDF",COS,0);
  RooLegendre L1("L1","Legendre l=1 PDF",COS,1);
  RooLegendre L2("L2","Legendre l=2 PDF",COS,2);
  RooLegendre L3("L3","Legendre l=3 PDF",COS,3);
  RooLegendre L4("L4","Legendre l=4 PDF",COS,4);
  RooLegendre L5("L5","Legendre l=5 PDF",COS,5);
  RooLegendre L6("L6","Legendre l=6 PDF",COS,6);
  RooLegendre L7("L7","Legendre l=7 PDF",COS,7);
  RooLegendre L8("L8","Legendre l=8 PDF",COS,8);
  RooLegendre L9("L9","Legendre l=9 PDF",COS,9);

  // --- Construct total PDF ---
  RooRealVar n0("n0","l=0",1000.);
  /*
  RooRealVar n1("n1","l=1",amplidute[0]);
  RooRealVar n2("n2","l=2",amplidute[1]);
  RooRealVar n3("n3","l=3",amplidute[2]);
  RooRealVar n4("n4","l=4",amplidute[3]);
  RooRealVar n5("n5","l=5",amplidute[4]);
  RooRealVar n6("n6","l=6",amplidute[5]);
  RooRealVar n7("n7","l=8",amplidute[6]);
  RooRealVar n8("n8","l=7",amplidute[7]);
  RooRealVar n9("n9","l=9",amplidute[8]);
  */
  //RooAddPdf total("total","total pdf",RooArgList(L0,L1,L2,L3,L4,L5,L6,L7,L8,L9),RooArgList(n0,n1,n2,n3,n4,n5,n6,n7,n8,n9));
  RooAddPdf total("total","total pdf",RooArgList(L0),RooArgList(n0));
  TF1 *f1 = total.asTF(RooArgList(L0) , RooArgList(n0) );

  RooPlot* plot = COS.frame();
  plot -> GetXaxis() -> SetTitle(" COS ");
  total.plotOn(plot, LineColor(2));
  total.plotOn(plot, Components(L0), LineColor(3), LineStyle(kDashed));
  plot -> GetYaxis() -> SetTitle(" Counts ");

  can[legendre_class] -> cd();
  //plot -> Draw();
  f1 -> Draw();

  legendre_class++;
}
