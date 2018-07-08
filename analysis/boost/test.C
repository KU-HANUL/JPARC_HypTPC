void test() 
{
  // example of use of TGenPhaseSpace
  //Author: Valerio Filippini

  //   if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");
  TVector3* Dvec  = new TVector3();
  TVector3* Pivec = new TVector3();
  TFile *file = new TFile("Results.root","recreate");
  TTree *tree = new TTree("tree","GeV, rad");
  tree->Branch("D",Dvec);
  tree->Branch("Pi",Pivec);

  TLorentzVector target(0.0, 0.0, 0.0, 0.938); //proton
  TLorentzVector beam(0.0, 0.0, 0.8, sqrt(pow(.939,2)+pow(0.8,2))); // 0.8 GeV/c n beam
  TLorentzVector W = beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[2] = {1.8, 0.135}; // d, pi0 

  TGenPhaseSpace event;
  event.SetDecay(W, 2, masses);

  for (Int_t n=0;n<100000;n++)
  {
    Double_t weight = event.Generate();

    TLorentzVector *pD = event.GetDecay(0);
    TLorentzVector *pPi = event.GetDecay(1);

    *Dvec = pD->Vect();
    *Pivec = pPi->Vect();

    tree->Fill();

  }

  file->cd();
  tree->Write();
}
