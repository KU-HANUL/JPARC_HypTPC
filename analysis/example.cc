void example(){

  TH1D *hist =new TH1D("hist","",10,5,15);
  TRandom *eventgen = new TRandom();
  double r;
  for(int i=0;i<100000;i++){
    r=eventgen->Gaus(10,2);
    hist -> Fill(r);
  }
  hist -> Draw();
}
