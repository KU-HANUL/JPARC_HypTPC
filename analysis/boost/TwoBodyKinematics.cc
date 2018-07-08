#include <TMath.h>

const double m_gamma = 0.;
const double m_proton = 0.938272;
const double m_neutron = 0.939565;
const double m_deuteron = 0.9315016 * 2.01;
const double m_pi0 = 0.134977;

double r = 1;
double length1 = 3;
double length2 = 2;

void Decay2(TLorentzVector lv_mother,
	    TLorentzVector &lv1, TLorentzVector &lv2,
	    double m1, double m2);
double triangle(double a, double b, double c);
double Getcos(TLorentzVector lv1, TLorentzVector lv2);
double GetAng(TLorentzVector lv1, TLorentzVector lv2);

void TwoBodyKinematics() {
  
  double p_neutron = 0.8;
  double e_neutron = sqrt(p_neutron*p_neutron + m_neutron*m_neutron);

  TLorentzVector lv_beam(0.0, 0.0, p_neutron, e_neutron);
  TLorentzVector lv_target(0.0, 0.0 , 0.0 , m_proton);
  TLorentzVector lv_total;
  lv_total  = lv_beam + lv_target;

  TLorentzVector lv_deuteron;
  TLorentzVector lv_pi0;

  const int nevent = 100000;

  TH1D* h1 = new TH1D("h1","np #rightarrow d#pi0 ; cos#theta ",220, 0.0, 1.1);
  TH1D* h_1 = new TH1D("h_1","np #rightarrow d#pi0 ; #theta",180,0,20);

  TH1D* h2 = new TH1D("h2","#pi0 #rightarrow 2#gamma ; cos#theta",220, -1.1, 1.1);
  TH1D* h_2 = new TH1D("h_2","#pi0 #rightarrow 2#gamma ; #theta",190,-5,185);
  
  TH1D* h3 = new TH1D("h3","#pi0 #rightarrow 2#gamma ; cos#theta", 220, -1.1, 1.1);
  TH1D* h_3 = new TH1D("h_3","#pi0 #rightarrow 2#gamma ; #theta",180,0,180); 

  TH1D* h4 = new TH1D("h4","np #rightarrow d#pi0 #rightarrow d2#gamma ; cos#theta",220,-1.1,1.1);
  TH1D* h4_1 = new TH1D("h4_1","np #rightarrow d#pi0 #rightarrow d2#gamma ; cos#theta", 220, -1.1, 1.1);
  
  TH1D* h_4 = new TH1D("h_4","np #rightarrow d#pi0 #rightarrow d2#gamma ; theta", 180, 0, 180);
  TH1D* h_4_1 = new TH1D("h_4_1","np #rightarrow d#pi0 #rightarrow d2#gamma ; theta", 180, 0, 180);
  
  TH2D* hist1 = new TH2D("hist1","z=+1.5 gamma1",60,-1.5,1.5,60,-1.5,1.5);
  TH2D* hist2 = new TH2D("hist2","z=+1.5 gamma2",60,-1.5,1.5,60,-1.5,1.5);
  
  TH2D* histo1 = new TH2D("histo1","z=+1 gamma1",60,-1.5,1.5,60,-1.5,1.5);
  TH2D* histo2 = new TH2D("histo2","z=+1 gamma2",60,-1.5,1.5,60,-1.5,1.5);

  TH1D* histo1_x = new TH1D("histo1_x","z=+1 gamma1 X", 60, -1.5, 1.5 );
  TH1D* histo1_y = new TH1D("histo1_y","z=+1 gamma1 Y", 60, -1.5, 1.5 );

  TH1D* histo2_x = new TH1D("histo2_x","z=+1 gamma2 X", 60, -1.5, 1.5 );
  TH1D* histo2_y = new TH1D("histo2_y","z=+1 gamma2 Y", 60, -1.5, 1.5 );
  
  for (int i = 0 ; i<nevent ; i++) {
    Decay2(lv_total, lv_deuteron, lv_pi0, m_deuteron, m_pi0);

    TLorentzVector lv_gamma1;
    TLorentzVector lv_gamma2;
   
    double theta1 = GetAng(lv_total,lv_deuteron);
    
    h1->Fill(Getcos(lv_deuteron,lv_total));
    h_1->Fill(theta1);
    //    h_1->Fill(lv_deuteron.Theta()*180/TMath::Pi());
    
    Decay2(lv_pi0, lv_gamma1, lv_gamma2, m_gamma, m_gamma);

    h2->Fill(Getcos(lv_pi0,lv_gamma1));
   
    double theta2 = GetAng(lv_pi0,lv_gamma1);

    //  h_2->Fill(theta2);
    h_2->Fill(lv_gamma1.Theta()*180/TMath::Pi());
    
    double theta_2 = GetAng(lv_pi0,lv_gamma2);
    
    h3->Fill(Getcos(lv_pi0,lv_gamma2));
    h_3->Fill(theta_2);

    h4->Fill(Getcos(lv_total, lv_gamma1));
    h4_1->Fill(Getcos(lv_total,lv_gamma2));

    double theta3 = GetAng(lv_total,lv_gamma1);
    double theta3_1 = GetAng(lv_total,lv_gamma2);
    
    h_4->Fill(theta3);
    h_4_1->Fill(theta3_1);
    
    double p_t1 = sqrt(lv_gamma1.Pz()*lv_gamma1.Pz());
    double p_t2 = sqrt(lv_gamma2.Pz()*lv_gamma2.Pz());

    double t1 = (length1/2)/p_t1;
    double t2 = (length1/2)/p_t2;

    double t_1 = (length2/2)/p_t1;
    double t_2 = (length2/2)/p_t2;
    
    TVector3 hit1(t1*lv_gamma1.Px(),t1*lv_gamma1.Py(),t1*lv_gamma1.Pz());
    TVector3 hit2(t2*lv_gamma2.Px(),t2*lv_gamma2.Py(),t2*lv_gamma2.Pz());

    TVector3 hit_1(t_1*lv_gamma1.Px(),t_1*lv_gamma1.Py(),t_1*lv_gamma1.Pz());
    TVector3 hit_2(t_2*lv_gamma2.Px(),t_2*lv_gamma2.Py(),t_2*lv_gamma2.Pz());
    
    //    cout << hit1.X() << "\t" << hit1.Y() << "\t" << hit1.Z() <<endl;
    //    cout << hit2.X() << "\t" << hit2.Y() << "\t" << hit2.Z() << endl;

    double len1 = sqrt(hit1.X()*hit1.X() + hit1.Y()*hit1.Y());
    double len2 = sqrt(hit2.X()*hit2.X() + hit2.Y()*hit2.Y());

    double len_1 = sqrt(hit_1.X()*hit_1.X() + hit_1.Y()*hit_1.Y());
    double len_2 = sqrt(hit_2.X()*hit_2.X() + hit_2.Y()*hit_2.Y());
    
    if ( len1 < 1 && hit1.Z() == 1.5) 
      hist1->Fill(hit1.X(),hit1.Y());
    
    if (len_1 < 1 && hit_1.Z() == 1) {
      histo1->Fill(hit_1.X(),hit_1.Y());
      histo1_x->Fill(hit_1.X());
      histo1_y->Fill(hit_1.Y());
    }
    if (len2 < 1 && hit2.Z() == 1.5)
      hist2->Fill(hit2.X(),hit2.Y());

    if (len_2 < 1 && hit_2.Z() == 1) {
      histo2->Fill(hit_2.X(),hit_2.Y());
      histo2_x->Fill(hit_2.X());
      histo2_y->Fill(hit_2.Y());
    }
  }
  
  auto *c1 = new TCanvas("c1","c1",600,1200);
  auto *c2 = new TCanvas("c2","c2",600,1200);
  auto *c12 = new TCanvas("c12","c12",600,1200);
  auto *c22 = new TCanvas("c22","c22",1200,1200);
  auto *c3 = new TCanvas("c3","c3",600,1200);
  auto *c4_1 = new TCanvas("c3_1","c3_1",1200,1200);
  auto *c4 = new TCanvas("c4","c4",600,1200);
  // auto *c4_2 = new TCanvas("c4_1","c4_1",600,1200);
  
  c1->Divide(1,2);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h_1->Draw();

  c2->Divide(1,2);
  c2->cd(1);
  h2->Draw();
  c2->cd(2);
  h_2->Draw();

  c12->Divide(1,2);
  c12->cd(1);
  h3->Draw();
  c12->cd(2);
  h_3->Draw();

  c22->Divide(2,2);
  c22->cd(1);
  h4->Draw();
  c22->cd(2);
  h4_1->Draw();
  c22->cd(3);
  h_4->Draw();
  c22->cd(4);
  h_4_1->Draw();
  
  c3->Divide(1,2);
  c3->cd(1);
  hist1->Draw("colz");
  c3->cd(2);
  hist2->Draw("colz");

  c4_1->Divide(2,2);
  c4_1->cd(1);
  histo1_x->Draw();
  c4_1->cd(2);
  histo1_y->Draw();
  c4_1->cd(3);
  histo2_x->Draw();
  c4_1->cd(4);
  histo2_y->Draw();
  
  c4->Divide(1,2);
  c4->cd(1);
  histo1->Draw("colz");
  c4->cd(2);
  histo2->Draw("colz");
}

void Decay2(TLorentzVector lv_mother, TLorentzVector &lv1, TLorentzVector &lv2,
	   double m1, double m2) {
  TRandom3 ran_gen(0);
  double p;
  double m_mother;
  m_mother = lv_mother.Mag();
  p = (m_mother/2.0) * sqrt(triangle(1,m1/m_mother,m2/m_mother));

  double xx,yy,zz;
  ran_gen.Sphere(xx,yy,zz,1);
  lv1.SetXYZT(p*xx, p*yy, p*zz, sqrt(p*p+m1*m1));
  lv2.SetXYZT(-p*xx, -p*yy, -p*zz, sqrt(p*p+m2*m2));

  TVector3 beta_lab(lv_mother.X()/lv_mother.E(),
		    lv_mother.Y()/lv_mother.E(),
		    lv_mother.Z()/lv_mother.E());
  lv1.Boost(beta_lab);
  lv2.Boost(beta_lab);
}

double triangle(double a, double b, double c) {
  double value = ((a*a)-((b+c)*(b+c)))*((a*a)-((b-c)*(b-c)));
  return value;
}

double Getcos(TLorentzVector lv1, TLorentzVector lv2) {
  double dot_product = lv1.Px()*lv2.Px() + lv1.Py()*lv2.Py() + lv1.Pz()*lv2.Pz();
  double cos = dot_product / (lv1.P() * lv2.P());
  return cos;
}

    
double GetAng(TLorentzVector lv1, TLorentzVector lv2) {
  //  double rad = 57.2958;
  double ang = (180*TMath::ACos(Getcos(lv1,lv2)))/TMath::Pi();
  return ang;
}
