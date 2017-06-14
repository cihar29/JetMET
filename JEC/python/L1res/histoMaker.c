//Chad Harrington 6/13/2017
//Execute as root -l -b -q 'histoMaker.c (mcName, dataName, ISMC, rCone)'

using namespace std;

const int nEta = 82;
float etabins[nEta+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

void FillHist1D(const TString& histName, const Double_t& value, double weight);
void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, double weight);
void getGeometry(double (&geo)[nEta][nEta], float rCone);
double areaS(double R, double x1, double x2);
double dist(double R, double x1, double x2);

map<TString, TH1*> m_Histos1D;
map<TString, TH2*> m_Histos2D;
map<TString, TProfile*> m_Profiles;
map<TString, TProfile2D*> m_Profiles2D;

const int MAXNPU = 60;
const int MAXNPV = 60;
const int MAXRHO = 60;

void histoMaker(TString mcName="/eos/cms/store/group/phys_jetmet/schoef/L1res/SingleNeutrino/tree.root",
                TString dataName="/eos/cms/store/group/phys_jetmet/schoef/L1res/ZeroBias_Run2016G-03Feb2017-v1/tree.root",
                bool ISMC=true, float rCone=0.4) {

  //Open Files//

  TString inName = ISMC ? mcName : dataName;
  TString outName = ISMC ? "histoMC" : "histoData";
  outName += "_R" + to_string( int(rCone*10) ) + ".root";
  TFile* inFile = TFile::Open(inName);

  TTree* tree = (TTree*) inFile->Get("tree");
  Long64_t nEntries = tree->GetEntries();
  cout << nEntries << " Events" << endl;

  //Declare Histos//

  enum Flavor{ all = 0, chm, chu, nh, ne, hfh, hfe, lep, untrk, numFlavors};
  TString ids[] = {"all", "chm", "chu", "nh", "ne", "hfh", "hfe", "lep", "untrk"};
  TString hname;

  for (int i_id=0; i_id<numFlavors; i_id++){
    for (int i_nPU=0; i_nPU<MAXNPU; i_nPU++){
      hname = Form("p_offset_eta_nPU%i_", i_nPU) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, nEta, etabins);
    }
    for (int i_nPV=0; i_nPV<MAXNPV; i_nPV++){
      hname = Form("p_offset_eta_nPV%i_", i_nPV) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, nEta, etabins);
    }
  }

  hname = "nPV";
  m_Histos1D[hname] = new TH1F(hname,hname,MAXNPV,0,MAXNPV);
  hname = "nPU";
  m_Histos1D[hname] = new TH1F(hname,hname,2*MAXNPU,0,MAXNPU);
  hname = "rho";
  m_Histos1D[hname] = new TH1F(hname,hname,2*MAXRHO,0,MAXRHO);

  hname = "p_nPV_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,2*MAXNPU,0,MAXNPU);
  hname = "p_rho_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,2*MAXNPU,0,MAXNPU);
  hname = "p_rho_nPV";
  m_Profiles[hname] = new TProfile(hname,hname,MAXNPV,0,MAXNPV);

  //Get Areas//

  double geo[nEta][nEta];
  getGeometry(geo, rCone);

  //Weighting//

  TH1F* h_weights = 0;
  if (ISMC) {

    TFile* dataFile = TFile::Open(dataName);
    TTree* dTree = (TTree*) dataFile->Get("tree");

    h_weights = new TH1F("h_weights","h_weights",2*MAXNPU,0,MAXNPU);
    dTree->Draw("mu>>h_weights");

    TH1F* h_muMC = new TH1F("h_muMC","h_muMC",2*MAXNPU,0,MAXNPU);
    tree->Draw("mu>>h_muMC");

    h_weights->Divide(h_muMC);
    h_weights->Scale( 1/ h_weights->GetMaximum() );
  }

  //Set Branches//

  float offenergy[numFlavors][nEta];
  float energy[nEta];
  float mu, rho;
  int nPV;

  tree->SetBranchAddress("Offset_e_all",  offenergy[all]);
  tree->SetBranchAddress("Offset_e_ch",   offenergy[chm]);
  tree->SetBranchAddress("Offset_e_unm",  offenergy[chu]);
  tree->SetBranchAddress("Offset_e_nh",   offenergy[nh]);
  tree->SetBranchAddress("Offset_e_ph",   offenergy[ne]);
  tree->SetBranchAddress("Offset_e_hhf",  offenergy[hfh]);
  tree->SetBranchAddress("Offset_e_ehf",  offenergy[hfe]);
  tree->SetBranchAddress("Offset_e_ele",  offenergy[lep]);
  tree->SetBranchAddress("Offset_e_lost", offenergy[untrk]);
  tree->SetBranchAddress("mu", &mu);
  tree->SetBranchAddress("rho", &rho);
  tree->SetBranchAddress("nVert", &nPV);  

  //Loop Over Entries//

  for (Long64_t n=0; n<nEntries; n++) {
    tree->GetEntry(n);

    float weight = ISMC ? h_weights->GetBinContent( h_weights->FindBin(mu) ) : 1.;

    FillHist1D("nPU", mu, weight);
    FillHist1D("nPV", nPV, weight);
    FillHist1D("rho", rho, weight);
    FillProfile("p_nPV_nPU", mu, nPV, weight);
    FillProfile("p_rho_nPU", mu, rho, weight);
    FillProfile("p_rho_nPV", nPV, rho, weight);

    int intmu = mu + 0.5;

    for (int ieta=0; ieta<nEta; ieta++){
      double eta = 0.5*(etabins[ieta] + etabins[ieta+1]);

      for (int i_id=0; i_id<numFlavors; i_id++){
        double offpt = 0;

        for(int jeta = ieta-10; jeta <= ieta+10; jeta++){
          if( jeta<0 || jeta+1 >= nEta) continue;

          offpt += offenergy[i_id][jeta] * geo[ieta][jeta];
        }
        hname = Form("p_offset_eta_nPU%i_", intmu) + ids[i_id];
        FillProfile(hname, eta, offpt, weight);

        hname = Form("p_offset_eta_nPV%i_", nPV) + ids[i_id];
        FillProfile(hname, eta, offpt, weight);
      }
    }
  } //end event loop

  //Write Histos//

  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

  for (int i_id=0; i_id<numFlavors; i_id++){
    outFile->mkdir( "offset_nPU/" + ids[i_id] );
    outFile->mkdir( "offset_nPV/" + ids[i_id] );
  }

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++)
    hid->second->Write();

  for (map<TString, TProfile*>::iterator hid = m_Profiles.begin(); hid != m_Profiles.end(); hid++){
    outFile->cd();
    hname = hid->first;

    if ( hname.Contains("p_offset_eta_nPU") ){
      TString id = hname( hname.Last('_')+1, hname.Length() );
      outFile->cd(outName + ":/offset_nPU/" + id);
    }
    else if ( hname.Contains("p_offset_eta_nPV") ){
      TString id = hname( hname.Last('_')+1, hname.Length() );
      outFile->cd(outName + ":/offset_nPV/" + id);
    }
    hid->second->Write();
  }

  outFile->Write();
  delete outFile;
  outFile = 0;
}

//Functions//

void FillHist1D(const TString& histName, const Double_t& value, double weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, double weight) 
{
  map<TString, TProfile*>::iterator hid=m_Profiles.find(histName);
  if (hid==m_Profiles.end())
    cout << "%FillProfile -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

//Sum the eta strips within a cone of rCone centered around each etabin
void getGeometry(double (&geo)[nEta][nEta], float rCone){

  for (int ieta=0; ieta<nEta; ieta++){
    double eta = 0.5*(etabins[ieta] + etabins[ieta+1]);

    for(int jeta = ieta-10; jeta <= ieta+10; jeta++){   

      if( jeta<0 || jeta+1 >= nEta) continue;
    
      double etaL = etabins[jeta];                    // left  edge of the eta strip
      double etaR = etabins[jeta+1];                  // right edge of the eta strip
      double etaC = 0.5*(etaL+etaR);                  // center of the eta strip
      double A = areaS(rCone, etaL-eta, etaR-eta);    // area of the eta strip inside the cone

      // the next lines make sure that we do vector addition in phi direction;
      // We would be doing scalar addition with coef = 1.

      double dphi = dist(rCone, etaL-eta, etaR-eta);
      double coef = dphi > 0.000001 ? TMath::Sin(dphi) / dphi : 1.;

      //geometric factor and factor to convert energy to pT
      geo[ieta][jeta] = A * coef / (2*TMath::Pi()*(etaR-etaL)) / cosh(etaC);
    }
  }
}

double areaS(double R, double x1, double x2){
//
// Area inside a shape delineated by a circle of radius R
// centered at (0,0) and two parallel lines at x=x1 and x=x2
//

   if( R<=0. || x1==x2 ) return 0.;
   if(x1*x2>0 && fabs(x1) > R && fabs(x2) > R) return 0. ;

   double d1 = fabs(x1) ;
   double d2 = fabs(x2) ;
   if (d1>d2){
     double d = d1;
     d1 = d2;
     d2 = d ;
   }

   // area of segment at distance d1 from the center
   double theta1 = 2.* TMath::ACos(d1/R) ;
   double A1 = (0.5*R*R*(theta1-TMath::Sin(theta1))) ;

   // area of segment at distance d2 from the center
   double A2 = 0. ;
   if(d2<=R){
     double theta2 = 2.* TMath::ACos(d2/R) ;
     A2 = (0.5*R*R*(theta2-TMath::Sin(theta2))) ;
   }

   if(x1*x2>=0){ // both lines on the same side from the center
     return A1 - A2 ;
   } else { // the lines on the opposite side from the center
    return TMath::Pi()*R*R - A1 - A2 ;
   }
}

double dist(double R, double x1, double x2){
//
// Take a circle of radius R centered at (0,0).
// This function calculates distance between a 
// midpoint of x=x1 and x=x2 and the point on circle rim 
// along vertical line.
//

   if( R<=0.) return 0.;
   if(x1*x2>0 && fabs(x1) >= R && fabs(x2) >= R) return 0. ; // both lines outside the circle 
                                                             // and on the same side of origin. 
  if(fabs(x1)>R) x1 = TMath::Sign(R, x1);
  if(fabs(x2)>R) x2 = TMath::Sign(R, x2); 

   double x = 0.5*( x1 + x2) ;

   return TMath::Sqrt(R*R-x*x);
}
