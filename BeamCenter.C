#include "GLSFitter.h"
#include "UtilFns.h"

// Solve the linear system y0[i] = -m[i]*bc0 + bc1 for i..ntrk-1 tracks.
// m[i] is the slope (tan(phi)), y0[i] is the y-intercept, and (bc0,bc1)
// is the least-squares beam (x,y) position.

TGraph *DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm, 
                TH1D* hr=0, int ntracks=50000);

void BeamCenter()
{
  TFile *f = new TFile("rootfiles/411768_cluster.1.root");
  TNtuple *t = (TNtuple *)f->Get("trktree");
  gStyle->SetOptStat(0);
  int ntrk = 5000;
  TMatrixD M(ntrk,2);
  TVectorD y0(ntrk);
  TMatrixD L(ntrk, ntrk);
  L.UnitMatrix();
  L *= 0.01;
  TMatrixD cov(2,2);

  Printf("Computing east arm beam center...");
  GetData(f,t,y0,M,"east");
  TVectorD bce = SolveGLS(M, y0, L, cov);
  Printf("East x,y (%.3f +- %.3f, %.3f +- %.3f)",
         bce(0), TMath::Sqrt(cov(0,0)),
         bce(1), TMath::Sqrt(cov(1,1)));
  cov.Print();
  TH1D* he = new TH1D("he", "he", 100,0,0.5);
  TGraph* ge = DcaDist(f,t,bce,"east",he);

  Printf("Computing west arm beam center...");
  GetData(f,t,y0,M,"west");
  TVectorD bcw = SolveGLS(M, y0, L, cov);
  Printf("West x,y (%.3f +- %.3f, %.3f +- %.3f)",
         bcw(0), TMath::Sqrt(cov(0,0)),
         bcw(1), TMath::Sqrt(cov(1,1)));
  cov.Print();
  TH1D* hw = new TH1D("hw", "hw", 100,0,0.5);
  TGraph* gw = DcaDist(f,t,bcw,"west",hw);


  TCanvas* c = new TCanvas("c", "dca2d", 500, 500);
  TH1F* h = c->DrawFrame(-0.3, -0.3, 0.3, 0.3, 
                         "DCA #color[861]{east}, #color[800]{west}");
  h->GetXaxis()->SetTitle("x [cm]");
  h->GetYaxis()->SetTitle("y [cm]");
  TEllipse ell(0,0,0.2,0.2);
  ell.Draw("same");
  SetGraphProps(ge,kNone,kNone,kAzure+1,kFullDotMedium);
  SetGraphProps(gw,kNone,kNone,kOrange,kFullDotMedium);
  ge->Draw("psame");
  gw->Draw("psame");

  c->Print("pdfs/dca.pdf");
  TCanvas* cr = new TCanvas("cr", "dcar", 500, 500);
  SetHistProps(he,kAzure+1,kNone,kAzure+1);
  SetHistProps(hw,kOrange,kNone,kOrange);
  he->SetLineWidth(2);
  hw->SetLineWidth(2);
  hw->SetTitle("DCA #color[861]{east}, #color[800]{west};DCA [cm];");
  hw->Draw("");
  he->Draw("same");
  cr->SetLogy();
  cr->Print("pdfs/dcar.pdf");

  Printf("Fraction outside 2000um: %.3f (e) %.3f (w)", 
         1.0 - he->Integral(1,he->FindBin(0.199))/he->Integral(),
         1.0 - hw->Integral(1,hw->FindBin(0.199))/hw->Integral());
  return;
}

void
GetData(TFile *f, TNtuple *t, TVectorD &y0, TMatrixD &M, TString arm)
{
  assert(y0.GetNrows()==M.GetNrows());
  assert(M.GetNcols()==2);

  TTreeReader r(t->GetName(), f);
  TTreeReaderValue<float> ty0(r, "y0");
  TTreeReaderValue<float> phi(r, "phi");

  int i=0;
  while (i<y0.GetNrows() && r.Next())
  {
    double m = TMath::Tan(*phi);
    bool east = (*phi > 0.5*TMath::Pi() && *phi < 1.5*TMath::Pi());
    if ((arm=="east" && east) || (arm=="west" && !east))
    {
      M(i, 0) = -m;
      M(i, 1) = 1;
      y0(i) = *ty0;
      i++;
    }
  }

  return;
}

TGraph *
DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm, TH1D* hr, int ntracks)
{
  TTreeReader r(t->GetName(), f);
  TTreeReaderValue<float> ty0(r, "y0");
  TTreeReaderValue<float> phi(r, "phi");
  TGraph *g = new TGraph();
  g->SetMarkerStyle(kFullDotMedium);
  // g->SetMarkerSize(0.5);
  int i=0;
  while (r.Next())
  {
    double m = TMath::Tan(*phi);
    bool east = (*phi > 0.5*TMath::Pi() && *phi < 1.5*TMath::Pi());
    if ((arm=="east" && east) || (arm=="west" && !east))
    {
      TVectorD a(2); a(1) = *ty0;
      TVectorD n(2); n(0) = TMath::Cos(*phi); n(1) = TMath::Sin(*phi);
      TVectorD d = a - bc - ((a - bc)*n)*n;
      
      if (i<ntracks)
        g->SetPoint(i, d(0), d(1));
      if (hr)
        hr->Fill(TMath::Sqrt(d*d));

      i++;
    }
  }
  return g;
}

// TVectorD
// DCA(TVectorD &a, TVectorD &n, TVectorD &p)
// {
// Compute distance from point p to line x = a + tn where
// a is a point on vector x and n is a unit vector directing x.
//   return a - p - ((a - p)*n)*n;
// }
