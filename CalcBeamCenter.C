#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "UtilFns.h"
#include "BeamCenterFuncs.h"
#include "ParameterDefs.h"
#include "GLSFitter.h"
#include "VtxIO.h"

#include <TEllipse.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGeoManager.h>

// Solve the linear system y0[i] = -m[i]*bc0 + bc1 for i..ntrk-1 tracks.
// m[i] is the slope (tan(phi)), y0[i] is the y-intercept, and (bc0,bc1)
// is the least-squares beam (x,y) position.

TGraph *DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=50000);
TGraph *DcaDist(geoTracks &tracks, TVectorD &bc, TString arm, 
                TH1D *hr=0, int ntracks=50000);

void CalcBeamCenter()
{
  // Initial declarations & assignments
  TFile *f = new TFile("rootfiles/411768_cluster.root");
  TNtuple *t = (TNtuple *)f->Get("seg_clusntuple");
  gStyle->SetOptStat(0);
  int ntrk = 5000;

  SvxTGeo *geo = new SvxTGeo;
  geo->ReadParFile("geom/svxPISA-411768.par");
  geo->MakeTopVolume(100, 100, 100);
  geo->AddSensors();
  geo->GeoManager()->CloseGeometry();

  geoTracks tracks;
  GetTracksFromTree(t, geo, tracks);
  FitTracks(tracks);

  Printf("Computing east arm beam center...");
  TVectorD bce = BeamCenter(tracks, 15000, "east");
  TH1D *he = new TH1D("he", "he", 100,0,0.5);
  TGraph *ge = DcaDist(tracks,bce,"east",he);

  Printf("Computing west arm beam center...");
  TVectorD bcw = BeamCenter(tracks, 15000, "west");
  TH1D *hw = new TH1D("hw", "hw", 100,0,0.5);
  TGraph *gw = DcaDist(tracks,bcw,"west",hw);

  // ---- ---- ---- ---- ---- ---- ----
  //                Plot
  // ---- ---- ---- ---- ---- ---- ----
  TCanvas *c = new TCanvas("c", "dca2d", 500, 500);
  TH1F *h = c->DrawFrame(-0.3, -0.3, 0.3, 0.3,
                         "Beam center DCA #color[861]{east}, #color[800]{west}");
  h->GetXaxis()->SetTitle("x [cm]");
  h->GetYaxis()->SetTitle("y [cm]");
  TEllipse ell(0,0,0.1,0.1);
  ell.Draw("same");
  SetGraphProps(ge,kNone,kNone,kAzure+1,kFullDotMedium);
  SetGraphProps(gw,kNone,kNone,kOrange,kFullDotMedium);
  ge->Draw("psame");
  gw->Draw("psame");

  c->Print("pdfs/dca.pdf");
  TCanvas *cr = new TCanvas("cr", "dcar", 500, 500);
  SetHistProps(he,kAzure+1,kNone,kAzure+1);
  SetHistProps(hw,kOrange,kNone,kOrange);
  he->SetLineWidth(2);
  hw->SetLineWidth(2);
  hw->SetTitle("Beam center DCA #color[861]{east}, #color[800]{west};DCA [cm];");
  hw->Draw("");
  he->Draw("same");
  cr->SetLogy();
  cr->Print("pdfs/dcar.pdf");

  Printf("Fraction outside 1000um: %.3f (e) %.3f (w)",
         1.0 - he->Integral(1,he->FindBin(0.099))/he->Integral(),
         1.0 - hw->Integral(1,hw->FindBin(0.099))/hw->Integral());
  return;
}

TGraph *
DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  // This function is no longer being used.
  TTreeReader r(t->GetName(), f);
  TTreeReaderValue<float> ty0(r, "y0");
  TTreeReaderValue<float> phi(r, "phi");
  TGraph *g = new TGraph();
  g->SetMarkerStyle(kFullDotMedium);

  int i=0;
  while (r.Next())
  {
    double m = TMath::Tan(*phi);
    bool east = (*phi > 0.5*TMath::Pi() && *phi < 1.5*TMath::Pi());
    if ((arm=="east" && east) || (arm=="west" && !east))
    {
      TVectorD a(2); a(1) = *ty0;
      TVectorD n(2); n(0) = TMath::Cos(*phi); n(1) = TMath::Sin(*phi);
      TVectorD d = IPVec(a,n,bc);  // d = a - bc - ((a - bc)*n)*n;

      if (i<ntracks)
        g->SetPoint(i, d(0), d(1));
      if (hr)
        hr->Fill(TMath::Sqrt(d*d));

      i++;
    }
  }
  return g;
}

TGraph *
DcaDist(geoTracks &tracks, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  TGraph *g = new TGraph();
  g->SetMarkerStyle(kFullDotMedium);
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    double phi = tracks[i].phi0;
    double m = TMath::Tan(phi);
    bool east = (phi > 0.5*TMath::Pi() && phi < 1.5*TMath::Pi());
    if ((arm=="east" && east) || (arm=="west" && !east))
    {
      TVectorD a(2); a(1) = tracks[i].vy;
      TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
      TVectorD d = IPVec(a,n,bc);  // d = a - bc - ((a - bc)*n)*n;

      if (i<(unsigned int)ntracks)
        g->SetPoint(i, d(0), d(1));
      if (hr)
        hr->Fill(TMath::Sqrt(d*d));
    }
  }
  return g;
}

