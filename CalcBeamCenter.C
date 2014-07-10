#include "UtilFns.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
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
TGraph *DcaDist(geoEvents &events, TVectorD &bc, TString arm,
                TH1D *hr, int ntracks=50000);
TH2D *DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
               const char *name, const char *title);
TH2D *DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
               const char *name, const char *title);

void CalcBeamCenter()
{
  TFile *f = new TFile("rootfiles/411768_july3_parv1_small.root");
  TNtuple *t = (TNtuple *)f->Get("seg_clusntuple");
  gStyle->SetOptStat(0);
  int nbc = 10000;

  SvxTGeo *geo = new SvxTGeo;
  geo->ReadParFile("geom/svxPISA-411768.par");
  geo->MakeTopVolume(100, 100, 100);
  geo->AddSensors();
  geo->GeoManager()->CloseGeometry();

#if 0
  TH1D *hne = new TH1D("hne", "East;tracks/event", 100, 0, 100);
  TH1D *hnw = new TH1D("hnw", "West;tracks/event", 100, 0, 100);
  geoEvents events;
  GetEventsFromTree(t, geo, events);
  FitTracks(events);

  // for (unsigned int ev=0; ev<events.size(); ev++)
  // {
  //   int ntrk = events[ev].size();
  //   int ne=0, nw=0;
  //   for (int t=0; t<ntrk; t++)
  //   {
  //     FitTrack(events[ev][t]);
  //     SvxGeoTrack trk = events[ev][t];
  //     if (East(trk.phi0)) ne++;
  //     else nw++;

  //     // Printf("%d %f %f %f %f",
  //     //        trk.nhits,
  //     //        trk.phi0,
  //     //        trk.the0,
  //     //        trk.vy,
  //     //        trk.vz);
  //   }
  //   hne->Fill(ne);
  //   hnw->Fill(nw);
  // }
  // DrawObject(hne); gPad->SetLogx(); gPad->SetLogy();
  // DrawObject(hnw); gPad->SetLogx(); gPad->SetLogy();

  Printf("Computing east arm beam center...");
  TVectorD bce = BeamCenter(events, nbc, "east");
  TH1D *he = new TH1D("he", "he", 100, 0, 0.5);
  TGraph *ge = DcaDist(events,bce,"east",he);

  Printf("Computing west arm beam center...");
  TVectorD bcw = BeamCenter(events, nbc, "west");
  TH1D *hw = new TH1D("hw", "hw", 100, 0, 0.5);
  TGraph *gw = DcaDist(events,bcw,"west",hw);

  TH2D *hbp = DcaVsPhi(events, bce, bcw, "hbp", "DCA to beam center vs  #phi");
  DrawObject(hbp, "colz", "dcavsphi");
  gPad->Print("pdfs/dcavsphi.pdf");
#endif
  
#if 1
  geoTracks tracks;
  GetTracksFromTree(t, geo, tracks);
  FitTracks(tracks);

  Printf("Computing east arm beam center...");
  TVectorD bce = BeamCenter(tracks, nbc, "east");
  TH1D *he = new TH1D("he", "he", 100, 0, 0.5);
  TGraph *ge = DcaDist(tracks,bce,"east",he);

  Printf("Computing west arm beam center...");
  TVectorD bcw = BeamCenter(tracks, nbc, "west");
  TH1D *hw = new TH1D("hw", "hw", 100, 0, 0.5);
  TGraph *gw = DcaDist(tracks,bcw,"west",hw);

  TH2D *hbp = DcaVsPhi(tracks, bce, bcw, "hbp", "DCA to beam center vs  #phi");
  DrawObject(hbp, "colz", "dcavsphi");
  gPad->Print("pdfs/dcavsphi.pdf");
#endif

  // ---- ---- ---- ---- ---- ---- ----
  //                Plot
  // ---- ---- ---- ---- ---- ---- ----
  TCanvas *c = new TCanvas("c", "dca2d", 500, 500);
  TH1F *h = c->DrawFrame(-0.3, -0.3, 0.3, 0.3,
                         "Beam center DCA #color[861]{ east}, #color[800]{ west}");
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
  hw->SetTitle("Beam center DCA #color[861]{ east}, #color[800]{ west};DCA [cm];");
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
    if ((arm=="east" && East(*phi)) || (arm=="west" && !East(*phi)))
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
    if ((arm=="east" && East(phi)) || (arm=="west" && !East(phi)))
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

TGraph *
DcaDist(geoEvents &events, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  int i = 0;
  TGraph *g = new TGraph(ntracks);
  g->SetMarkerStyle(kFullDotMedium);

  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      double phi = trk.phi0;
      if ((arm=="east" && East(phi)) || (arm=="west" && !East(phi)))
      {
        TVectorD a(2); a(1) = trk.vy;
        TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
        TVectorD d = IPVec(a,n,bc);  // d = a - bc - ((a - bc)*n)*n;

        if (i<ntracks)
        {
          g->SetPoint(i, d(0), d(1));
          i++;
        }

        if (hr)
          hr->Fill(TMath::Sqrt(d*d));
      }
    }

  return g;
}

TH2D *
DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
         const char *name, const char *title)
{
  TH2D *h = new TH2D(name, title, 100, 0, TMath::TwoPi(), 100, -0.0, 0.1);

  for (unsigned int i=0; i<tracks.size(); i++)
  {
    double phi = tracks[i].phi0;
    TVectorD a(2); a(1) = tracks[i].vy;
    TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
    TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
    h->Fill(fmod(TMath::PiOver2()+phi,TMath::TwoPi()), TMath::Sqrt(d*d));
  }

  return h;
}

TH2D *
DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
         const char *name, const char *title)
{
  TH2D *h = new TH2D(name, title, 100, 0, TMath::TwoPi(), 100, -0.0, 0.1);

  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      double phi = trk.phi0;
      TVectorD a(2); a(1) = trk.vy;
      TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
      TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
      h->Fill(fmod(TMath::PiOver2()+phi,TMath::TwoPi()), TMath::Sqrt(d*d));
    }

  return h;
}
