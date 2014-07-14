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
#include <TProfile.h>

// Solve the linear system y0[i] = -m[i]*bc0 + bc1 for i..ntrk-1 tracks.
// m[i] is the slope (tan(phi)), y0[i] is the y-intercept, and (bc0,bc1)
// is the least-squares beam (x,y) position.

TGraph *DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=10000);
TGraph *DcaDist(geoTracks &tracks, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=10000);
TGraph *DcaDist(geoEvents &events, TVectorD &bc, TString arm,
                TH1D *hr, int ntracks=10000);
TH2D *DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
               const char *name, const char *title);
TH2D *DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
               const char *name, const char *title);
void DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
              TH2D *hist, TProfile *prof);

void CalcBeamCenter(int run = 411768, int iter = 2)
{
  int nbc = 10000;
  TString pisaFileIn = Form("geom/svxPISA-%d.par", run);
  if (iter > 0)
    pisaFileIn += Form(".%d", iter);

  TString fnames[] = {"",
                      "july3_parv1_small",
                      "july11_v3_500kevents"
                     };
  TFile *f = new TFile(Form("rootfiles/%d_%s.root", run, fnames[iter].Data()));
  TNtuple *t = (TNtuple *)f->Get("seg_clusntuple");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(56, 0, 0.5); // Inverted "radiator", 50% transparency
  TObjArray *cList = new TObjArray();

  SvxTGeo *geo = new SvxTGeo;
  geo->ReadParFile(pisaFileIn.Data());
  geo->MakeTopVolume(100, 100, 100);
  geo->AddSensors();
  geo->GeoManager()->CloseGeometry();

#if 1
  geoEvents events;
  GetEventsFromTree(t, geo, events);
  FitTracks(events);

  // Fill event multiplicity histograms
  TH1D *hne = new TH1D("hne", "East;tracks/event", 100, 0, 100);
  TH1D *hnw = new TH1D("hnw", "West;tracks/event", 100, 0, 100);
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ne=0, nw=0;
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (East(trk.phi0)) ne++;
      else nw++;
      if (false)
        Printf("%d %f %f %f %f", trk.nhits, trk.phi0, trk.the0, trk.vy, trk.vz);
    }
    hne->Fill(ne);
    hnw->Fill(nw);
  }

  Printf("Computing east arm beam center...");
  TVectorD bce = BeamCenter(events, nbc, "east", "print");
  TH1D *he = new TH1D("he", "he", 100, 0, 0.5);
  TGraph *ge = DcaDist(events,bce,"east",he);

  Printf("Computing west arm beam center...");
  TVectorD bcw = BeamCenter(events, nbc, "west", "print");
  TH1D *hw = new TH1D("hw", "hw", 100, 0, 0.5);
  TGraph *gw = DcaDist(events,bcw,"west",hw);

  // "DCA to beam center vs  #phi",
  double dmin = 0.0, dmax = 0.1;
  TH2D *hbp = new TH2D("hbp", "", 100, 0, TMath::TwoPi(), 100, dmin, dmax);
  TProfile *prof = new TProfile("prof", "", 100, 0, TMath::TwoPi(), dmin, dmax);
  prof->SetMarkerStyle(kFullCircle);

  // Error calculation options:
  // "" : error = sigma/sqrt(N) = std. dev of mean.
  // "s": error = sigma = std. dev of distribution in dmin, dmax.
  // "g": error = 1/sqrt(sum(w)) where w is (weighted) sum of entries.
  prof->BuildOptions(dmin, dmax, "g");
  DcaVsPhi(events, bce, bcw, hbp, prof);
#endif

#if 0
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
  DrawObject(hne, "", "east_mult", cList);
  gPad->SetLogx();
  gPad->SetLogy();
  DrawObject(hnw, "", "west_mult", cList);
  gPad->SetLogx();
  gPad->SetLogy();

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
  cList->Add(c);

  TCanvas *cr = new TCanvas("cr", "dcar", 500, 500);
  SetHistProps(he,kAzure+1,kNone,kAzure+1);
  SetHistProps(hw,kOrange,kNone,kOrange);
  he->SetLineWidth(2);
  hw->SetLineWidth(2);
  hw->SetTitle("Beam center DCA #color[861]{ east}, #color[800]{ west};DCA [cm];");
  hw->Draw("");
  he->Draw("same");
  cr->SetLogy();
  cList->Add(cr);

  DrawObject(hbp, "colz", "dcavsphi", cList, 900, 500);
  hbp->SetTitle("DCA to beam center");
  TAxis *ax = hbp->GetXaxis();
  TAxis *ay = hbp->GetYaxis();
  ax->SetTitle("#Leftarrow  west      #phi + #pi/2      east  #Rightarrow");
  ay->SetTitle("DCA [cm]");
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(408);
  ay->SetNdivisions(205);
  gPad->SetGridy();
  prof->Draw("same");

  Printf("Fraction outside 1000um: %.3f (e) %.3f (w)",
         1.0 - he->Integral(1,he->FindBin(0.099))/he->Integral(),
         1.0 - hw->Integral(1,hw->FindBin(0.099))/hw->Integral());

  PrintPDFs(cList, Form("pdfs/iter%d", iter), "");
  PrintPDF(cList, Form("pdfs/beam-center-iter%d", iter));

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

void
DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
         TH2D *hist, TProfile *prof)
{
  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      double phi = trk.phi0;
      TVectorD a(2); a(1) = trk.vy;
      TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
      TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
      double dmag = TMath::Sqrt(d*d);
      double phir = fmod(TMath::PiOver2()+phi,TMath::TwoPi());
      hist->Fill(phir, dmag);
      prof->Fill(phir, dmag);
    }

  return;
}
