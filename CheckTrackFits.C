#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"
#include "VtxIO.h"
#include "VtxVis.h"

void CheckTrackFits(int run = 411768,
                    int prod = 6,
                    int subiter = 20)
{
  gStyle->SetOptStat(0);
  TObjArray *cList = new TObjArray();
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);

  // Inputs:
  TString rootFileIn = Form("rootfiles/%d-%d-%d.root", run, prod, subiter);
  TString pisaFileIn = Form("geom/%d-%d-%d.par", run, prod, subiter);
  TString bcFileIn   = Form("rootfiles/bc-%d-%d-%d.root", run, prod, subiter);

  TFile *f = new TFile(rootFileIn.Data());
  TNtuple *vtxhits = (TNtuple *)f->Get("vtxhits");

  TFile *bcf = new TFile(bcFileIn.Data());
  TGraphErrors *gbc = (TGraphErrors *)bcf->Get("gbc");
  assert(gbc);

  SvxTGeo *geo = VTXModel(pisaFileIn.Data());

  geoEvents events;
  int nevents = 1000;
  GetEventsFromTree(vtxhits, geo, events, nevents);
  
  FitTracks(events, gbc);
  // FitTracks(events, 0);

  TF1 *linfit = new TF1("linfit", "[0] + [1]*x", -20, 20);
  linfit->SetLineColor(kBlack);
  linfit->SetNpx(500);
  TCanvas *c = DrawXY(geo, "trk_fits", "", "faint");
  c->Draw();

  TGraphErrors *g = new TGraphErrors(4);
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerColor(kRed);

  TRandom3 ran;
  ran.SetSeed(0);
  for (int i=0; i<100; i++)
  {
    int ev = ran.Integer(nevents);
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (TMath::Abs(TMath::Sin(trk.phi0)) < 0.7) 
        continue;
      int arm = (trk.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.
      double bcx = gbc->GetX()[arm];
      double bcy = gbc->GetY()[arm];

      // Reset graph
      while (g->GetN() > 0)
        g->RemovePoint(g->GetN()-1);

      g->SetPoint(0, bcx, bcy);
      for (int ihit=0; ihit<trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);
        g->SetPoint(ihit+1, hit.x, hit.y);
      }

      if (arm==0)
        linfit->SetRange(-20., bcx);
      if (arm==1)
        linfit->SetRange(bcx, +20.);
      linfit->SetParameters(trk.vy, TMath::Tan(trk.phi0));
  
      linfit->Draw("same");
      g->Draw("p,same");
      ltx.SetText(0.1, 0.92, Form("Event %d track %d", ev, t));
      ltx.Draw();
      gPad->Update();
      gSystem->Sleep(800);

    }
  }


  return;
}