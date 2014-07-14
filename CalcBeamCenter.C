#include "UtilFns.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "DcaFunctions.h"
#include "BeamCenterFuncs.h"
#include "ParameterDefs.h"
#include "GLSFitter.h"
#include "VtxIO.h"

#include <TEllipse.h>

// Solve the linear system y0[i] = -m[i]*bc0 + bc1 for i..ntrk-1 tracks.
// m[i] is the slope (tan(phi)), y0[i] is the y-intercept, and (bc0,bc1)
// is the least-squares beam (x,y) position.

void CalcBeamCenter(int run = 411768, int iter = 1)
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

  // DCA to beam center vs phi
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

  DrawObject(hbp, "colz", "bcdcavsphi", cList, 900, 500);
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
