#include "UtilFns.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "DcaFunctions.h"
#include "BeamCenterFuncs.h"
#include "VertexFinder.h"
#include "ParameterDefs.h"
#include "GLSFitter.h"
#include "VtxIO.h"

#include <TEllipse.h>
#include <TLatex.h>

// Solve the linear system y0[i] = -m[i]*bc0 + bc1 for i..ntrk-1 tracks.
// m[i] is the slope (tan(phi)), y0[i] is the y-intercept, and (bc0,bc1)
// is the least-squares beam (x,y) position.

void CalcBeamCenter(int run = 411768, int iter = 1)
{
  int nbc = 10000;
  TLatex ltx;
  ltx.SetNDC();

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

  geoEvents events;
  GetEventsFromTree(t, geo, events, -1);
  const int N = events.size();
  FitTracks(events);

  // Event multiplicity histograms
  TH1D *hne = new TH1D("hne", "East;tracks/event", 100, 0, 100);
  TH1D *hnw = new TH1D("hnw", "West;tracks/event", 100, 0, 100);

  Printf("Filling multiplicity distributions...");
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

  int minmult = 20;
  int nfilled = 0;
  double vxe[N], vye[N], vze[N];  // Vertex from east arm
  double vxw[N], vyw[N], vzw[N];  // Vertex from west arm
  double dvx[N], dvy[N], dvz[N];  // West - east offsets

  Printf("Filling vertex arrays...");
  FillVertexArrays(events, minmult, vxe, vye, vze, vxw, vyw, vzw, dvx, dvy, dvz,
                   nfilled);

  Printf("Computing quantiles for {e,w}{x,y}");
  const int nq = 5;
  double probs[nq] = {0.01, 0.32, 0.50, 0.68, 0.99};
  double qxe[nq] = {0}, qxw[nq] = {0}, qye[nq] = {0}, qyw[nq] = {0};
  TMath::Quantiles(nfilled, nq, vxe, qxe, probs, false);
  TMath::Quantiles(nfilled, nq, vxw, qxw, probs, false);
  TMath::Quantiles(nfilled, nq, vye, qye, probs, false);
  TMath::Quantiles(nfilled, nq, vyw, qyw, probs, false);

  TVectorD bce(2); bce(0) = qxe[2]; bce(1) = qye[2];
  TVectorD bcw(2); bcw(0) = qxw[2]; bcw(1) = qyw[2];

  Printf("East-west offsets...");
  TH2D *hve = new TH2D("hve", Form("East vertex - iteration %d;x [cm];y [cm]",iter),
                       100, qxe[0], qxe[4], 100, qye[0], qye[4]);
  TH2D *hvw = new TH2D("hvw", Form("West vertex - iteration %d;x [cm];y [cm]",iter),
                       100, qxw[0], qxw[4], 100, qyw[0], qyw[4]);

  hve->FillN(nfilled, vxe, vye, NULL, 1);
  hvw->FillN(nfilled, vxw, vyw, NULL, 1);

  // West - East offset histograms
  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k=0; k<3; k++)
    hdv[k] = new TH1D(Form("hdv%d",k),
                      Form("West - East vertex difference #Delta%s;"
                           "#Delta%s [cm];events",
                           xyzstr[k], xyzstr[k]),
                      200, -1, 1);
  hdv[0]->FillN(nfilled, dvx, NULL, 1);
  hdv[1]->FillN(nfilled, dvy, NULL, 1);
  hdv[2]->FillN(nfilled, dvz, NULL, 1);

  // Printf("Computing east arm beam center...");
  // TVectorD bce  = BeamCenter(events, nbc, "east", "print");
  TH1D *he = new TH1D("he", "he", 100, 0, 0.5);
  TGraph *ge = DcaDist(events,bce,"east",he);

  // Printf("Computing west arm beam center...");
  // TVectorD bcw = BeamCenter(events, nbc, "west", "print");
  TH1D *hw = new TH1D("hw", "hw", 100, 0, 0.5);
  TGraph *gw = DcaDist(events,bcw,"west",hw);

  Printf("DCA to beam center vs phi");
  double dmin = 0.0, dmax = 0.1;
  TH2D *hbp = new TH2D("hbp", "", 100, 0, TMath::TwoPi(), 100, dmin, dmax);
  TProfile *bprof = new TProfile("bprof", "", 100, 0, TMath::TwoPi(), dmin, dmax);
  bprof->SetMarkerStyle(kFullCircle);

  // Error calculation options:
  // "" : error = sigma/sqrt(N) = std. dev of mean.
  // "s": error = sigma = std. dev of distribution in dmin, dmax.
  // "g": error = 1/sqrt(sum(w)) where w is (weighted) sum of entries.
  bprof->BuildOptions(dmin, dmax, "g");
  DcaVsPhi(events, bce, bcw, hbp, bprof);

  Printf("DCA vs phi...");
  TH2D *hdp = new TH2D("hdp", "", 100, 0, TMath::TwoPi(), 100, dmin, dmax);
  TProfile *dprof = new TProfile("dprof", "", 100, 0, TMath::TwoPi(), dmin, dmax);
  dprof->SetMarkerStyle(kFullCircle);
  dprof->BuildOptions(dmin, dmax, "g");

  // This is inefficient, repeating all this work...improve later
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");
      DcaVsPhi(events[ev], ve, vw, hdp, dprof);
    }
  }

  // ---- ---- ---- ---- ---- ---- ----
  //                Plot
  // ---- ---- ---- ---- ---- ---- ----
  Printf("Plotting...");
  DrawObject(hne, "", "east_mult", cList);
  gPad->SetLogx();
  gPad->SetLogy();
  DrawObject(hnw, "", "west_mult", cList);
  gPad->SetLogx();
  gPad->SetLogy();

  TCanvas *cve = new TCanvas("cve", "cve", 500, 500);
  hve->Draw("col");
  ltx.DrawLatex(0.15, 0.85, Form("BC = %.3f, %.3f", bce(0), bce(1)));
  cList->Add(cve);
  TCanvas *cvw = new TCanvas("cvw", "cvw", 500, 500);
  hvw->Draw("col");
  ltx.DrawLatex(0.15, 0.85, Form("BC = %.3f, %.3f", bcw(0), bcw(1)));
  cList->Add(cvw);

  for (int k=0; k<3; k++)
  {
    DrawObject(hdv[k], "", Form("cdv%d",k), cList);
    ltx.DrawLatex(0.15, 0.85, Form("Mean %.3f", hdv[k]->GetMean()));
    ltx.DrawLatex(0.15, 0.80, Form("Std dev %.3f", hdv[k]->GetRMS()));
  }

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
  bprof->Draw("same");

  DrawObject(hdp, "colz", "dcavsphi", cList, 900, 500);
  hdp->SetTitle("Zero-field DCA");
  ax = hdp->GetXaxis();
  ay = hdp->GetYaxis();
  ax->SetTitle("#Leftarrow  west      #phi + #pi/2      east  #Rightarrow");
  ay->SetTitle("DCA [cm]");
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(408);
  ay->SetNdivisions(205);
  gPad->SetGridy();
  dprof->Draw("same");


  Printf("Fraction outside 1000um: %.3f (e) %.3f (w)",
         1.0 - he->Integral(1,he->FindBin(0.099))/he->Integral(),
         1.0 - hw->Integral(1,hw->FindBin(0.099))/hw->Integral());

  PrintPDFs(cList, Form("pdfs/nonideal/iter%d", iter), "");
  PrintPDF(cList, Form("pdfs/nonideal/beam-center-iter%d", iter));

  return;
}
