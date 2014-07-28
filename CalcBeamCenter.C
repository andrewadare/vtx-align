#include "VtxAlignBase.h"

#include "GLSFitter.h"
#include "VertexFinder.h"
#include "DcaFunctions.h"
#include "BeamCenterFuncs.h"
#include "ParameterDefs.h"
#include "VtxIO.h"


void CalcBeamCenter(int run = 411768, 
                    int prod = 0,
                    int subiter = 1)
{
  TString rootFileIn  = Form("rootfiles/%d-%d-%d.root", run, prod, subiter);
  TString pisaFileIn  = Form("geom/%d-%d-%d.par", run, prod, subiter);

  TFile *f = new TFile(rootFileIn.Data());
  TNtuple *t = (TNtuple *)f->Get("vtxhits");

  gStyle->SetOptStat(0);
  gStyle->SetPalette(56, 0, 0.5); // Inverted "radiator", 50% transparency
  TObjArray *cList = new TObjArray();
  TLatex ltx;
  ltx.SetNDC();

  SvxTGeo *geo = VTXModel(pisaFileIn.Data());

  geoEvents events;
  GetEventsFromTree(t, geo, events, -1);
  FitTracks(events);

  // Event multiplicity histograms
  TH1D *hne = new TH1D("hne", "East;tracks/event", 100, 0, 100);
  TH1D *hnw = new TH1D("hnw", "West;tracks/event", 100, 0, 100);

  Printf("Filling multiplicity distributions...");
  int minmult = 20;
  int N = 0; // Pre-allocation size
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ne=0, nw=0;
    int mult = events[ev].size();
    if (mult >= minmult)
      N++;
    for (int t=0; t<mult; t++)
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

  Printf("Allocating vertex arrays...");
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
  double eastrange = TMath::Max(qxe[nq-1]-qxe[0], qye[nq-1]-qye[0]);
  double westrange = TMath::Max(qxw[nq-1]-qxw[0], qyw[nq-1]-qyw[0]);

  double x0 = qxe[2] - eastrange/2;
  double y0 = qye[2] - eastrange/2;
  double x1 = qxe[2] + eastrange/2;
  double y1 = qye[2] + eastrange/2;
  TH2D *hve = new TH2D("hve", Form("East vertex - prod %d sub %d;"
                       "x [cm];y [cm]", prod, subiter),
                       100, x0, x1, 100, y0, y1);
  x0 = qxw[2] - westrange/2;
  y0 = qyw[2] - westrange/2;
  x1 = qxw[2] + westrange/2;
  y1 = qyw[2] + westrange/2;
  TH2D *hvw = new TH2D("hvw", Form("West vertex - prod %d sub %d;"
                       "x [cm];y [cm]", prod, subiter),
                       100, x0, x1, 100, y0, y1);
  hve->FillN(nfilled, vxe, vye, NULL, 1);
  hvw->FillN(nfilled, vxw, vyw, NULL, 1);

  // West - East offset histograms
  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k=0; k<3; k++)
    hdv[k] = new TH1D(Form("hdv%d",k),
                      Form("West - East vertex difference #Delta%s "
                           "- prod %d sub %d; #Delta%s [cm];events", 
                           xyzstr[k], prod, subiter, xyzstr[k]),
                      200, -1, 1);
  hdv[0]->FillN(nfilled, dvx, NULL, 1);
  hdv[1]->FillN(nfilled, dvy, NULL, 1);
  hdv[2]->FillN(nfilled, dvz, NULL, 1);

  x0=-0.05;
  y0=-0.05;
  x1=+0.05;
  y1=+0.05;
  TH2D *hxye = new TH2D("hxye", "DCA (east);x [cm];y [cm]",
                        200, x0, x1, 200, y0, y1);
  TH2D *hxyw = new TH2D("hxyw", "DCA (west);x [cm];y [cm]",
                        200, x0, x1, 200, y0, y1);
  TH1D *hre = new TH1D("hre", "hre", 100, 0, 0.5);
  TH1D *hrw = new TH1D("hrw", "hrw", 100, 0, 0.5);
  hxye->SetLineColorAlpha(kRed+2, 0.4);
  hxyw->SetLineColorAlpha(kRed+2, 0.4);
  DcaDist(events, "east", hxye, hre);
  DcaDist(events, "west", hxyw, hrw);

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

  TCanvas *c = new TCanvas("c", "dca2d", 1000, 500);
  c->Divide(2,1,0,0);
  c->cd(1);
  TAxis *ax = hxye->GetXaxis();
  TAxis *ay = hxye->GetYaxis();
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(208);
  ay->SetNdivisions(208);
  hxye->Draw("cont3");
  c->cd(2);
  ax = hxyw->GetXaxis();
  ay = hxyw->GetYaxis();
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(208);
  ay->SetNdivisions(208);
  hxyw->Draw("cont3");
  cList->Add(c);

  TCanvas *cr = new TCanvas("cr", "dcar", 500, 500);
  SetHistProps(hre,kAzure+1,kNone,kAzure+1);
  SetHistProps(hrw,kOrange,kNone,kOrange);
  hre->SetLineWidth(2);
  hrw->SetLineWidth(2);
  hrw->SetTitle("Beam center DCA #color[861]{ east}, #color[800]{ west};DCA [cm];");
  hrw->Draw("");
  hre->Draw("same");
  cr->SetLogy();
  cList->Add(cr);

  DrawObject(hbp, "colz", "bcdcavsphi", cList, 900, 500);
  hbp->SetTitle("DCA to beam center");
  ax = hbp->GetXaxis();
  ay = hbp->GetYaxis();
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
         1.0 - hre->Integral(1,hre->FindBin(0.099))/hre->Integral(),
         1.0 - hrw->Integral(1,hrw->FindBin(0.099))/hrw->Integral());

  PrintPDFs(cList, Form("pdfs/prod%d/subiter%d", prod, subiter), "");
  PrintPDF(cList, Form("pdfs/beam-center-pro%d-sub%d", prod, subiter));

  return;
}
