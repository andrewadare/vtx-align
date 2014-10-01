#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"
#include "VtxIO.h"
#include "VtxVis.h"

void ResidPlots(geoEvents &ev1, geoEvents &ev2, TObjArray *cList=0);
void ModifyPad(TVirtualPad *pad, TH1D *h1, TH1D *h2, TString coord = "ds");
double MeanPhi(SvxGeoTrack &trk, double *bc = 0);

void TrackInfo(int run = 411768,
               int prod = 0,
               int subiter = 1)
{
  bool write = true;
  gStyle->SetOptStat(0);

  TObjArray *cList = new TObjArray();

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

  geoEvents vtx_evts;
  GetEventsFromTree(vtxhits, geo, vtx_evts, 100000);
  FitTracks(vtx_evts, gbc);

  geoEvents evs_nobc = vtx_evts;
  FitTracks(evs_nobc);

  TH1D *phic[4][2];
  for (int lyr=0; lyr<4; lyr++)
    for (int stage=0; stage<2; stage++)
    {
      phic[lyr][stage] = new TH1D(Form("phic%d%d",lyr,stage),"",
                                  200, -0.08, 0.08);
      TH1D *h = phic[lyr][stage];
      h->SetTitle(Form("Layer %d", lyr));
      h->SetLineColor(stage ? kRed : kBlack);
      TAxis *ax = h->GetXaxis();
      TAxis *ay = h->GetYaxis();
      ay->SetNdivisions(205);
      ay->SetLabelSize(0.05);

      ax->SetNdivisions(205);
      ax->SetTitle("#phi_{c} - #LT#phi_{c}#GT [cm]");
      ax->SetLabelSize(0.06);
      ax->SetTitleSize(0.06);
      ax->CenterTitle();
    }

  TH1D *hpol = new TH1D("hpol", "Track #theta distribution", 100, 0., TMath::Pi());
  TH1D *hphi = new TH1D("hphi", "Track #phi distribution", 100, 0., TMath::TwoPi());
  TH1D *hphibc = new TH1D("hphibc", "Track #phi distribution", 100, 0., TMath::TwoPi());
  TH1D *hz0 = new TH1D("hz0", "Track z-vertex distribution", 100, -12, 12);
  TH1D *hzslope = new TH1D("hzslope", "Track z vs r slope", 100, -2, 2);
  TH1D *hxys_bc = new TH1D("hxys_bc", "Track y vs x slope", 100, -3, 3);
  TH1D *hxys_no = new TH1D("hxys_no", "Track y vs x slope", 100, -3, 3);
  hxys_no->SetLineColor(kBlack);
  hxys_bc->SetLineColor(kRed);

  TH1D *hy0 = new TH1D("hy0", "Track y-intercept distribution", 100, -1, 1);
  TH1D *hy0bc = new TH1D("hy0bc", "Track y-intercept distribution", 100, -1, 1);
  hy0->SetLineColor(kBlack);
  hy0bc->SetLineColor(kRed);

  hphi->SetLineColor(kBlack);
  hphibc->SetLineColor(kRed);

  // Track loop: no beam center used in fits
  for (unsigned int ev=0; ev<evs_nobc.size(); ev++)
    for (unsigned int t=0; t<evs_nobc[ev].size(); t++)
    {
      SvxGeoTrack trk = evs_nobc[ev][t];
      hphi->Fill(trk.phi0);
      hxys_no->Fill(TMath::Tan(trk.phi0));
      hy0->Fill(trk.yp0);

      double mphi = MeanPhi(trk);
      for (int ihit=0; ihit<trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);
        double phi = TMath::ATan2(hit.y, hit.x);
        double dphi = fmod(phi - mphi, TMath::TwoPi());
        phic[hit.layer][0]->Fill(dphi);
      }

    }

  // Track loop: beam center included
  for (unsigned int ev=0; ev<vtx_evts.size(); ev++)
    for (unsigned int t=0; t<vtx_evts[ev].size(); t++)
    {
      SvxGeoTrack trk = vtx_evts[ev][t];
      hpol->Fill(trk.the0);
      hphibc->Fill(trk.phi0);
      hz0->Fill(trk.vz);
      hzslope->Fill(TMath::Tan(trk.the0 - TMath::PiOver2()));
      hxys_bc->Fill(TMath::Tan(trk.phi0));
      hy0bc->Fill(trk.yp0);

      int arm = (trk.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.
      double bc[2] = {gbc->GetX()[arm], gbc->GetY()[arm]};
      double mphi = MeanPhi(trk, bc);
      for (int ihit=0; ihit<trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);
        double phi = TMath::ATan2(hit.y - bc[1], hit.x - bc[0]);
        double dphi = fmod(phi - mphi, TMath::TwoPi());
        phic[hit.layer][1]->Fill(dphi);
      }
    }

  // ---- ---- ---- ---- ---- ---- ----
  //                Plot
  // ---- ---- ---- ---- ---- ---- ----
  DrawObject(hpol, "", "theta_dist", cList);
  DrawObject(hz0, "", "z0_dist", cList);
  DrawObject(hphi, "", "trk_phi_dist", cList);
  hphibc->Draw("same");
  DrawObject(hzslope, "", "zvsr_slope_dist", cList);
  DrawObject(hxys_no, "", "az_slope_dist", cList);
  hxys_bc->Draw("same");
  DrawObject(hy0bc, "", "y0_dist", cList);
  hy0->Draw("same");

  ResidPlots(vtx_evts, evs_nobc, cList);

  TCanvas *cphi = new TCanvas("cphi", "Cluster phi dispersion", 1200, 400);
  cphi->Divide(4, 1, 0.001, 0.001);
  for (int lyr=0; lyr<4; lyr++)
  {
    cphi->cd(lyr+1);
    phic[lyr][1]->Draw();
    phic[lyr][0]->Draw("same");
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextSize(0.07);
    ltx.SetTextFont(42);
    ltx.DrawLatex(0.23, 0.91, Form("stdev %.0f #mum",
                                   1e4*phic[lyr][1]->GetRMS()));
  }
  cList->Add(cphi);

  if (write)
  {
    PrintPDFs(cList, Form("pdfs"), "");
    PrintPDF(cList, Form("pdfs/trackinfo-%d-%d", prod, subiter));
  }

  return;
}

double
MeanPhi(SvxGeoTrack &trk, double *bc)
{
  double meanphi = 0.;
  for (int ihit=0; ihit<trk.nhits; ihit++)
  {
    SvxGeoHit hit = trk.GetHit(ihit);
    if (bc)
      meanphi += 1./trk.nhits * TMath::ATan2(hit.y - bc[1], hit.x - bc[0]);
    else
      meanphi += 1./trk.nhits * TMath::ATan2(hit.y, hit.x);
  }
  meanphi = fmod(meanphi, TMath::TwoPi());
  if (meanphi < 0.)
    meanphi += TMath::TwoPi();
  return meanphi;
}

void
ResidPlots(geoEvents &ev1, geoEvents &ev2, TObjArray *cList)
{
  const int nlayers = 4;
  const int nladders = 24;
  const int ntrees = 2;
  const int narms = 2; // 0=east, 1=west
  int nLadders[4] = {10,20,16,24}; // ladders/layer

  // Create histograms
  TH1D *hr[nlayers][nladders][2];
  TH1D *hls[nlayers][narms][2];
  TH1D *hlz[nlayers][narms][2];

  // Residuals for individual VTX ladders
  for (int lyr=0; lyr<nlayers; lyr++)
    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      for (int stage=0; stage<2; stage++)
      {
        int nbins   = 200;
        double xmax = 0.06;
        hr[lyr][ldr][stage] = new TH1D(Form("hr%d%d%d", lyr, ldr, stage),
                                       Form("B%dL%d", lyr, ldr),
                                       nbins, -xmax, xmax);
      }

  // Residuals for half-layer units
  const char *armlabel[2] = {"E","W"};
  for (int lyr=0; lyr<nlayers; lyr++)
    for (int arm=0; arm<narms; arm++)
      for (int stage=0; stage<2; stage++)
      {
        int nbins   = 200;
        double xmax = 0.06;
        hls[lyr][arm][stage] = new TH1D(Form("hls%d%d%d", lyr, arm, stage),
                                        Form("B%d%s", lyr, armlabel[arm]),
                                        nbins, -xmax, xmax);
        hlz[lyr][arm][stage] = new TH1D(Form("hlz%d%d%d", lyr, arm, stage),
                                        Form("B%d%s", lyr, armlabel[arm]),
                                        nbins, -xmax, xmax);
      }


  // Fill residuals from ev1
  for (unsigned int ev=0; ev<ev1.size(); ev++)
    for (unsigned int t=0; t<ev1[ev].size(); t++)
    {
      SvxGeoTrack trk = ev1[ev][t];
      int arm = East(trk.phi0) ? 0 : 1;
      for (int ihit=0; ihit < trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);
        hr[hit.layer][hit.ladder][0]->Fill(hit.ds);
        if (trk.nhits == 4)
        {
          hls[hit.layer][arm][0]->Fill(hit.ds);
          hlz[hit.layer][arm][0]->Fill(hit.dz);
        }
      }
    }
  // Fill residuals from ev2
  for (unsigned int ev=0; ev<ev2.size(); ev++)
    for (unsigned int t=0; t<ev2[ev].size(); t++)
    {
      SvxGeoTrack trk = ev2[ev][t];
      int arm = East(trk.phi0) ? 0 : 1;
      for (int ihit=0; ihit < trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);
        hr[hit.layer][hit.ladder][1]->Fill(hit.ds);
        if (trk.nhits == 4)
        {
          hls[hit.layer][arm][1]->Fill(hit.ds);
          hlz[hit.layer][arm][1]->Fill(hit.dz);
        }
      }
    }

  // Draw individual residual distributions on sub-pads
  for (int lyr=0; lyr<4; lyr++)
  {
    Printf("Drawing residuals in layer %d", lyr);
    int nl = nLadders[lyr];
    int nx = lyr ? nl/4 : nl/2; // Number of pads along x
    int ny = lyr ? 4 : 2;       // Number of pads along y
    int ph = 200;               // pad height in px
    int pw = 200;               // pad width in px
    TCanvas *cs = new TCanvas(Form("cs%d", lyr), Form("ds layer %d", lyr),
                              pw*nx, ph*ny);
    cs->Divide(nx, ny, 0.001, 0.001);

    for (int ldr=0; ldr<nl; ldr++)
    {
      cs->cd(ldr+1);
      TH1D *h1 = hr[lyr][ldr][0];
      TH1D *h2 = hr[lyr][ldr][1];
      SetYMax(h1,h2);
      h1->Draw("");
      h2->Draw("same");
      ModifyPad(gPad, h1, h2);
    }
    if (cList)
      cList->Add(cs);
  }


  TCanvas *cls = new TCanvas(Form("cls%d", 0), Form("%d", 0), 1200, 800);
  cls->Divide(4, 2, 0.001, 0.001);
  for (int lyr=0; lyr<4; lyr++)
    for (int arm=0; arm<2; arm++)
    {
      cls->cd(arm*4 + lyr+1);
      TH1D *h1 = hls[lyr][arm][0];
      TH1D *h2 = hls[lyr][arm][1];
      SetYMax(h1,h2);
      h1->Draw("");
      h2->Draw("same");
      ModifyPad(gPad, h1, h2, "ds");
    }
  cList->Add(cls);

  TCanvas *clz = new TCanvas(Form("clz%d", 0), Form("%d", 0), 1200, 800);
  clz->Divide(4, 2, 0.001, 0.001);
  for (int lyr=0; lyr<4; lyr++)
    for (int arm=0; arm<2; arm++)
    {
      clz->cd(arm*4 + lyr+1);
      TH1D *h1 = hlz[lyr][arm][0];
      TH1D *h2 = hlz[lyr][arm][1];
      SetYMax(h1,h2);
      h1->Draw("");
      h2->Draw("same");
      ModifyPad(gPad, h1, h2, "dz");
    }
  cList->Add(clz);

  return;
}

void
ModifyPad(TVirtualPad *pad, TH1D *h1, TH1D *h2, TString coord)
{
  pad->SetBottomMargin(0.15);
  pad->SetLeftMargin(0.2);
  pad->SetRightMargin(0.01);
  pad->SetTopMargin(0.01);

  // Modify histograms too - not really appropriate here - don't care for now
  h1->SetLineColor(kGray+3);
  h2->SetLineColor(kRed+1);
  TAxis *ax = h1->GetXaxis();
  TAxis *ay = h1->GetYaxis();
  ay->SetRangeUser(0, 1.3*h1->GetMaximum());
  ax->SetNdivisions(205);
  ax->SetLabelSize(0.06);
  ay->SetNdivisions(205);
  ay->SetLabelSize(0.05);
  ax->SetTitleSize(0);

  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextSize(0.07);
  ltx.DrawLatex(0.23, 0.92, h1->GetTitle());
  ltx.DrawLatex(0.23, 0.85, Form("%s [cm]", coord.Data()));
  ltx.SetTextSize(0.06);
  ltx.DrawLatex(0.62,0.93, Form("#color[%d]{%.0f (%.0f) #mum}",
                                h1->GetLineColor(),
                                1e4*h1->GetMean(),
                                1e4*h1->GetRMS()));
  ltx.DrawLatex(0.62,0.87, Form("#color[%d]{%.0f (%.0f) #mum}",
                                h2->GetLineColor(),
                                1e4*h2->GetMean(),
                                1e4*h2->GetRMS()));
  return;
}

