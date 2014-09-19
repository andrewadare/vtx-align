#include "VtxAlignBase.h"
#include "BadLadders.h"
#include "VtxVis.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

const int nlayers = 4;
const int narms = 2; // 0=east, 1=west
const int nladders = 24;
const int ntrees = 2;

void FillHists(TH1D *hs[nlayers][nladders][ntrees],
               TH1D *hz[nlayers][nladders][ntrees],
               TH1D *hls[nlayers][narms][ntrees],
               TH1D *hlz[nlayers][narms][ntrees],
               TFile *f, const char *treename, int stage);

void SetupHist(TH1D *h, int stage);
void ModifyPad(TVirtualPad *pad, TH1D *h1, TH1D *h2, TString coord);

void DrawResults(int run = 411768,
                 int prod1 = 0,
                 int subit1 = 0,
                 int prod2 = 0,
                 int subit2 = 1,
                 const char *treename = "vtxhits")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int nLadders[4] = {10,20,16,24}; // ladders/layer

  TLatex ltx;
  ltx.SetNDC();
  TObjArray *cList = new TObjArray();

  TString rootFileIn1 = Form("rootfiles/%d-%d-%d.root", run, prod1, subit1);
  TString rootFileIn2 = Form("rootfiles/%d-%d-%d.root", run, prod2, subit2);
  TFile *inFile1  = new TFile(rootFileIn1.Data(), "read");
  TFile *inFile2  = new TFile(rootFileIn2.Data(), "read");
  assert(inFile1);
  assert(inFile2);


  TGraphErrors *gdead[nlayers];
  for (int lyr=0; lyr<nlayers; lyr++)
  {
    gdead[lyr] = new TGraphErrors();
    SetGraphProps(gdead[lyr],kNone, kOrange-8, kNone,kDot);
    gdead[lyr]->SetFillColorAlpha(kOrange-8, 0.5);
    gdead[lyr]->SetLineWidth(0);
    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
    {
      double dx = Dead(lyr,ldr) ? 0.5 : 0.0;
      double dy = Dead(lyr,ldr) ? 0.089 : 0.0;
      gdead[lyr]->SetPoint(ldr, ldr, 0.0);
      gdead[lyr]->SetPointError(ldr, dx, dy);
    }
  }

  TH1D *hds[nlayers][ntrees];
  TH1D *hdz[nlayers][ntrees];
  TH1D *hs[nlayers][nladders][ntrees];
  TH1D *hz[nlayers][nladders][ntrees];
  for (int lyr=0; lyr<nlayers; lyr++)
    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      for (int stage=0; stage<ntrees; stage++)
      {
        int nbins   = 100;
        double xmax = 0.12;
        hs[lyr][ldr][stage] = new TH1D(Form("hs%d%d%d", lyr, ldr, stage),
                                       Form("B%dL%d", lyr, ldr),
                                       nbins, -xmax, xmax);
        hz[lyr][ldr][stage] = new TH1D(Form("hz%d%d%d", lyr, ldr, stage),
                                       Form("B%dL%d", lyr, ldr),
                                       nbins, -xmax, xmax);
      }

  // Residuals for half-layer units
  const char *armlabel[narms] = {"E","W"};
  TH1D *hls[nlayers][narms][ntrees];
  TH1D *hlz[nlayers][narms][ntrees];
  for (int lyr=0; lyr<nlayers; lyr++)
    for (int arm=0; arm<narms; arm++)
      for (int stage=0; stage<ntrees; stage++)
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

  FillHists(hs, hz, hls, hlz, inFile1, treename, 0);
  FillHists(hs, hz, hls, hlz, inFile2, treename, 1);

  for (int lyr=0; lyr<nlayers; lyr++)
    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      for (int stage=0; stage<ntrees; stage++)
      {
        SetupHist(hs[lyr][ldr][stage], stage);
        SetupHist(hz[lyr][ldr][stage], stage);
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
      SetupHist(h1, 0);
      SetupHist(h2, 1);
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
      SetupHist(h1, 0);
      SetupHist(h2, 1);
      ModifyPad(gPad, h1, h2, "dz");
    }
  cList->Add(clz);

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
    TCanvas *cz = new TCanvas(Form("cz%d", lyr), Form("dz layer %d", lyr),
                              pw*nx, ph*ny);
    cs->Divide(nx, ny, 0.001, 0.001);
    cz->Divide(nx, ny, 0.001, 0.001);

    for (int ldr=0; ldr<nl; ldr++)
    {
      cs->cd(ldr+1);
      SetYMax(hs[lyr][ldr][0], hs[lyr][ldr][1]);
      for (int stage=0; stage<ntrees; stage++)
        hs[lyr][ldr][stage]->Draw(stage==0?"":"same");
      ModifyPad(gPad, hs[lyr][ldr][0], hs[lyr][ldr][1], "ds");

      cz->cd(ldr+1);
      SetYMax(hz[lyr][ldr][0], hz[lyr][ldr][1]);
      for (int stage=0; stage<ntrees; stage++)
        hz[lyr][ldr][stage]->Draw(stage==0?"":"same");
      ModifyPad(gPad, hz[lyr][ldr][0], hz[lyr][ldr][1], "dz");

      if (Dead(lyr,ldr))
      {
        TLatex l;
        l.SetTextColor(kGray+2);
        l.SetTextSize(0.4);
        l.SetNDC();
        cs->cd(ldr+1);
        l.DrawLatex(0.5, 0.5, "#times");
        cz->cd(ldr+1);
        l.DrawLatex(0.5, 0.5, "#times");
      }

    }
    cList->Add((TCanvas *)cs);
    cList->Add((TCanvas *)cz);
  }


  for (int lyr=0; lyr<4; lyr++)
  {
    for (int stage=0; stage<ntrees; stage++)
    {
      hdz[lyr][stage] = new TH1D(Form("hdz%d_%d",lyr, stage),
                                 Form("z-resid layer %d", lyr),
                                 nLadders[lyr], -0.5, nLadders[lyr] - 0.5);
      hds[lyr][stage] = new TH1D(Form("hds%d_%d",lyr, stage),
                                 Form("s-resid layer %d", lyr),
                                 nLadders[lyr], -0.5, nLadders[lyr] - 0.5);
      SetHistProps(hds[lyr][stage], kRed, kNone, kRed, kFullCircle, 1.5);
      SetHistProps(hdz[lyr][stage], kBlue, kNone, kBlue, kFullCircle, 1.5);
      if (stage==0)
      {
        SetHistProps(hds[lyr][stage], kGray, kNone, kGray, kFullCircle, 1.5);
        SetHistProps(hdz[lyr][stage], kGray, kNone, kGray, kFullCircle, 1.5);
        hds[lyr][0]->SetLineWidth(2);
        hdz[lyr][0]->SetLineWidth(2);
      }

      hds[lyr][stage]->SetXTitle("Ladder");
      hdz[lyr][stage]->SetXTitle("Ladder");

      hds[lyr][stage]->SetYTitle("#Deltas [cm]");
      hdz[lyr][stage]->SetYTitle("#Deltaz [cm]");

      hds[lyr][stage]->GetYaxis()->SetRangeUser(-0.045, 0.045);
      hdz[lyr][stage]->GetYaxis()->SetRangeUser(-0.045, 0.045);

      hds[lyr][stage]->GetXaxis()->SetNdivisions(210);
      hdz[lyr][stage]->GetXaxis()->SetNdivisions(210);

      for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      {
        hds[lyr][stage]->SetBinContent(ldr+1, hs[lyr][ldr][stage]->GetMean());
        hdz[lyr][stage]->SetBinContent(ldr+1, hz[lyr][ldr][stage]->GetMean());
        hds[lyr][stage]->SetBinError(ldr+1, hs[lyr][ldr][stage]->GetRMS());
        hdz[lyr][stage]->SetBinError(ldr+1, hz[lyr][ldr][stage]->GetRMS());
      }
    }
    ltx.SetTextSize(0.07);
    DrawObject(hds[lyr][0], "e0p", Form("ds_lyr%d",lyr), cList);
    gdead[lyr]->Draw("e5p,same");
    hds[lyr][1]->Draw("e0p,same");
    gPad->RedrawAxis();
    ltx.DrawLatex(0.23, 0.92, hds[lyr][0]->GetTitle());
    DrawObject(hdz[lyr][0], "e0p", Form("dz_lyr%d",lyr), cList);
    gdead[lyr]->Draw("e5p,same");
    hdz[lyr][1]->Draw("e0p,same");
    gPad->RedrawAxis();
    ltx.DrawLatex(0.23, 0.92, hdz[lyr][0]->GetTitle());
  }
  ltx.SetTextSize(0.06);
  ltx.SetTextFont(42);

  PrintPDFs(cList, Form("pdfs/run%d-pro%dsub%d-vs-pro%dsub%d",
                        run, prod1, subit1, prod2, subit2), "");
  PrintPDF(cList, Form("pdfs/run%d-pro%dsub%d-vs-pro%dsub%d",
                       run, prod1, subit1, prod2, subit2), "");

  gStyle->SetOptTitle(1);
  const char *parfile1 = Form("geom/%d-%d-%d.par", run, prod1, subit1);
  const char *parfile2 = Form("geom/%d-%d-%d.par", run, prod2, subit2);
  gROOT->Macro(Form("DiffGeometry.C(\"%s\", \"%s\")", parfile1, parfile2));

  return;
}

void
SetupHist(TH1D *h, int stage)
{
  if (stage==0)
  {
    h->SetLineColor(kGray+1);
    h->SetFillColor(kGray);
  }
  else if (TString(h->GetName()).Contains("hs"))
    h->SetLineColor(kGreen+3);
  else if (TString(h->GetName()).Contains("hz"))
    h->SetLineColor(kAzure+2);
  else if (TString(h->GetName()).Contains("hls"))
    h->SetLineColor(kGreen+3);
  else if (TString(h->GetName()).Contains("hlz"))
    h->SetLineColor(kAzure+2);

  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();
  ay->SetRangeUser(0, 1.3*h->GetMaximum());
  ax->SetNdivisions(205);
  ax->SetLabelSize(0.06);
  ay->SetNdivisions(205);
  ay->SetLabelSize(0.05);
  ax->SetTitleSize(0);
  return;
}

void
FillHists(TH1D *hs[nlayers][nladders][ntrees],
          TH1D *hz[nlayers][nladders][ntrees],
          TH1D *hls[nlayers][narms][ntrees],
          TH1D *hlz[nlayers][narms][ntrees],
          TFile *f, const char *treename, int stage)
{
  TTreeReader r(treename, f);
  TTreeReaderValue<float> ds(r, "res_s");
  TTreeReaderValue<float> dz(r, "res_z");
  TTreeReaderValue<float> layer(r, "layer");
  TTreeReaderValue<float> ladder(r, "ladder");
  TTreeReaderValue<float> gx(r, "gx");

  while (r.Next())
  {
    int lyr = *layer;
    int ldr = *ladder;
    float s = *ds;
    float z = *dz;
    float x = *gx;
    int arm = (x < 0.) ? 0 : 1; // 0 = East, 1 = West.

    // Printf("%s %d %d %d %f %f %x",
    //        treename, lyr, ldr, stage, s, z, hs[lyr][ldr][stage]);
    if (s>-0.2 && s<0.2 && z>-0.2 && z<0.2)
    {
      hs[lyr][ldr][stage]->Fill(s);
      hz[lyr][ldr][stage]->Fill(z);
      hls[lyr][arm][stage]->Fill(s);
      hlz[lyr][arm][stage]->Fill(z);
    }
  }

  return;
}

void
ModifyPad(TVirtualPad *pad, TH1D *h1, TH1D *h2, TString coord)
{
  pad->SetBottomMargin(0.15);
  pad->SetLeftMargin(0.2);
  pad->SetRightMargin(0.01);
  pad->SetTopMargin(0.01);

  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextSize(0.07);
  ltx.DrawLatex(0.23, 0.92, h1->GetTitle());
  ltx.DrawLatex(0.23, 0.85, Form("%s [cm]", coord.Data()));
  ltx.SetTextSize(0.06);
  ltx.DrawLatex(0.62,0.93, Form("#color[%d]{%.0f (%.0f) #mum}",
                                kBlack /*h1->GetLineColor()*/,
                                1e4*h1->GetMean(),
                                1e4*h1->GetRMS()));
  ltx.DrawLatex(0.62,0.87, Form("#color[%d]{%.0f (%.0f) #mum}",
                                h2->GetLineColor(),
                                1e4*h2->GetMean(),
                                1e4*h2->GetRMS()));
  return;
}
