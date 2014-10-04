#include "VtxAlignBase.h"
#include "BadLadders.h"
#include "VtxVis.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

int nLaddersPerLayer[4] = {10, 20, 16, 24}; // ladders/layer

void DrawHalfLayerResidPlots(int ntrees, TObjArray *cList = 0);
void DrawLadderResidPlots(int ntrees, TObjArray *cList = 0);
void DrawSummaryPlots(int ntrees, TObjArray *cList = 0);
void InitResidHists(const char *treename, int ntrees);
void SetupHist(TH1D *h, int stage);
void FillHists(TFile *f, const char *treename, int stage, TObjArray *cList = 0);
void ModifyPad(TVirtualPad *pad, TH1D *h1, TH1D *h2, TString coord);
bool CheckValue(ROOT::TTreeReaderValueBase *value);
TH1D *Hist(const char *prefix, int layer, int ladder_or_arm, int stage);
TH1D *Hist(const char *prefix, int layer, int stage);

// To plot variables from only one TTree, either:
// 1. set prod2 and subit2 to any int < 0, or
// 2. set prod2 = prod1, subit2 = subit1.
void DrawResults(int run = 123456,
                 int prod1 = 0,
                 int subit1 = 0,
                 int prod2 = -1,
                 int subit2 = -1,
                 const char *treename = "vtxtrks")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TObjArray *cList = new TObjArray();

  int ntrees = 2;
  if (prod1 == prod2 && subit1 == subit2)
    ntrees = 1;
  if (prod2 < 0 && subit2 < 0)
    ntrees = 1;

  vecs rootFileIn;
  rootFileIn.push_back(Form("rootfiles/%d-%d-%d.root", run, prod1, subit1));
  if (ntrees==2)
    rootFileIn.push_back(Form("rootfiles/%d-%d-%d.root", run, prod2, subit2));

  vector<TFile *> inFiles;
  for (int stage = 0; stage < ntrees; stage++)
  {
    inFiles.push_back(new TFile(rootFileIn[0].c_str(), "read"));
    assert(inFiles[stage]);
  }

  InitResidHists(treename, ntrees);
  for (int stage = 0; stage < ntrees; stage++)
  {
    FillHists(inFiles[stage], treename, stage, cList);
  }

  // TGraphErrors *gdead[nlayers];
  // MarkDeadLadders(gdead);
  DrawHalfLayerResidPlots(ntrees, cList);
  DrawLadderResidPlots(ntrees, cList);
  DrawSummaryPlots(ntrees, cList);

  // ltx.SetTextSize(0.06);
  // ltx.SetTextFont(42);

  PrintPDFs(cList, Form("pdfs/run%d-pro%dsub%d-vs-pro%dsub%d-%s",
                        run, prod1, subit1, prod2, subit2, treename), "");
  PrintPDF(cList, Form("pdfs/run%d-pro%dsub%d-vs-pro%dsub%d-%s",
                       run, prod1, subit1, prod2, subit2, treename), "");

  if (ntrees == 2)
  {
    gStyle->SetOptTitle(1);
    const char *parfile1 = Form("geom/%d-%d-%d.par", run, prod1, subit1);
    const char *parfile2 = Form("geom/%d-%d-%d.par", run, prod2, subit2);
    gROOT->Macro(Form("DiffGeometry.C(\"%s\", \"%s\")", parfile1, parfile2));
  }
  return;
}

void
FillHists(TFile *f, const char *treename, int stage, TObjArray *cList)
{
  TH1D *xydcae = new TH1D(Form("xydcae%d",stage), Form("%d",stage), 200, -0.05, 0.05);
  TH1D *xydcaw = new TH1D(Form("xydcaw%d",stage), Form("%d",stage), 200, -0.05, 0.05);
  TH1D *zdcae = new TH1D(Form("zdcae%d",stage), Form("%d",stage), 200, -0.05, 0.05);
  TH1D *zdcaw = new TH1D(Form("zdcaw%d",stage), Form("%d",stage), 200, -0.05, 0.05);
  TH1D *hvze = new TH1D(Form("hvze%d",stage), Form("%d",stage), 200, -5, 5);
  TH1D *hvzw = new TH1D(Form("hvzw%d",stage), Form("%d",stage), 200, -5, 5);
  TH2D *xydcaphi  = new TH2D(Form("xydcaphi%d",stage), ";#phi [rad];x-y DCA [cm]",
                             100, -TMath::PiOver2(), 3*TMath::PiOver2(),
                             100, -0.05, +0.05);
  TH2D *zdcatheta  = new TH2D(Form("zdcatheta%d",stage), ";#theta [rad];z DCA [cm]",
                              100, 0.0*TMath::Pi(), 1.0*TMath::Pi(),
                              100, -0.05, +0.05);

  // x-y vertex distributions
  int prod = 9999;
  int subiter = 9999;
  double x0 = -0.5, y0 = -0.5, x1 = +0.5, y1 = +0.5;
  TH2D *hve = new TH2D(Form("hve%d",stage),
                       Form("East vertex - prod %d step %d;"
                            "x [cm];y [cm]", prod, subiter),
                       200, x0, x1, 200, y0, y1);
  TH2D *hvw = new TH2D(Form("hvw%d",stage),
                       Form("West vertex - prod %d step %d;"
                            "x [cm];y [cm]", prod, subiter),
                       200, x0, x1, 200, y0, y1);

  // West - East offset histograms
  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k=0; k<3; k++)
    hdv[k] = new TH1D(Form("hdv_%d_%d",stage, k),
                      Form("West - East vertex difference #Delta%s "
                           "- prod %d step %d; #Delta%s [cm];events",
                           xyzstr[k], prod, subiter, xyzstr[k]),
                      200, -1, 1);

  TTreeReader r(treename, f);

  TTreeReaderValue<float> xydca(r, "xydca");
  TTreeReaderValue<float>  zdca(r, "zdca");
  TTreeReaderValue<float>   phi(r, "trkphi");
  TTreeReaderValue<float> theta(r, "trktheta");

  TTreeReaderArray<float>   ds(r, "res_s");
  TTreeReaderArray<float>   dz(r, "res_z");
  TTreeReaderArray<int>  layer(r, "layer");
  TTreeReaderArray<int> ladder(r, "ladder");
  TTreeReaderArray<float>   gx(r, "gx");
  TTreeReaderArray<float> vertex(r, "primvtx");

  Info("", "Filling histograms from %s", f->GetName());
  while (r.Next())
  {
    int arm = (gx[0] < 0.) ? 0 : 1; // 0 = East, 1 = West.

    // Hit-level variables
    for (int i = 0; i < (int)ds.GetSize(); ++i)
    {
      int lyr = layer[i];
      int ldr = ladder[i];
      float s = ds[i];
      float z = dz[i];
      float x = gx[i];

      // Printf("%s %d %d %d %f %f", treename, lyr, ldr, stage, s, z);

      if (s > -0.2 && s < 0.2 && z > -0.2 && z < 0.2)
      {
        // Residuals by ladder
        Hist("s",lyr,ldr,stage) ->Fill(s);
        Hist("z",lyr,ldr,stage) ->Fill(z);
        // Residuals by halflayer
        Hist("ls",lyr,arm,stage)->Fill(s);
        Hist("lz",lyr,arm,stage)->Fill(z);
      }

    }

    float phiwrap = *phi;
    if (phiwrap > 1.5*TMath::Pi())
      phiwrap -= TMath::TwoPi();
    // Track/event level variables
    TVectorD ve(3);
    TVectorD vw(3);
    if (arm == 0)
    {
      ve(0) = vertex[0];
      ve(1) = vertex[1];
      ve(2) = vertex[2];
      xydcae->Fill(*xydca);
      zdcae->Fill(*zdca);
      hve->Fill(vertex[0], vertex[1]);
    }
    if (arm == 1)
    {
      vw(0) = vertex[0];
      vw(1) = vertex[1];
      vw(2) = vertex[2];
      xydcaw->Fill(*xydca);
      zdcaw->Fill(*zdca);
      hvw->Fill(vertex[0], vertex[1]);
    }

    for (int k=0; k<3; k++)
      hdv[k]->Fill(vw(k) - ve(k));

    xydcaphi->Fill(phiwrap, *xydca);
    zdcatheta->Fill(*theta, *zdca);
  }


  // Plot--------------------------------------------------------
  Info("", "Plotting histograms from %s", f->GetName());
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);
  TCanvas *c = 0;
  c = new TCanvas(Form("cdca%d",stage),
                  Form("dca_east_west_%d",stage), 1000, 500);
  SetYMax(xydcae, xydcaw);
  c->Divide(2,1,0,0);
  c->cd(1);
  xydcae->Draw();
  ltx.DrawLatex(0.15, 0.90, Form("East: (Mean, Std Dev) = (%.0f, %.0f) #mum",
                                 1e4*xydcae->GetMean(), 1e4*xydcae->GetRMS()));
  c->cd(2);
  xydcaw->Draw();
  ltx.DrawLatex(0.15, 0.90, Form("West: (Mean, Std Dev) = (%.0f, %.0f) #mum",
                                 1e4*xydcaw->GetMean(), 1e4*xydcaw->GetRMS()));
  cList->Add(c);

  DrawObject(zdcae, "",  Form("zdcae_%d", stage), cList);
  DrawObject(zdcaw, "",  Form("zdcae_%d", stage), cList);

  DrawObject(xydcaphi, "col",  Form("xydca_vs_phi_%d",stage), cList);
  DrawObject(zdcatheta, "col", Form("zdca_vs_theta_%d",stage), cList);

  hve->GetXaxis()->SetRangeUser(hve->GetMean(1)-0.1,hve->GetMean(1)+0.1);
  hve->GetYaxis()->SetRangeUser(hve->GetMean(2)-0.1,hve->GetMean(2)+0.1);
  // DrawObject(hve, "col", Form("xy_vertex%d",stage), cList, 500, 500);
  // hvw->Draw("same");

  for (int k=0; k<3; k++)
  {
    DrawObject(hdv[k], "", Form("cdv%d",k), cList);
    ltx.DrawLatex(0.15, 0.85, Form("Mean %.3f", hdv[k]->GetMean()));
    ltx.DrawLatex(0.15, 0.80, Form("Std dev %.3f", hdv[k]->GetRMS()));
  }

  return;
}

TH1D *
Hist(const char *prefix, int layer, int ladder_or_arm, int stage)
{
  const char *name = Form("h%s%d%d%d", prefix, layer, ladder_or_arm, stage);
  TH1D *h = (TH1D *)gDirectory->GetList()->FindObject(name);
  if (!h)
  {
    Error("Hist() in DrawResults.C", "%s not found. "
          "Available objects in current directory:", name);
    gDirectory->ls();
  }
  assert(h);
  return h;
}

TH1D *
Hist(const char *prefix, int layer, int stage)
{
  const char *name = Form("h%s%d%d", prefix, layer, stage);
  TH1D *h = (TH1D *)gDirectory->Get(name);
  if (!h)
  {
    Error("Hist() in DrawResults.C", "%s not found. "
          "Available objects in current directory:", name);
    gDirectory->ls();
  }
  assert(h);
  return h;
}

void
InitResidHists(const char *treename, int ntrees)
{
  const int nlayers = 4;
  const int narms = 2; // 0=east, 1=west
  const int nladders = 24;

  Info("", "Initializing histograms");
  // Residuals for ladder units
  TH1D *hs[nlayers][nladders][ntrees];
  TH1D *hz[nlayers][nladders][ntrees];
  for (int lyr = 0; lyr < nlayers; lyr++)
    for (int ldr = 0; ldr < nLaddersPerLayer[lyr]; ldr++)
      for (int stage = 0; stage < ntrees; stage++)
      {
        int nbins   = 100;
        double xmax = 0.12;
        if (string(treename).find("cnt") != string::npos)
          xmax = 0.3;
        hs[lyr][ldr][stage] = new TH1D(Form("hs%d%d%d", lyr, ldr, stage),
                                       Form("B%dL%d", lyr, ldr),
                                       nbins, -xmax, xmax);
        hz[lyr][ldr][stage] = new TH1D(Form("hz%d%d%d", lyr, ldr, stage),
                                       Form("B%dL%d", lyr, ldr),
                                       nbins, -xmax, xmax);
        SetupHist(hs[lyr][ldr][stage], stage);
        SetupHist(hz[lyr][ldr][stage], stage);

        // assert(gDirectory->Get(Form("hs%d%d%d", lyr, ldr, stage)));
        // assert(gDirectory->Get(Form("hz%d%d%d", lyr, ldr, stage)));
      }

  // Residuals for half-layer units
  const char *armlabel[narms] = {"E", "W"};
  TH1D *hls[nlayers][narms][ntrees];
  TH1D *hlz[nlayers][narms][ntrees];
  for (int lyr = 0; lyr < nlayers; lyr++)
    for (int arm = 0; arm < narms; arm++)
      for (int stage = 0; stage < ntrees; stage++)
      {
        int nbins   = 200;
        double xmax = 0.06;
        if (string(treename).find("cnt") != string::npos)
          xmax = 0.2;
        hls[lyr][arm][stage] = new TH1D(Form("hls%d%d%d", lyr, arm, stage),
                                        Form("B%d%s", lyr, armlabel[arm]),
                                        nbins, -xmax, xmax);
        hlz[lyr][arm][stage] = new TH1D(Form("hlz%d%d%d", lyr, arm, stage),
                                        Form("B%d%s", lyr, armlabel[arm]),
                                        nbins, -xmax, xmax);
      }

  // Residual means and widths vs ladder
  TH1D *hds[nlayers][ntrees];
  TH1D *hdz[nlayers][ntrees];
  for (int lyr = 0; lyr < 4; lyr++)
  {
    for (int stage = 0; stage < ntrees; stage++)
    {
      hdz[lyr][stage] = new TH1D(Form("hdz%d%d", lyr, stage),
                                 Form("z-resid layer %d", lyr),
                                 nLaddersPerLayer[lyr], -0.5,
                                 nLaddersPerLayer[lyr] - 0.5);
      hds[lyr][stage] = new TH1D(Form("hds%d%d", lyr, stage),
                                 Form("s-resid layer %d", lyr),
                                 nLaddersPerLayer[lyr], -0.5,
                                 nLaddersPerLayer[lyr] - 0.5);
      SetHistProps(hds[lyr][stage], kRed, kNone, kRed, kFullCircle, 1.5);
      SetHistProps(hdz[lyr][stage], kBlue, kNone, kBlue, kFullCircle, 1.5);
      if (stage == 0)
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

    }
  }

  return;
}

void
SetupHist(TH1D *h, int stage)
{
  if (stage == 0)
  {
    h->SetLineColor(kGray + 1);
    h->SetFillColor(kGray);
  }
  else if (TString(h->GetName()).Contains("hs"))
    h->SetLineColor(kGreen + 3);
  else if (TString(h->GetName()).Contains("hz"))
    h->SetLineColor(kAzure + 2);
  else if (TString(h->GetName()).Contains("hls"))
    h->SetLineColor(kGreen + 3);
  else if (TString(h->GetName()).Contains("hlz"))
    h->SetLineColor(kAzure + 2);

  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();
  ay->SetRangeUser(0, 1.3 * h->GetMaximum());
  ax->SetNdivisions(205);
  ax->SetLabelSize(0.06);
  ay->SetNdivisions(205);
  ay->SetLabelSize(0.05);
  ax->SetTitleSize(0);
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
  ltx.DrawLatex(0.62, 0.93, Form("#color[%d]{%.0f (%.0f) #mum}",
                                 kBlack /*h1->GetLineColor()*/,
                                 1e4 * h1->GetMean(),
                                 1e4 * h1->GetRMS()));
  if (h2)
    ltx.DrawLatex(0.62, 0.87, Form("#color[%d]{%.0f (%.0f) #mum}",
                                   h2->GetLineColor(),
                                   1e4 * h2->GetMean(),
                                   1e4 * h2->GetRMS()));
  return;
}

bool
CheckValue(ROOT::TTreeReaderValueBase *value)
{
  if (value->GetSetupStatus() < 0)
  {
    std::cerr << "Error " << value->GetSetupStatus()
              << " setting up reader for " << value->GetBranchName() << '\n';
    return false;
  }
  return true;
}

// void
// MarkDeadLadders(TGraphErrors **gdead)
// {
//   for (int lyr = 0; lyr < nlayers; lyr++)
//   {
//     gdead[lyr] = new TGraphErrors();
//     SetGraphProps(gdead[lyr], kNone, kOrange - 8, kNone, kDot);
//     gdead[lyr]->SetFillColorAlpha(kOrange - 8, 0.5);
//     gdead[lyr]->SetLineWidth(0);
//     for (int ldr = 0; ldr < nLaddersPerLayer[lyr]; ldr++)
//     {
//       // BadLadders.h lists dead ladders.
//       double dx = Dead(lyr, ldr) ? 0.5 : 0.0;
//       double dy = Dead(lyr, ldr) ? 0.089 : 0.0;
//       gdead[lyr]->SetPoint(ldr, ldr, 0.0);
//       gdead[lyr]->SetPointError(ldr, dx, dy);
//     }
//   }
//   return;
// }

void
DrawHalfLayerResidPlots(int ntrees, TObjArray *cList)
{
  // Draw halflayer s and z residuals on sub pads
  TCanvas *cls = new TCanvas(Form("cls%d", 0), Form("%d", 0), 1200, 800);
  TCanvas *clz = new TCanvas(Form("clz%d", 0), Form("%d", 0), 1200, 800);
  cls->Divide(4, 2, 0.001, 0.001);
  clz->Divide(4, 2, 0.001, 0.001);
  for (int lyr = 0; lyr < 4; lyr++)
  {
    Info("", "Drawing halflayer residuals in layer %d", lyr);
    for (int arm = 0; arm < 2; arm++)
    {
      for (int stage = 0; stage < ntrees; stage++)
      {
        cls->cd(arm * 4 + lyr + 1);
        TH1D *h = Hist("ls",lyr,arm,stage);
        h->Draw(stage == 0 ? "" : "same");

        clz->cd(arm * 4 + lyr + 1);
        h = Hist("lz",lyr,arm,stage);
        h->Draw(stage == 0 ? "" : "same");
      }
      if (ntrees > 1)
      {
        cls->cd(arm * 4 + lyr + 1);
        TH1D *h0 = Hist("ls",lyr,arm,0);
        TH1D *h1 = Hist("ls",lyr,arm,1);
        SetYMax(h0, h1);
        ModifyPad(gPad, h0, h1, "ds");

        clz->cd(arm * 4 + lyr + 1);
        h0 = Hist("lz",lyr,arm,0);
        h1 = Hist("lz",lyr,arm,1);
        SetYMax(h0, h1);
        ModifyPad(gPad, h0, h1, "dz");
      }
    }
    cList->Add(cls);
    cList->Add(clz);
  }
  return;
}

void
DrawLadderResidPlots(int ntrees, TObjArray *cList)
{
  // Draw individual s and z ladder residual distributions on sub-pads
  for (int lyr = 0; lyr < 4; lyr++)
  {
    Info("", "Drawing ladder residuals in layer %d", lyr);
    int nl = nLaddersPerLayer[lyr];
    int nx = lyr ? nl / 4 : nl / 2; // Number of pads along x
    int ny = lyr ? 4 : 2;       // Number of pads along y
    int ph = 200;               // pad height in px
    int pw = 200;               // pad width in px
    TCanvas *cs = new TCanvas(Form("cs%d", lyr), Form("ds layer %d", lyr),
                              pw * nx, ph * ny);
    TCanvas *cz = new TCanvas(Form("cz%d", lyr), Form("dz layer %d", lyr),
                              pw * nx, ph * ny);
    cs->Divide(nx, ny, 0.001, 0.001);
    cz->Divide(nx, ny, 0.001, 0.001);

    for (int ldr = 0; ldr < nl; ldr++)
    {
      for (int stage = 0; stage < ntrees; stage++)
      {
        cs->cd(ldr + 1);
        TH1D *h = Hist("s",lyr,ldr,stage);
        h->Draw(stage == 0 ? "" : "same");

        cz->cd(ldr + 1);
        h = Hist("z",lyr,ldr,stage);
        h->Draw(stage == 0 ? "" : "same");
      }
      if (ntrees > 1)
      {
        cs->cd(ldr + 1);
        TH1D *h0 = Hist("s",lyr,ldr,0);
        TH1D *h1 = Hist("s",lyr,ldr,1);
        SetYMax(h0, h1);
        ModifyPad(gPad, h0, h1, "ds");

        cz->cd(ldr + 1);
        h0 = Hist("z",lyr,ldr,0);
        h1 = Hist("z",lyr,ldr,1);
        SetYMax(h0, h1);
        ModifyPad(gPad, h0, h1, "dz");
      }
      //   if (Dead(lyr, ldr)) // Put a big "X" over dead ladders
      //   {
      //     TLatex l;
      //     l.SetTextColor(kGray + 2);
      //     l.SetTextSize(0.4);
      //     l.SetNDC();
      //     cs->cd(ldr + 1);
      //     l.DrawLatex(0.5, 0.5, "#times");
      //     cz->cd(ldr + 1);
      //     l.DrawLatex(0.5, 0.5, "#times");
      //   }
      // }
    }
    cList->Add((TCanvas *)cs);
    cList->Add((TCanvas *)cz);
  }
  return;
}

void
DrawSummaryPlots(int ntrees, TObjArray *cList)
{
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);

  // Draw summary plots of residual mean vs ladder index for each layer
  for (int lyr = 0; lyr < 4; lyr++)
  {
    Info("", "Drawing residual summary plots in layer %d", lyr);
    // Mean of means, mean of std devs
    double smm = 0, smrms = 0;
    double zmm = 0, zmrms = 0;
    for (int stage = 0; stage < ntrees; stage++)
    {
      for (int ldr = 0; ldr < nLaddersPerLayer[lyr]; ldr++)
      {
        // Get mean and width from individual ladder residual dists
        double sm   = Hist("s",lyr,ldr,stage)->GetMean();
        double zm   = Hist("z",lyr,ldr,stage)->GetMean();
        double srms = Hist("s",lyr,ldr,stage)->GetRMS();
        double zrms = Hist("z",lyr,ldr,stage)->GetRMS();

        if (stage == 1)
        {
          smm   += sm / nLaddersPerLayer[lyr];
          zmm   += zm / nLaddersPerLayer[lyr];
          smrms += srms / nLaddersPerLayer[lyr];
          zmrms += zrms / nLaddersPerLayer[lyr];
        }

        // Assign ladder residual means to bins of summary histograms
        Hist("ds",lyr,stage)->SetBinContent(ldr + 1, sm);
        Hist("dz",lyr,stage)->SetBinContent(ldr + 1, zm);
        Hist("ds",lyr,stage)->SetBinError(ldr + 1, srms);
        Hist("dz",lyr,stage)->SetBinError(ldr + 1, zrms);
      }
    }

    DrawObject(Hist("ds",lyr,0), "e0p", Form("ds_lyr%d", lyr), cList);
    // gdead[lyr]->Draw("e5p,same");
    if (ntrees > 1)
      Hist("ds",lyr,1)->Draw("e0p,same");
    gPad->RedrawAxis();
    ltx.SetTextSize(0.07);
    ltx.DrawLatex(0.23, 0.92, Hist("ds",lyr,0)->GetTitle());
    ltx.SetTextSize(0.05);
    ltx.DrawLatex(0.6, 0.85, Form("#LTmean#GT    %.3f", smm));
    ltx.DrawLatex(0.6, 0.80, Form("#LTStd dev#GT %.3f", smrms));

    DrawObject(Hist("dz",lyr,0), "e0p", Form("dz_lyr%d", lyr), cList);
    // gdead[lyr]->Draw("e5p,same");
    if (ntrees > 1)
      Hist("dz",lyr,1)->Draw("e0p,same");
    gPad->RedrawAxis();
    ltx.SetTextSize(0.07);
    ltx.DrawLatex(0.23, 0.92, Hist("dz",lyr,0)->GetTitle());
    ltx.SetTextSize(0.05);
    ltx.DrawLatex(0.6, 0.85, Form("#LTmean#GT    %.3f", zmm));
    ltx.DrawLatex(0.6, 0.80, Form("#LTStd dev#GT %.3f", zmrms));
  }

  return;
}
