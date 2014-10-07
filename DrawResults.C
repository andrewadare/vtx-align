#include "VtxAlignBase.h"
#include "BadLadders.h"
#include "VtxVis.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

bool openPDF = true; // Mac only, but simple to adapt to other viewers with CLI

int nLaddersPerLayer[4] = {10, 20, 16, 24}; // ladders/layer

void DrawHalfLayerResidPlots(int ntrees, TObjArray *cList = 0);
void DrawLadderResidPlots(int ntrees, TObjArray *cList = 0);
void DrawSummaryPlots(int ntrees, TObjArray *cList = 0);
void InitResidHists(const char *treename, int ntrees);
void SetupHist(TH1D *h, int stage);
void SetAxisProps(TH1 *h, int xdivs = 208, int ydivs = 208,
                  double labelSize = 0.04, double titleSize = 0.04,
                  double xTitleOffset = 1.1, double yTitleOffset = 1.5,
                  bool centerX = true, bool centerY = true);

void FillHists(TFile *f, const char *treename, int stage, int prod, int subiter,
               TObjArray *cList = 0);
void ModifyPad(TVirtualPad *pad);
void PrintMeanToPad(TH1D *h1, TH1D *h2, TString coord);
bool CheckValue(ROOT::TTreeReaderValueBase *value);
TH1D *Hist(const char *prefix, int layer, int ladder_or_arm, int stage);
TH1D *Hist(const char *prefix, int layer, int stage);
TGraphErrors *DeadLadderGraph(int lyr);

// To plot variables from only one TTree, either:
// 1. set prod2 and subit2 to any int < 0, or
// 2. set prod2 = prod1, subit2 = subit1.
void DrawResults(int run = 411768,
                 int prod1 = 0,
                 int subit1 = 0,
                 int prod2 = 0,
                 int subit2 = 1,
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
    inFiles.push_back(new TFile(rootFileIn[stage].c_str(), "read"));
    assert(inFiles[stage]);
  }

  InitResidHists(treename, ntrees);
  FillHists(inFiles[0], treename, 0, prod1, subit1, cList);
  if (ntrees == 2)
    FillHists(inFiles[1], treename, 1, prod2, subit2, cList);

  DrawHalfLayerResidPlots(ntrees, cList);
  DrawLadderResidPlots(ntrees, cList);
  DrawSummaryPlots(ntrees, cList);

  if (ntrees == 1)
  {
    PrintPDFs(cList, Form("pdfs/run%d-pro%dsub%d-%s",
                          run, prod1, subit1, treename), "");
    const char *pdfName = Form("pdfs/run%d-pro%dsub%d-%s",
                               run, prod1, subit1, treename);
    PrintPDF(cList, pdfName, "");
    if (openPDF)
      gSystem->Exec(Form("open %s.pdf", pdfName));
  }
  if (ntrees == 2)
  {
    PrintPDFs(cList, Form("pdfs/run%d-pro%dsub%d-vs-pro%dsub%d-%s",
                          run, prod1, subit1, prod2, subit2, treename), "");
    const char *pdfName = Form("pdfs/run%d-pro%dsub%d-vs-pro%dsub%d-%s",
                               run, prod1, subit1, prod2, subit2, treename);
    PrintPDF(cList, pdfName, "");
    if (openPDF)
      gSystem->Exec(Form("open %s.pdf", pdfName));
  }
  if (ntrees == 2)
  {
    gStyle->SetOptTitle(1);
    const char *parfile1 = Form("geom/%d-%d-%d.par", run, prod1, subit1);
    const char *parfile2 = Form("geom/%d-%d-%d.par", run, prod2, subit2);
    gROOT->Macro(Form("DiffGeometry.C(\"%s\", \"%s\")", parfile1, parfile2));

    if (openPDF)
      gSystem->Exec(Form("open pdfs/%d-%d-%d-vs-%d-%d-%d.pdf",
                         run, prod1, subit1, run, prod2, subit2));
  }
  return;
}

void
FillHists(TFile *f, const char *treename, int stage, int prod, int subiter,
          TObjArray *cList)
{
  TH1D *xydcae = new TH1D(Form("xydcae_%d_%d",prod,subiter),
                          Form(";east arm x-y DCA [cm];tracks"), 200, -0.1, 0.1);
  TH1D *xydcaw = new TH1D(Form("xydcaw_%d_%d",prod,subiter),
                          Form(";west arm x-y DCA [cm];tracks"), 200, -0.1, 0.1);
  TH1D *zdcae = new TH1D(Form("zdcae_%d_%d",prod,subiter),
                         ";east z DCA [cm];tracks", 200, -0.25, 0.25);
  TH1D *zdcaw = new TH1D(Form("zdcaw_%d_%d",prod,subiter),
                         ";west z DCA [cm];tracks", 200, -0.25, 0.25);
  TH1D *hvze = new TH1D(Form("hvze_%d_%d",prod,subiter),
                        Form(";east arm z vertex [cm];tracks"), 200, -15, 15);
  TH1D *hvzw = new TH1D(Form("hvzw_%d_%d",prod,subiter),
                        Form(";west arm z vertex [cm];tracks"), 200, -15, 15);
  TH2D *xydcaphi  = new TH2D(Form("xydcaphi_%d_%d",prod,subiter),
                             ";#phi [rad];x-y DCA [cm]",
                             100, -TMath::PiOver2(), 3*TMath::PiOver2(),
                             100, -0.1, +0.1);
  TH2D *ezdcatheta  = new TH2D(Form("ezdcatheta_%d_%d",prod,subiter),
                               ";#theta [rad];z DCA [cm]",
                               100, 0.22*TMath::Pi(), 0.78*TMath::Pi(),
                               100, -0.25, +0.25);
  TH2D *wzdcatheta  = new TH2D(Form("wzdcatheta_%d_%d",prod,subiter),
                               ";#theta [rad];z DCA [cm]",
                               100, 0.22*TMath::Pi(), 0.78*TMath::Pi(),
                               100, -0.25, +0.25);
  SetAxisProps(xydcae);
  SetAxisProps(xydcaw);
  SetAxisProps(zdcae, 208, 208, 0.04, 0.04, 1.1, 1.7);
  SetAxisProps(zdcaw, 208, 208, 0.04, 0.04, 1.1, 1.7);
  SetAxisProps(hvze);
  SetAxisProps(hvzw);
  SetAxisProps(xydcaphi);
  SetAxisProps(ezdcatheta);
  SetAxisProps(wzdcatheta);

  // x-y vertex distributions
  double x0 = -0.5, y0 = -0.5, x1 = +0.5, y1 = +0.5;
  TH2D *hve = new TH2D(Form("hve_%d_%d",prod,subiter),
                       Form("East vertex - prod %d step %d;"
                            "x [cm];y [cm]", prod, subiter),
                       500, x0, x1, 500, y0, y1);
  SetAxisProps(hve, 205, 205);

  TH2D *hvw = new TH2D(Form("hvw_%d_%d",prod,subiter),
                       Form("West vertex - prod %d step %d;"
                            "x [cm];y [cm]", prod, subiter),
                       500, x0, x1, 500, y0, y1);
  SetAxisProps(hvw, 205, 205);

  // West - East offset histograms
  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k=0; k<3; k++)
  {
    double lim = 0.41;
    if (k == 2)
      lim *= 2;

    hdv[k] = new TH1D(Form("hdv_%d__%d_%d",prod,subiter, k),
                      Form("West - East vertex difference #Delta%s "
                           "- prod %d step %d; (W-E) #Delta%s [cm];tracks",
                           xyzstr[k], prod, subiter, xyzstr[k]),
                      200, -lim, lim);
    SetAxisProps(hdv[k]);
  }

  TTreeReader r(treename, f);

  TTreeReaderValue<int>   event(r, "event");
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
  int prevev = -1;
  int ne = 0, nw = 0;
  TVectorD ve(3);
  TVectorD vw(3);
  while (r.Next())
  {
    // TEMP //////////////////////////////////////////////////
    if (fabs(*xydca) < 1e-12)
      continue;
    // TEMP //////////////////////////////////////////////////

    // Check for a new event
    if (*event != prevev)
    {

      // Fill West - East vertex difference histograms
      if (ne > 0 && nw > 0)
      {
        for (int k=0; k<3; k++)
          hdv[k]->Fill(vw(k) - ve(k));
      }

      // Reset for next event
      prevev = *event;
      ve *= 0;
      vw *= 0;
      ne = 0;
      nw = 0;
    }

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
    if (arm == 0)
    {
      ++ne;
      if (ne == 1)
      {
        ve(0) = vertex[0];
        ve(1) = vertex[1];
        ve(2) = vertex[2];
      }
      xydcae->Fill(*xydca);
      zdcae->Fill(*zdca);
      hve->Fill(vertex[0], vertex[1]);
      hvze->Fill(vertex[2]);
      ezdcatheta->Fill(*theta, *zdca);
    }
    if (arm == 1)
    {
      ++nw;
      if (nw == 1)
      {
        vw(0) = vertex[0];
        vw(1) = vertex[1];
        vw(2) = vertex[2];
      }
      xydcaw->Fill(*xydca);
      zdcaw->Fill(*zdca);
      hvw->Fill(vertex[0], vertex[1]);
      hvzw->Fill(vertex[2]);
      wzdcatheta->Fill(*theta, *zdca);
    }

    xydcaphi->Fill(phiwrap, *xydca);
  }


  // Plot--------------------------------------------------------
  Info("", "Plotting histograms from %s", f->GetName());
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);
  TCanvas *c = 0;

  // East and West xy vertex distributions
  hve->GetXaxis()->SetRangeUser(hve->GetMean(1)-0.1,hve->GetMean(1)+0.1);
  hve->GetYaxis()->SetRangeUser(hve->GetMean(2)-0.1,hve->GetMean(2)+0.1);
  hvw->GetXaxis()->SetRangeUser(hvw->GetMean(1)-0.1,hvw->GetMean(1)+0.1);
  hvw->GetYaxis()->SetRangeUser(hvw->GetMean(2)-0.1,hvw->GetMean(2)+0.1);
  c = new TCanvas(Form("xy_vertex_%d",stage),
                  Form("xy_vertex_ew_%d",stage), 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  hve->Draw("col");
  ltx.DrawLatex(0.25, 0.92, "East");
  ltx.DrawLatex(0.45, 0.92, Form("mean (%.0f, %.0f) #mum",
                                 1e4*hve->GetMean(1), 1e4*hve->GetMean(2)));
  ltx.DrawLatex(0.45, 0.87, Form("std dev (%.0f, %.0f) #mum",
                                 1e4*hve->GetRMS(1), 1e4*hve->GetRMS(2)));
  c->cd(2);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  hvw->Draw("col");
  ltx.DrawLatex(0.25, 0.92, "West");
  ltx.DrawLatex(0.45, 0.92, Form("mean (%.0f, %.0f) #mum",
                                 1e4*hvw->GetMean(1), 1e4*hvw->GetMean(2)));
  ltx.DrawLatex(0.45, 0.87, Form("std dev (%.0f, %.0f) #mum",
                                 1e4*hvw->GetRMS(1), 1e4*hvw->GetRMS(2)));
  cList->Add(c);

  // East and West z vertex distributions
  c = new TCanvas(Form("z_vertex_%d",stage),
                  Form("z_vertex_ew_%d",stage), 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  hvze->GetYaxis()->SetRangeUser(0, 1.2*hvze->GetMaximum());
  hvze->Draw("");
  ltx.DrawLatex(0.25, 0.92, "East");
  ltx.DrawLatex(0.45, 0.92, Form("mean %.3f cm", hvze->GetMean(1)));
  ltx.DrawLatex(0.45, 0.87, Form("std dev %.3f cm", hvze->GetRMS(1)));
  c->cd(2);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  hvzw->GetYaxis()->SetRangeUser(0, 1.2*hvze->GetMaximum());
  hvzw->Draw("");
  ltx.DrawLatex(0.25, 0.92, "West");
  ltx.DrawLatex(0.45, 0.92, Form("mean %.3f cm", hvzw->GetMean(1)));
  ltx.DrawLatex(0.45, 0.87, Form("std dev %.3f cm", hvzw->GetRMS(1)));
  cList->Add(c);

  // East-west offsets in x, y, and z
  c = new TCanvas(Form("ew_offsets_%d",stage),
                  Form("ew_offsets_%d",stage), 1000, 500);
  c->Divide(3,1, 0.001, 0.001);
  for (int k=0; k<3; k++)
  {
    c->cd(k+1);
    gPad->SetMargin(0.15, 0.02, 0.12, 0.02); // L, R, B, T
    hdv[k]->GetYaxis()->SetRangeUser(0, 1.2*hdv[k]->GetMaximum());
    hdv[k]->Draw();
    ltx.DrawLatex(0.25, 0.95, Form("Mean %.3f", hdv[k]->GetMean()));
    ltx.DrawLatex(0.25, 0.90, Form("Std dev %.3f", hdv[k]->GetRMS()));
  }
  cList->Add(c);

  // XY DCA - East and West
  c = new TCanvas(Form("cxydca_%d_%d",prod,subiter),
                  Form("xy_dca_east_west__%d_%d",prod,subiter), 1000, 500);
  SetYMax(xydcae, xydcaw);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.01, 0.12, 0.01); // L, R, B, T
  xydcae->Draw();
  ltx.DrawLatex(0.16, 0.93, Form("Mean %.0f #mum    Std dev %.0f #mum",
                                 1e4*xydcae->GetMean(), 1e4*xydcae->GetRMS()));
  c->cd(2);
  gPad->SetMargin(0.12, 0.01, 0.12, 0.01); // L, R, B, T
  xydcaw->Draw();
  ltx.DrawLatex(0.16, 0.93, Form("Mean %.0f #mum    Std dev %.0f #mum",
                                 1e4*xydcaw->GetMean(), 1e4*xydcaw->GetRMS()));
  cList->Add(c);

  // XY DCA vs phi
  DrawObject(xydcaphi, "col",  Form("xydca_vs_phi_%d",stage), cList);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  TProfile *dca2dprof = xydcaphi->ProfileX(Form("dca2dprof%d_%d",prod,subiter), 
                                           1, -1, "d,same");
  dca2dprof->SetMarkerStyle(kFullCircle);

  // Z DCA - East and West
  c = new TCanvas(Form("czdca_%d_%d",prod,subiter),
                  Form("z_dca_east_west_%d",stage), 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.01, 0.12, 0.01); // L, R, B, T
  zdcae->GetYaxis()->SetRangeUser(0, 1.2*zdcae->GetMaximum());
  zdcae->Draw();
  ltx.DrawLatex(0.15, 0.93, Form("Mean %.0f #mum    Std dev %.0f #mum",
                                 1e4*zdcae->GetMean(), 1e4*zdcae->GetRMS()));
  c->cd(2);
  gPad->SetMargin(0.12, 0.01, 0.12, 0.01); // L, R, B, T
  zdcaw->GetYaxis()->SetRangeUser(0, 1.2*zdcaw->GetMaximum());
  zdcaw->Draw();
  ltx.DrawLatex(0.15, 0.93, Form("Mean %.0f #mum    Std dev %.0f #mum",
                                 1e4*zdcaw->GetMean(), 1e4*zdcaw->GetRMS()));
  cList->Add(c);


  // Z DCA vs theta - East and West
  c = new TCanvas(Form("z_dca_vs_theta_%d",stage),
                  Form("z_dca_vs_theta_ew_%d",stage), 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  ezdcatheta->Draw("col");
  TProfile *zeprof = ezdcatheta->ProfileX(Form("zeprof%d_%d",prod,subiter), 1, -1, "d,same");
  zeprof->SetMarkerStyle(kFullCircle);
  ltx.DrawLatex(0.25, 0.92, "East");
  c->cd(2);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  wzdcatheta->Draw("col");
  TProfile *zwprof = wzdcatheta->ProfileX(Form("zwprof%d_%d",prod,subiter), 1, -1, "d,same");
  zwprof->SetMarkerStyle(kFullCircle);
  ltx.DrawLatex(0.25, 0.92, "West");
  cList->Add(c);

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

  h->GetYaxis()->SetRangeUser(0, 1.3 * h->GetMaximum());
  SetAxisProps(h, 205, 205, 0.06, 0.0);
  return;
}

void
SetAxisProps(TH1 *h, int xdivs, int ydivs,
             double labelSize, double titleSize,
             double xTitleOffset, double yTitleOffset,
             bool centerX, bool centerY)
{
  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();
  ax->SetNdivisions(xdivs);
  ay->SetNdivisions(ydivs);

  ax->SetLabelSize(labelSize);
  ay->SetLabelSize(labelSize);

  ax->SetTitleSize(titleSize);
  ay->SetTitleSize(titleSize);

  ax->SetTitleOffset(xTitleOffset);
  ay->SetTitleOffset(yTitleOffset);

  ax->CenterTitle(centerX);
  ay->CenterTitle(centerY);

  return;
}

void
ModifyPad(TVirtualPad *pad)
{
  pad->SetBottomMargin(0.15);
  pad->SetLeftMargin(0.2);
  pad->SetRightMargin(0.01);
  pad->SetTopMargin(0.01);

  return;
}

void
PrintMeanToPad(TH1D *h1, TH1D *h2, TString coord)
{
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

TGraphErrors *
DeadLadderGraph(int lyr)
{
  static int ncalls = 0;
  TGraphErrors *g = new TGraphErrors(nLaddersPerLayer[lyr]);
  g->SetName(Form("gdead%d", ncalls));

  SetGraphProps(g, kNone, kOrange - 8, kNone, kDot);
  g->SetFillColorAlpha(kOrange - 8, 0.5);
  g->SetLineWidth(0);
  for (int ldr = 0; ldr < nLaddersPerLayer[lyr]; ldr++)
  {
    // BadLadders.h lists dead ladders.
    double dx = Dead(lyr, ldr) ? 0.5 : 0.0;
    double dy = Dead(lyr, ldr) ? 0.089 : 0.0;
    g->SetPoint(ldr, ldr, 0.0);
    g->SetPointError(ldr, dx, dy);
  }

  return g;
}

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
        SetupHist(h, stage);
        ModifyPad(gPad);
        if (ntrees == 1)
          PrintMeanToPad(h, 0, "ds");

        clz->cd(arm * 4 + lyr + 1);
        h = Hist("lz",lyr,arm,stage);
        h->Draw(stage == 0 ? "" : "same");
        ModifyPad(gPad);
        SetupHist(h, stage);
        if (ntrees == 1)
          PrintMeanToPad(h, 0, "dz");
      }
      if (ntrees > 1)
      {
        cls->cd(arm * 4 + lyr + 1);
        TH1D *h0 = Hist("ls",lyr,arm,0);
        TH1D *h1 = Hist("ls",lyr,arm,1);
        SetYMax(h0, h1);
        PrintMeanToPad(h0, h1, "ds");

        clz->cd(arm * 4 + lyr + 1);
        h0 = Hist("lz",lyr,arm,0);
        h1 = Hist("lz",lyr,arm,1);
        SetYMax(h0, h1);
        PrintMeanToPad(h0, h1, "dz");
      }
    }
  }
  cList->Add(cls);
  cList->Add(clz);
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
        ModifyPad(gPad);
        SetupHist(h, stage);
        if (ntrees == 1)
          PrintMeanToPad(h, 0, "ds");

        cz->cd(ldr + 1);
        h = Hist("z",lyr,ldr,stage);
        h->Draw(stage == 0 ? "" : "same");
        ModifyPad(gPad);
        SetupHist(h, stage);
        if (ntrees == 1)
          PrintMeanToPad(h, 0, "dz");
      }

      if (ntrees > 1)
      {
        cs->cd(ldr + 1);
        TH1D *h0 = Hist("s",lyr,ldr,0);
        TH1D *h1 = Hist("s",lyr,ldr,1);
        SetYMax(h0, h1);
        PrintMeanToPad(h0, h1, "ds");

        cz->cd(ldr + 1);
        h0 = Hist("z",lyr,ldr,0);
        h1 = Hist("z",lyr,ldr,1);
        SetYMax(h0, h1);
        PrintMeanToPad(h0, h1, "dz");
      }
      if (Dead(lyr, ldr)) // Put a big "X" over dead ladders
      {
        TLatex l;
        l.SetTextColor(kGray + 2);
        l.SetTextSize(0.4);
        l.SetNDC();
        cs->cd(ldr + 1);
        l.DrawLatex(0.5, 0.5, "#times");
        cz->cd(ldr + 1);
        l.DrawLatex(0.5, 0.5, "#times");
      }
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
  ltx.SetTextSize(0.05);
  TLatex ltx2;
  ltx2.SetNDC();
  ltx2.SetTextFont(42);
  ltx2.SetTextSize(0.07);

  // Draw summary plots of residual mean vs ladder index for each layer
  for (int lyr = 0; lyr < 4; lyr++)
  {
    Info("", "Drawing residual summary plots in layer %d", lyr);

    // Precompute mean of means, mean of std devs
    double smm[2] = {0}, smrms[2] = {0};
    double zmm[2] = {0}, zmrms[2] = {0};
    for (int stage = 0; stage < ntrees; stage++)
    {
      for (int ldr = 0; ldr < nLaddersPerLayer[lyr]; ldr++)
      {
        // Get mean and width from individual ladder residual dists
        double sm   = Hist("s",lyr,ldr,stage)->GetMean();
        double zm   = Hist("z",lyr,ldr,stage)->GetMean();
        double srms = Hist("s",lyr,ldr,stage)->GetRMS();
        double zrms = Hist("z",lyr,ldr,stage)->GetRMS();

        smm[stage]   += sm / nLaddersPerLayer[lyr];
        zmm[stage]   += zm / nLaddersPerLayer[lyr];
        smrms[stage] += srms / nLaddersPerLayer[lyr];
        zmrms[stage] += zrms / nLaddersPerLayer[lyr];

        // Assign ladder residual means to bins of summary histograms
        Hist("ds",lyr,stage)->SetBinContent(ldr + 1, sm);
        Hist("dz",lyr,stage)->SetBinContent(ldr + 1, zm);
        Hist("ds",lyr,stage)->SetBinError(ldr + 1, srms);
        Hist("dz",lyr,stage)->SetBinError(ldr + 1, zrms);
      }
    }

    for (int stage = 0; stage < ntrees; stage++)
    {
      // TGraphErrors *gdead = DeadLadderGraph(lyr);
      if (stage == 0)
        DrawObject(Hist("ds",lyr,stage), "e0p", Form("ds_lyr%d", lyr), cList);
      else
        Hist("ds",lyr,stage)->Draw("e0p,same");

      // gdead->Draw("e5p,same");
      ltx.SetTextColor(Hist("ds",lyr,stage)->GetMarkerColor()+1);
      ltx.DrawLatex(stage?0.6:0.2, 0.85,
                    Form("#LTmean#GT    %.0f #mum", 1e4*smm[stage]));
      ltx.DrawLatex(stage?0.6:0.2, 0.80,
                    Form("#LTStd dev#GT %.0f #mum", 1e4*smrms[stage]));
      ltx2.DrawLatex(0.3, 0.92, Form("Layer %d #Deltas vs ladder", lyr));
      gPad->RedrawAxis();
      gPad->SetRightMargin(0.02);
    }

    for (int stage = 0; stage < ntrees; stage++)
    {
      if (stage == 0)
        DrawObject(Hist("dz",lyr,stage), "e0p", Form("dz_lyr%d", lyr), cList);
      else
        Hist("dz",lyr,stage)->Draw("e0p,same");

      // gdead->Draw("e5p,same");
      ltx.SetTextColor(Hist("dz",lyr,stage)->GetMarkerColor()+1);
      ltx.DrawLatex(stage?0.6:0.2, 0.85,
                    Form("#LTmean#GT    %.0f #mum", 1e4*zmm[stage]));
      ltx.DrawLatex(stage?0.6:0.2, 0.80,
                    Form("#LTStd dev#GT %.0f #mum", 1e4*zmrms[stage]));
      ltx2.DrawLatex(0.3, 0.92, Form("Layer %d #Deltaz vs ladder", lyr));
      gPad->RedrawAxis();
      gPad->SetRightMargin(0.02);
    }

  }

  return;
}
