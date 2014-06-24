#include "UtilFns.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

const int nlayers = 4;
const int nladders = 24;
const int ntrees = 2;

bool Dead(int layer, int ladder);
void FillHists(TH1D *hs[nlayers][nladders][ntrees],
               TH1D *hz[nlayers][nladders][ntrees],
               TFile *f, const char *treename, int stage);
void SetupHist(TH1D *h, int stage);

void DrawResults()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int nLadders[4] = {10,20,16,24}; // ladders/layer

  TFile *inFile = new TFile("rootfiles/411768_cluster.1.root", "read");
  TObjArray *cList = new TObjArray();

  TNtuple *ht[ntrees] = {0};
  ht[0] = (TNtuple *) inFile->Get("ht1");
  ht[1] = (TNtuple *) inFile->Get("ht2");

  TLatex ltx;
  ltx.SetNDC();
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

  FillHists(hs, hz, inFile, "ht1", 0);
  FillHists(hs, hz, inFile, "ht2", 1);

  for (int lyr=0; lyr<nlayers; lyr++)
    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      for (int stage=0; stage<ntrees; stage++)
      {
        SetupHist(hs[lyr][ldr][stage], stage);
        SetupHist(hz[lyr][ldr][stage], stage);
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
    TCanvas *cz = new TCanvas(Form("cz%d", lyr), Form("dz layer %d", lyr),
                              pw*nx, ph*ny);
    cs->Divide(nx, ny, 0.001, 0.001);
    cz->Divide(nx, ny, 0.001, 0.001);

    for (int ldr=0; ldr<nl; ldr++)
    {
      cs->cd(ldr+1);
      for (int stage=0; stage<ntrees; stage++)
        hs[lyr][ldr][stage]->Draw(stage==0?"":"same");

      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.01);
      gPad->SetTopMargin(0.01);
      ltx.SetTextSize(0.07);
      ltx.DrawLatex(0.23, 0.92, hs[lyr][ldr][1]->GetTitle());
      ltx.DrawLatex(0.23, 0.85, "ds [cm]");
      ltx.SetTextSize(0.06);
      ltx.DrawLatex(0.62,0.93, Form("%.0f (%.0f) #mum",
                                    1e4*hs[lyr][ldr][0]->GetMean(),
                                    1e4*hs[lyr][ldr][0]->GetRMS()));
      ltx.DrawLatex(0.62,0.87, Form("#color[419]{%.0f (%.0f) #mum}",
                                    1e4*hs[lyr][ldr][1]->GetMean(),
                                    1e4*hs[lyr][ldr][1]->GetRMS()));
      cz->cd(ldr+1);
      for (int stage=0; stage<ntrees; stage++)
        hz[lyr][ldr][stage]->Draw(stage==0?"":"same");

      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.01);
      gPad->SetTopMargin(0.01);
      ltx.SetTextSize(0.07);
      ltx.DrawLatex(0.23, 0.92, hz[lyr][ldr][1]->GetTitle());
      ltx.DrawLatex(0.23, 0.85, "dz [cm]");
      ltx.SetTextSize(0.06);
      ltx.DrawLatex(0.62,0.93, Form("%.0f (%.0f) #mum",
                                    1e4*hz[lyr][ldr][0]->GetMean(),
                                    1e4*hz[lyr][ldr][0]->GetRMS()));
      ltx.DrawLatex(0.62,0.87, Form("#color[419]{%.0f (%.0f) #mum}",
                                    1e4*hz[lyr][ldr][1]->GetMean(),
                                    1e4*hz[lyr][ldr][1]->GetRMS()));
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
      hdz[lyr][stage] = new TH1D(Form("hdz%d_%d",lyr, stage), Form("z-resid layer %d", lyr), nLadders[lyr], 0, float(nLadders[lyr]));
      hds[lyr][stage] = new TH1D(Form("hds%d_%d",lyr, stage), Form("s-resid layer %d", lyr), nLadders[lyr], 0, float(nLadders[lyr]));
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

      hds[lyr][stage]->GetYaxis()->SetRangeUser(-0.09, 0.09);
      hdz[lyr][stage]->GetYaxis()->SetRangeUser(-0.09, 0.09);

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
    hds[lyr][1]->Draw("e0p,same");
    ltx.DrawLatex(0.23, 0.92, hds[lyr][0]->GetTitle());
    DrawObject(hdz[lyr][0], "e0p", Form("dz_lyr%d",lyr), cList);
    hdz[lyr][1]->Draw("e0p,same");
    ltx.DrawLatex(0.23, 0.92, hdz[lyr][0]->GetTitle());
  }
  ltx.SetTextSize(0.06);
  ltx.SetTextFont(42);

  const int nxyp = 2;
  const char *xyplots[nxyp] = {"vtx_xy","millepede_dp"};
  for (int i=0; i<nxyp; i++)
  {
    TCanvas *c = (TCanvas *) inFile->Get(xyplots[i]);
    assert(c);
    c->Draw();
    ltx.DrawLatex(0.15, 0.92, gPad->GetTitle());
    cList->Add(c);
  }

  PrintPDFs(cList, "pdfs", "");
  PrintPDF(cList, "pdfs/self-align");
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

bool
Dead(int layer, int ladder)
{
  if (layer==1 &&
      ladder==11) return true;
  if (layer==3 && ladder==10) return true;
  if (layer==3 && ladder==16) return true;
  if (layer==3 && ladder==23) return true;
  return false;
}

void
FillHists(TH1D *hs[nlayers][nladders][ntrees],
          TH1D *hz[nlayers][nladders][ntrees],
          TFile *f, const char *treename, int stage)
{
  TTreeReader r(treename, f);
  TTreeReaderValue<float> ds(r, "res_s");
  TTreeReaderValue<float> dz(r, "res_z");
  TTreeReaderValue<float> layer(r, "layer");
  TTreeReaderValue<float> ladder(r, "ladder");

  while (r.Next())
  {
    int lyr = *layer;
    int ldr = *ladder;
    float s = *ds;
    float z = *dz;

    // Printf("%s %d %d %d %f %f %x", 
    //        treename, lyr, ldr, stage, s, z, hs[lyr][ldr][stage]);
    if (s>-0.2 && s<0.2 && z>-0.2 && z<0.2)
    {
      hs[lyr][ldr][stage]->Fill(s);
      hz[lyr][ldr][stage]->Fill(z);
    }
  }

  return;
}
