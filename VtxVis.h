#ifndef __VTXVIS_H__
#define __VTXVIS_H__

#include "SvxTGeo.h"
#include "BadLadders.h"

#include <TH1.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>

typedef vector<double> vecd;

void GetLadderXYZ(SvxTGeo *tgeo, vecd &x, vecd &y, vecd &z);
TCanvas *DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt);
void DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2, TString opt = "sz");
TCanvas *DrawDiffsLinear(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2,
                         const char *name, const char *title, TString coord, double maxy=0);

void
GetLadderXYZ(SvxTGeo *tgeo, vecd &x, vecd &y, vecd &z)
{
  int ipt = 0;
  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      double xyz[3] = {0};
      tgeo->GetSensorXYZ(i,j,0,xyz);
      x.push_back(xyz[0]);
      y.push_back(xyz[1]);
      z.push_back(xyz[2]);

      if (false)
        Printf("xyz %f %f %f", x[ipt], y[ipt], z[ipt]);

      ipt++;
    }

  return;
}

TCanvas *
DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt)
{
  TCanvas *c = new TCanvas(name, title, 900, 900);
  TH1F *hf = c->DrawFrame(-20, -20, 20, 20, title);
  hf->SetXTitle("East                 x [cm]                 West");
  hf->GetXaxis()->CenterTitle();
  hf->SetYTitle("y [cm]");
  hf->GetYaxis()->CenterTitle();

  TLatex label;
  label.SetTextSize(0.02);
  double xyz[3] = {0};
  for (int i=0; i<geo->GetNLayers(); i++)
    for (int j=0; j<geo->GetNLadders(i); j++)
    {
      TPolyLine *s = geo->LadderOutlineXY(i,j);
      s->SetLineColor(opt.Contains("faint") ? kGray : kGray+2);
      s->SetLineWidth(opt.Contains("faint") ? 1 : 2);

      if (Locked(i,j))
      {
        s->SetLineColor(kOrange-3);
        s->SetLineWidth(2);
      }

      s->Draw("same");

      if (opt.Contains("dead") && Dead(i,j)) // Mark dead ladders
      {
        double xyz[3] = {0};
        geo->GetSensorXYZ(i,j,0,xyz);

        TLatex ltx;
        ltx.SetTextSize(0.04);
        ltx.SetTextAlign(22);
        ltx.SetTextColor(kRed-4);
        ltx.DrawLatex(xyz[0], xyz[1], "#times");
      }
      if (opt.Contains("faint"))
        continue;

      geo->GetSensorXYZ(i, j, 0, xyz);
      int horz = xyz[0] > 0 ? 1 : 3;
      int vert = xyz[1] > 0 ? 1 : 3;
      label.SetTextAlign(10*horz + vert);
      label.DrawLatex(xyz[0], xyz[1], Form(" %d ", j));
    }
  return c;
}

void
DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2, TString opt)
{
  int n = (int)x1.size();

  for (int i=0; i<n; i++)
  {
    double f = 20;
    double x = x1[i], y = y1[i];
    double dx = x2[i] - x;
    double dy = y2[i] - y;
    double ds = TMath::Sqrt(dx*dx + dy*dy);
    double dz = z2[i] - z1[i];

    // Draw points showing displacement in z coordinate
    TMarker mkr;
    mkr.SetMarkerStyle(kOpenCircle);
    mkr.SetMarkerColor(kGray+1);
    mkr.SetMarkerSize(0.5);
    //    mkr.DrawMarker(x,y);

    // Draw z displacements
    if (opt.Contains("z") && TMath::Abs(dz) > 5e-4) // Only label changes > 5 um
    {
      mkr.SetMarkerStyle(kOpenCircle);
      mkr.SetMarkerColor(dz>0 ? kRed-4 : kAzure-4); // Red/blue shift mnemonic
      mkr.SetMarkerSize(100*TMath::Abs(dz)); // 1 = 8px diam, 2 = 16px, ...
      mkr.DrawMarker(x,y);

      TLatex ltx;
      ltx.SetTextSize(0.018);
      ltx.SetTextAlign(22);
      ltx.SetTextColor(dz>0 ? kRed+2 : kAzure+3);
      ltx.DrawLatex(x, y, Form("%.0f", 1e4*dz));
    }

    // Draw arrows showing s displacements
    if (opt.Contains("s") && ds > 5e-4) // Only label changes > 5 um
    {
      TArrow a;
      a.SetLineWidth(2);
      a.DrawArrow(x, y, x + f*dx, y + f*dy, 0.005);

      TLatex ltx;
      ltx.SetTextSize(0.015);
      ltx.DrawLatex(x + f*dx, y + f*dy, Form("%.0f", 1e4*ds));
    }
  }

  return;
}

TCanvas *
DrawDiffsLinear(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2,
                const char *name, const char *title, TString coord, double maxy)
{
  int nLadders[4] = {10,20,16,24}; // ladders/layer
  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3

  int n = (int)x1.size();
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(43);
  ltx.SetTextSize(24);

  TGraphAsymmErrors *g = new TGraphAsymmErrors(n);
  g->SetMarkerStyle(kOpenCircle);
  g->SetMarkerColor(coord=="s" ? kRed+2 : kAzure+3);
  g->SetTitle(Form(";ladder index in layer {0,1,2,3};%s [cm]", 
              coord=="s" ? "#Deltas" : "#Deltaz"));
  double maxval = 0;
  for (int i=0; i<n; i++)
  {
    double x = x1[i], y = y1[i];
    double dx = x2[i] - x;
    double dy = y2[i] - y;
    double ds = TMath::Sqrt(dx*dx + dy*dy);
    if (TMath::ATan2(dy,dx) < 0)
      ds *= -1;

    double dz = z2[i] - z1[i];
    double val = (coord=="s") ? ds : dz;

    if (val > maxval)
      maxval = val;

    g->SetPoint(i, (double)i, val);
    if (val > 0.0)
      g->SetPointEYlow(i, val);
    if (val < 0.0)
      g->SetPointEYhigh(i, -val);

  }

  TCanvas *c = new TCanvas(name, title, 1400, 500);
  double xstart = 0.03;
  double x0 = xstart;

  g->Draw("ap");
  if (maxy > 0) maxval = maxy;
  g->GetYaxis()->SetRangeUser(-1.5*maxval, 1.5*maxval);
  g->GetYaxis()->SetTickLength(0.01);
  g->GetYaxis()->SetNdivisions(208);
  gPad->SetGridy();
  gPad->SetMargin(0.07, 0.02, 0.1, 0.01); // l,r,b,t

  g->GetXaxis()->Set(n, -0.5, n-0.5);
  for (int i=0; i<n; i++)
  {
    int layer = TMath::BinarySearch(4, start, i);
    int k = i - start[layer];
    g->GetXaxis()->SetBinLabel(i+1, Form("%d", k));

  }

  for (int layer=0; layer<4; layer++)
  {
    TLine l;
    l.DrawLine(start[layer]-0.5, -1.5*maxval, start[layer]-0.5, 1.5*maxval);
  }
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  return c;
}

#endif
