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

TCanvas *DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt);
void DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2);

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
DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2)
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

    if (TMath::Abs(dz) > 5e-4) // Only label changes > 5 um
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

    // Draw arrows showing (significant) displacements in xy plane
    if (ds > 5e-4) // Only label changes > 5 um
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

#endif
