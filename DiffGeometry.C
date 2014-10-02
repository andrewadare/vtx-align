#include "VtxAlignBase.h"

#include "ParameterDefs.h"
#include "BadLadders.h"
#include "VtxVis.h"
#include "UtilFns.h"

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z,
            float e2w[], float v2c[]);

void DiffGeometry(const char *pfa = "geom/svxPISA-ideal.par",
                  const char *pfb = "geom/taebong-run14p10.par",
                  float *e2wa = 0,
                  float *v2ca = 0,
                  float *e2wb = 0,
                  float *v2cb = 0)
{
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);
  ltx.SetTextSize(0.03);
  TObjArray *cList = new TObjArray();

  // Get ladder positions in "a" and "b" geometry
  vecd xa; vecd ya; vecd za;
  vecd xb; vecd yb; vecd zb;
  GetXYZ(pfa, xa, ya, za, e2wa, v2ca);
  GetXYZ(pfb, xb, yb, zb, e2wb, v2cb);

  SvxTGeo *geo = VTXModel(pfa);

  // Move ladders to coincide with arrow tails and circle centers
  if (e2wa)
  {
    geo->TranslateArm(0, e2wa[0], e2wa[1], e2wa[2]); // East arm only
  }
  if (v2ca)
  {
    geo->TranslateArm(0, v2ca[0], v2ca[1], v2ca[2]);
    geo->TranslateArm(1, v2ca[0], v2ca[1], v2ca[2]);
  }

  TString a = TString(gSystem->BaseName(pfa)).ReplaceAll(".par", "");
  TString b = TString(gSystem->BaseName(pfb)).ReplaceAll(".par", "");
  TString ab = Form("%s-vs-%s", a.Data(), b.Data());
  const char *fa = gSystem->BaseName(pfa);
  const char *fb = gSystem->BaseName(pfb);

  DrawXY(geo, "s_diff", "#Delta(x,y)_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "s");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}", fa, fb));
  ltx.SetTextSize(0.02);
  if (v2ca)
    ltx.DrawLatex(0.02, 0.98,
                  Form("v2c a: (%.4f, %.4f, %.4f)", v2ca[0], v2ca[1], v2ca[2]));
  if (v2cb)
    ltx.DrawLatex(0.02, 0.96,
                  Form("v2c b: (%.4f, %.4f, %.4f)", v2cb[0], v2cb[1], v2cb[2]));
  if (e2wa)
    ltx.DrawLatex(0.02, 0.94,
                  Form("e2w a: (%.4f, %.4f, %.4f)", e2wa[0], e2wa[1], e2wa[2]));
  if (e2wb)
    ltx.DrawLatex(0.02, 0.92,
                  Form("e2w b: (%.4f, %.4f, %.4f)", e2wb[0], e2wb[1], e2wb[2]));
  gPad->Print(Form("pdfs/dsxy%s.pdf", ab.Data()));
  cList->Add(gPad);

  DrawXY(geo, "z_diff", "#Deltaz_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "z");
  ltx.SetTextSize(0.03);
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}", fa, fb));
  ltx.SetTextSize(0.02);
  if (v2ca)
    ltx.DrawLatex(0.02, 0.98,
                  Form("v2c a: (%.4f, %.4f, %.4f)", v2ca[0], v2ca[1], v2ca[2]));
  if (v2cb)
    ltx.DrawLatex(0.02, 0.96,
                  Form("v2c b: (%.4f, %.4f, %.4f)", v2cb[0], v2cb[1], v2cb[2]));
  if (e2wa)
    ltx.DrawLatex(0.02, 0.94,
                  Form("e2w a: (%.4f, %.4f, %.4f)", e2wa[0], e2wa[1], e2wa[2]));
  if (e2wb)
    ltx.DrawLatex(0.02, 0.92,
                  Form("e2w b: (%.4f, %.4f, %.4f)", e2wb[0], e2wb[1], e2wb[2]));
  gPad->Print(Form("pdfs/dzxy%s.pdf", ab.Data()));
  cList->Add(gPad);

  ltx.SetTextSize(0.05);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "xdiff", "#Delta x", "x");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", fa, fb));
  /*
  //if desired, here is a good setup
  ltx.DrawLatex(0.48, 0.92, Form("#splitline{a: %s}{b: %s}", fa, fb));
  if (v2ca)
    ltx.DrawLatex(0.22, 0.945,
                  Form("v2c a: (%.4f, %.4f, %.4f)", v2ca[0], v2ca[1], v2ca[2]));
  if (v2cb)
    ltx.DrawLatex(0.22, 0.895,
                  Form("v2c b: (%.4f, %.4f, %.4f)", v2cb[0], v2cb[1], v2cb[2]));
  if (e2wa)
    ltx.DrawLatex(0.70, 0.945,
                  Form("e2w a: (%.4f, %.4f, %.4f)", e2wa[0], e2wa[1], e2wa[2]));
  if (e2wb)
    ltx.DrawLatex(0.70, 0.895,
                  Form("e2w b: (%.4f, %.4f, %.4f)", e2wb[0], e2wb[1], e2wb[2]));
  */
  gPad->Print(Form("pdfs/dx%s.pdf", ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "ydiff", "#Delta y", "y");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", fa, fb));
  gPad->Print(Form("pdfs/dy%s.pdf", ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "zdiff", "#Delta z", "z");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", fa, fb));
  gPad->Print(Form("pdfs/dz%s.pdf", ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "sdiff", "#Delta s", "s");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", fa, fb));
  gPad->Print(Form("pdfs/ds%s.pdf", ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "rdiff", "#Delta r", "r");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", fa, fb));
  gPad->Print(Form("pdfs/dr%s.pdf", ab.Data()));
  cList->Add(gPad);

  PrintPDF(cList, Form("pdfs/%s", ab.Data()), "");
  return;
}

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z,
            float e2w[], float v2c[])
{
  SvxTGeo *geo = VTXModel(parfile);

  // Apply E-W and VTX-CNT offsets, if provided.
  if (e2w)
  {
    geo->TranslateArm(0, e2w[0], e2w[1], e2w[2]); // East arm only
  }
  if (v2c)
  {
    geo->TranslateArm(0, v2c[0], v2c[1], v2c[2]);
    geo->TranslateArm(1, v2c[0], v2c[1], v2c[2]);
  }

  GetLadderXYZ(geo, x, y, z);
  delete geo;
  return;
}
