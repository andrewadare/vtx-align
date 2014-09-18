#include "VtxAlignBase.h"

#include "ParameterDefs.h"
#include "BadLadders.h"
#include "VtxVis.h"
#include "UtilFns.h"

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z);

void DiffGeometry(const char *pfa = "geom/svxPISA-ideal.par",
                  const char *pfb = "geom/taebong-run14p10.par")
{
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);
  ltx.SetTextSize(0.03);
  TObjArray *cList = new TObjArray();

  // Get ladder positions in "a" and "b" geometry
  vecd xa; vecd ya; vecd za;
  vecd xb; vecd yb; vecd zb;
  GetXYZ(pfa, xa, ya, za);
  GetXYZ(pfb, xb, yb, zb);

  SvxTGeo *geo = VTXModel(pfa);

  TString a = TString(gSystem->BaseName(pfa)).ReplaceAll(".par", "");
  TString b = TString(gSystem->BaseName(pfb)).ReplaceAll(".par", "");
  TString ab = Form("%s-vs-%s", a.Data(), b.Data());
  const char *fa = gSystem->BaseName(pfa);
  const char *fb = gSystem->BaseName(pfb);

  DrawXY(geo, "s_diff", "#Delta(x,y)_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "s");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/dsxy%s.pdf",ab.Data()));
  cList->Add(gPad);

  DrawXY(geo, "z_diff", "#Deltaz_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "z");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/dzxy%s.pdf",ab.Data()));
  cList->Add(gPad);

  ltx.SetTextSize(0.05);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "xdiff", "#Delta x", "x");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/dx%s.pdf",ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "ydiff", "#Delta y", "y");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/dy%s.pdf",ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "zdiff", "#Delta z", "z");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/dz%s.pdf",ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "sdiff", "#Delta s", "s");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/ds%s.pdf",ab.Data()));
  cList->Add(gPad);

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "rdiff", "#Delta r", "r");
  ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}",fa,fb));
  gPad->Print(Form("pdfs/dr%s.pdf",ab.Data()));
  cList->Add(gPad);

  PrintPDF(cList, Form("pdfs/%s", ab.Data()), "");
  return;
}

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z)
{
  SvxTGeo *geo = VTXModel(parfile);
  GetLadderXYZ(geo, x, y, z);
  delete geo;
  return;
}
