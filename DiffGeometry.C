#include "VtxAlignBase.h"

#include "ParameterDefs.h"
#include "BadLadders.h"
#include "VtxVis.h"

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z);

void DiffGeometry(const char *pfa = "geom/411768-2-0.par",
                  const char *pfb = "geom/411768-2-3.par")
{
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);
  ltx.SetTextSize(0.03);

  // Get ladder positions in "a" and "b" geometry
  vecd xa; vecd ya; vecd za;
  vecd xb; vecd yb; vecd zb;
  GetXYZ(pfa, xa, ya, za);
  GetXYZ(pfb, xb, yb, zb);

  SvxTGeo *geo = VTXModel(pfa);

  TString a = TString(gSystem->BaseName(pfa)).ReplaceAll(".par", "");
  TString b = TString(gSystem->BaseName(pfb)).ReplaceAll(".par", "");
  const char *ab = Form("%s-vs-%s", a.Data(), b.Data());

  DrawXY(geo, "s_diff", "#Deltas_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "s");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",pfa,pfb));
  gPad->Print(Form("pdfs/dsxy%s.pdf",ab));

  DrawXY(geo, "z_diff", "#Deltaz_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "z");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",pfa,pfb));
  gPad->Print(Form("pdfs/dzxy%s.pdf",ab));

  ltx.SetTextSize(0.05);
  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "sdiff", "#Delta s", "s", 0.083);
  ltx.DrawLatex(0.7, 0.92, Form("#splitline{a: %s}{b: %s}",pfa,pfb));
  gPad->Print(Form("pdfs/ds%s.pdf",ab));

  DrawDiffsLinear(xa, ya, za, xb, yb, zb, "zdiff", "#Delta z", "z", 0.028);
  ltx.DrawLatex(0.7, 0.92, Form("#splitline{a: %s}{b: %s}",pfa,pfb));
  gPad->Print(Form("pdfs/dz%s.pdf",ab));

  return;
}

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z)
{
  SvxTGeo *geo = VTXModel(parfile);
  GetLadderXYZ(geo, x, y, z);
  delete geo;
  return;
}
