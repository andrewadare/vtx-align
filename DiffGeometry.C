#include "VtxAlignBase.h"

#include "ParameterDefs.h"
#include "BadLadders.h"
#include "VtxVis.h"

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z);

void DiffGeometry(const char *pfa = "geom/svxPISA-ideal.par",
                  const char *pfb = "geom/svxPISA-411768.par",
                  const char* dspdf = "pdfs/geom-ds-ideal-v-db.pdf",
                  const char* dzpdf = "pdfs/geom-dz-ideal-v-db.pdf")
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

  DrawXY(geo, "s_diff", "#Deltas_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "s");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",pfa,pfb));
  gPad->Print(dspdf);

  DrawXY(geo, "z_diff", "#Deltaz_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "z");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",pfa,pfb));
  gPad->Print(dzpdf);

  return;
}

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z)
{
  SvxTGeo *geo = VTXModel(parfile);
  GetLadderXYZ(geo, x, y, z);
  delete geo;
  return;
}
