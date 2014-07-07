#include "SvxTGeo.h"
#include "ParameterDefs.h"
#include "BadLadders.h"
#include "VtxVis.h"

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z);

void DiffGeometry(const char *pfa = "geom/svxPISA-411768.par.1",
                  const char *pfb = "geom/svxPISA-411768.par.2")
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

  SvxTGeo *geo = new SvxTGeo;
  geo->ReadParFile(pfa);
  geo->MakeTopVolume(100, 100, 100);
  geo->AddSensors();
  geo->GeoManager()->CloseGeometry();
  DrawXY(geo, "s_diff", "#Deltas_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "s");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",pfa,pfb));

  DrawXY(geo, "z_diff", "#Deltaz_{ab}", "L, dead, faint");
  DrawDiffs(xa, ya, za, xb, yb, zb, "z");
  ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}",pfa,pfb));

  return;
}

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z)
{
  SvxTGeo *geo = new SvxTGeo;
  geo->ReadParFile(parfile);
  geo->MakeTopVolume(100, 100, 100);
  geo->AddSensors();
  geo->GeoManager()->CloseGeometry();
  GetLadderXYZ(geo, x, y, z);
  delete geo;
  return;
}
