#include "VtxIO.h"

void DiffGeometryConfig()
{
  // To assign the optional east to west and VTX to CNT offsets, either
  // type in the offsets directly (as done for geometry "a"),
  // or read them from a production config file (as done for "b"),
  // or use a combination of both.
  // float e2wa[3] = {0};
  // float v2ca[3] = {0};
  float e2wa[] = { -0.0183, -0.0071, -0.1475};
  float v2ca[] = { -0.45, 0.07, +0.1};
  float bcb[2] = {0};
  float e2wb[3] = {0};
  float v2cb[3] = {0};
  string geob = "";

  GetParamFromConfig("production/config/config-taebong-run14p10.txt",
                       bcb, e2wb, v2cb, geob);

  gROOT->LoadMacro("DiffGeometry.C");

  //D. McGlinchey - This must be wrong, but doesn't seem to work otherwise
  const char *geofileb = geob.c_str();
  DiffGeometry("geom/svxPISA-ideal.par",
               geofileb,
               e2wa, v2ca, e2wb, v2cb);
  return;
}
