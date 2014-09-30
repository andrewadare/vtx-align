#include "VtxIO.h"

void DiffGeometryConfig()
{
  // To assign the optional east to west and VTX to CNT offsets, either
  // type in the offsets directly (as done for geometry "a"),
  // or read them from a production config file (as done for "b"),
  // or use a combination of both.
  float e2wa[] = { -0.0183, -0.0071, -0.1475};
  float v2ca[] = { -0.45, 0.07, +0.1};

  float e2wb[3] = {0};
  float v2cb[3] = {0};
  GetOffsetsFromConfig("production/config/config-taebong-run14p10.txt",
                       e2wb, v2cb);

  string tmp = "geom/";
  tmp += GetGeomFromConfig("production/config/config-taebong-run14p10.txt");
  const char *geofile = tmp.c_str();

  gROOT->LoadMacro("DiffGeometry.C");

  DiffGeometry("geom/svxPISA-ideal.par",
               geofile,
               e2wa, v2ca, e2wb, v2cb);
  return;
}
