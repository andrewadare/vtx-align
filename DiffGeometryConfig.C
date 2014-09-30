#include "VtxIO.h"

void DiffGeometryConfig()
{
  // Either type in offsets here...
  float e2wa[] = {-0.0183, -0.0071, -0.1475};
  float v2ca[] = {-0.45, 0.07, +0.1};//-0.3095};

  // Or get them out of a production config file. 
  // Or do both, as done in this script.
  float e2wb[3] = {0};
  float v2cb[3] = {0};
  GetOffsetsFromConfig("production/config/config-taebong-run14p10.txt",
                       e2wb,v2cb);

  gROOT->LoadMacro("DiffGeometry.C");

  DiffGeometry("geom/svxPISA-ideal.par",
               "geom/taebong-run14p10.par",
               e2wa, v2ca, e2wb, v2cb);
  return;
}
