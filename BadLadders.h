#ifndef __BADLADDERS_H__
#define __BADLADDERS_H__

#include "VtxAlignBase.h"
#include "ParameterDefs.h"

bool Dead(int layer, int ladder);
bool RegularizedLadder(int layer, int ladder);
void RegularizedLadders(SvxTGeo* geo, vecs &gpars, veci &labels);

bool
Dead(int layer, int ladder)
{
  // return false; // TEMP - simulated data
  if (layer==1 && ladder==11) return true;
  if (layer==3 && ladder==10) return true;
  if (layer==3 && ladder==16) return true;
  if (layer==3 && ladder==23) return true;
  return false;
}

bool
RegularizedLadder(int layer, int ladder)
{
  // return false; // TEMP - simulated data
  if (layer==3 && ladder==13) return true; // Behind a dead pixel ladder
  if (layer==3 && ladder==17) return true; // Unstable, not sure why
  return false;
}

void
RegularizedLadders(SvxTGeo* geo, vecs &gpars, veci &labels)
{
  labels.clear();
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr); ldr++)
      for (unsigned int ic=0; ic<gpars.size(); ++ic)
        if (RegularizedLadder(lyr, ldr))
          labels.push_back(Label(lyr, ldr, gpars[ic]));
  return;
}

#endif
