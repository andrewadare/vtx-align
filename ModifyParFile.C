#include "VtxAlignBase.h"
#include "ParameterDefs.h"

void ModifyParFile(const char *pfa = "geom/svxPISA-hubert.par",
                   const char *pfb = "geom/svxPISA-hubert-mod.par")
{
  SvxTGeo *geo = VTXModel(pfa);

  geo->MoveLadderRadially(2, 2, -0.2950);
  geo->MoveLadderRadially(2, 5, -0.2950);
  geo->MoveLadderRadially(2, 10, -0.2950);
  geo->MoveLadderRadially(2, 13, -0.2950);

  Printf("Writing %s", pfb);
  geo->WriteParFile(pfb);

  return;
}
