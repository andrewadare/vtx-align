#include "ParameterDefs.h"
#include <vector>
#include <string>

#include "TROOT.h"

void CheckLadderPars();
void CheckHalfLayerPars();

void CheckParLabels()
{
  Printf("Ladder parameters:");
  Printf("  par | lyr ldr coord");
  CheckLadderPars();

  Printf("\n\n\n\n\n");
  Printf("Half-layer parameters:");
  CheckHalfLayerPars();
}

void
CheckLadderPars()
{
  int nLadders[4] = {10,20,16,24}; // ladders/layer

  vector<std::string> coords;
  coords.push_back("x");
  coords.push_back("y");
  coords.push_back("z");
  coords.push_back("s");
  coords.push_back("r");

  int lyr_back = -1, ldr_back = -1;
  string coord_back = "";

  for (int ic=0; ic<5; ic++)
    for (int lyr=0; lyr<4; lyr++)
      for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      {
        int l = Label(lyr, ldr, coords[ic]);
        ParInfo(l, lyr_back, ldr_back, coord_back);

        Printf("%5d | %d %3d %5s <--> %d %3d %5s",
               l, lyr, ldr, coords[ic].c_str(),
               lyr_back, ldr_back, coord_back.c_str());
      }
  return;
}

void
CheckHalfLayerPars()
{
  vector<std::string> coords;
  coords.push_back("x");
  coords.push_back("y");
  coords.push_back("z");
  coords.push_back("phi");

  int lyr_back = -1, arm_back = -1;
  string coord_back = "";

  for (int lyr=0; lyr<4; lyr++)
    for (int arm=0; arm<2; arm++)
      for (int ic=0; ic<4; ic++)
      {
        int l = HalfLayerLabel(lyr, arm, coords[ic]);
        HalfLayerParInfo(l, lyr_back, arm_back, coord_back);

        Printf("%d | %d %d %5s <--> %d %d %5s",
               l, lyr, arm, coords[ic].c_str(),
               lyr_back, arm_back, coord_back.c_str());
      }

  return;
}