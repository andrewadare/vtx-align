#ifndef __PARAMETERDEFS_H__
#define __PARAMETERDEFS_H__

#include <string>

int Label(int layer, int ladder, string coord);
int HalfLayerLabel(int layer, int arm, string coord);
void ParInfo(int label, int &layer, int &ladder, string &coord);
void HalfLayerParInfo(int label, int &layer, int &arm, string &coord);

int
Label(int layer, int ladder, string &coord)
{
  // Return a global parameter label for this layer, ladder, and coordinate.
  // In Millepede II, any unique integer > 0 will do.
  // Labels are not required to be sequential.
  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3
  int ic = 99;

  if (coord.compare("s") == 0)
    ic = 0;
  if (coord.compare("z") == 0)
    ic = 1;
  if (coord.compare("r") == 0)
    ic = 2;

  // Add other coordinates / degrees of freedom here as needed
  // Then ParInfo() would need corresponding modification.

  return 70*ic + start[layer] + ladder + 1;
}

int
HalfLayerLabel(int layer, int arm, string coord)
{
  // Return a global parameter label for this half-layer and coordinate.
  // Half-layer labels shouldn't conflict with ladder labels, since those
  // range into the hundreds and these have a +1000 offset.
  // Expecting layer \in {0-3}, arm \in {0=E,1=W} and coord \in {x,y,z,phi}.
  int nHalfLayers = 8;
  int ic = 99;

  if (coord.compare("x") == 0)
    ic = 0;
  if (coord.compare("y") == 0)
    ic = 1;
  if (coord.compare("z") == 0)
    ic = 2;
  if (coord.compare("phi") == 0)
    ic = 3;

  return 1000 + 100*layer + 10*arm + ic;
}

void
ParInfo(int label, int &layer, int &ladder, string &coord)
{
  // Get layer, ladder, and coordinate string from global parameter label.
  // Inverse of Label() function.

  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3

  // Assign layer
  int l = label % 70; // Global ladder index 1,2,...,69
  if (l == 0) l = 70; // Handle corner cases: label = 70, 140, 210 --> l = 70
  layer = 3;
  for (int i=3; i>=0; i--)
    if (l < start[i] + 1)
      layer = i - 1;

  // Assign ladder
  ladder = l - start[layer] - 1;

  // Assign coord
  coord = (label <= 70)? "s" : (label <= 140)? "z" : (label <= 210)? "r" : "";

  return;
}

void
HalfLayerParInfo(int label, int &layer, int &arm, string &coord)
{
  const char* coords[4] = {"x", "y", "z", "phi"};

  layer = (label%1000 - label%100)/100;
  arm   = (label%100 - label%10)/10;
  coord = string(coords[label%10]);
  return;
}

#endif
