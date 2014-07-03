#ifndef __PARAMETERDEFS_H__
#define __PARAMETERDEFS_H__

#include <string>

int Label(int layer, int ladder, string coord);
void ParInfo(int label, int &layer, int &ladder, string &coord);

int
Label(int layer, int ladder, string coord)
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
  // Add other coordinates / degrees of freedom here as needed
  // Then ParInfo() would need corresponding modification.

  return 70*ic + start[layer] + ladder + 1;
}

void
ParInfo(int label, int &layer, int &ladder, string &coord)
{
  // Get layer, ladder, and coordinate string from global parameter label.
  // Inverse of Label() function.

  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3
  int l = 0;

  if (label > 70)
  {
    coord = "z";
    l = label - 70;
  }
  else
  {
    coord = "s";
    l = label;
  }

  layer = 3;
  for (int i=3; i>=0; i--)
    if (l < start[i] + 1)
      layer = i - 1;

  ladder = l - start[layer] - 1;

  return;
}
#endif
