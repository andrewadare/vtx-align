#ifndef __BADLADDERS_H__
#define __BADLADDERS_H__

bool Dead(int layer, int ladder);
bool Locked(int layer, int ladder);

bool
Dead(int layer, int ladder)
{
  return false; // TEMP - simulated data
  if (layer==1 && ladder==11) return true;
  if (layer==3 && ladder==10) return true;
  if (layer==3 && ladder==16) return true;
  if (layer==3 && ladder==23) return true;
  return false;
}

bool
Locked(int layer, int ladder)
{
  return false; // TEMP - simulated data
  if (layer==3 && ladder==13) return true; // 66% // Behind a dead pixel ladder==15
  return false;
}

#endif
