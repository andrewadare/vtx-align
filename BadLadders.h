#ifndef __BADLADDERS_H__
#define __BADLADDERS_H__

bool Dead(int layer, int ladder);
bool Locked(int layer, int ladder);

bool
Dead(int layer, int ladder)
{
  if (layer==1 && ladder==11) return true;
  if (layer==3 && ladder==10) return true;
  if (layer==3 && ladder==16) return true;
  if (layer==3 && ladder==23) return true;
  return false;
}

bool
Locked(int layer, int ladder)
{
  // if (layer==2 && ladder==15) return true; // 50%
  // if (layer==3 && ladder==2)  return true; // 65%
  // //if (layer==3 && ladder==12) return true; // 70%
  if (layer==3 && ladder==13) return true; // 66% // Behind a dead pixel ladder==15
  // if (layer==3 && ladder==22) return true; // Not dead, but problematic

  // // if (layer==3 && ladder==9)  return true; // 66%
  // // if (layer==3 && ladder==15) return true; // 66%
  // // if (layer==3 && ladder==19) return true; // 68%
  // if (layer==3 && ladder==21) return true; // 50%
  return false;
}

#endif
