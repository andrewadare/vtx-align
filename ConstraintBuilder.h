#ifndef __CONSTRAINTBUILDER_H__
#define __CONSTRAINTBUILDER_H__

#include "SvxTGeo.h"
#include "ParameterDefs.h"
#include "BadLadders.h"

#include <fstream>
#include <string>

typedef vector<double> vecd;
typedef vector<int>    veci;

template<typename T> bool In(T val, vector<T> &v);
veci LadderLabels(SvxTGeo *geo, string arm, string coord);
void AddConstraint(veci &labels, vecd &coords, ofstream &fs,
                   string comment = "", double sumto = 0.0);
void AddConstraint(veci &labels, SvxTGeo *geo,
                   void(*fl)(veci &, SvxTGeo *, vecd &),
                   ofstream &fs, string comment = "", double sumto = 0.0);

// Position querying callback functions: store current positions into x for the
// provided labels.
// These are used to provide the f_l factors in the global parameter
// sum constraint equations of the form sum_l(f_l*p_l) = const.
// See sec. 6 of the Millepede manual for some documentation.
// They should all return void and take the same parameter list.
void      Ones(veci &labels, SvxTGeo *geo, vecd &x);
void     Radii(veci &labels, SvxTGeo *geo, vecd &x);
void PhiAngles(veci &labels, SvxTGeo *geo, vecd &x);
void      RPhi(veci &labels, SvxTGeo *geo, vecd &x);

template<typename T> 
bool In(T val, vector<T> &v)
{
  // Check for membership in a vector.
  // Just a less verbose wrapper for std::find().
  if (find(v.begin(), v.end(), val) != v.end())
    return true;
  return false;
}

veci
LadderLabels(SvxTGeo *geo, string arm, string coord)
{
  veci l;
  if (arm == "W" || arm == "w" || arm == "all")    // West arm
  {
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
        if (!Dead(lyr,ldr) && !Locked(lyr,ldr))
          l.push_back(Label(lyr, ldr, coord));
  }
  if (arm == "E" || arm == "e" || arm == "all")    // East arm
  {
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
        if (!Dead(lyr,ldr) && !Locked(lyr,ldr))
          l.push_back(Label(lyr, ldr, coord));
  }
  return l;
}

void
Ones(veci &labels, SvxTGeo * /*geo*/, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  for (unsigned int i=0; i<labels.size(); i++)
    x[i] = 1.0;
}

void
Radii(veci &labels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    x[i] = geo->SensorRadius(lyr, ldr, 0);
  }
}

void
PhiAngles(veci &labels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    float phi = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0),
                                        TMath::TwoPi());
    x[i] = phi;
  }
}

void
RPhi(veci &labels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    float phi = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0),
                                        TMath::TwoPi());
    float rphi = phi*geo->SensorRadius(lyr, ldr, 0);
    x[i] = rphi;
  }
}

void
AddConstraint(veci &labels, vecd &coords, ofstream &fs, string comment, double sumto)
{
  fs << "Constraint " << sumto << "  ! " << comment << endl;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    int lyr, ldr;
    string s;
    ParInfo(labels[i], lyr, ldr, s);
    fs << labels[i] << " " << coords[i]
       << " ! B" << lyr << "L" << ldr << "(" << s << ")" << endl;
  }
  fs << endl;

  return;
}

void
AddConstraint(veci &labels, SvxTGeo *geo,
              void(*fl)(veci &, SvxTGeo *, vecd &),
              ofstream &fs, string comment, double sumto)
{
  vecd x(labels.size(), 0.0);
  fl(labels, geo, x);

  AddConstraint(labels,x,fs,comment,sumto);
  return;
}

#endif
