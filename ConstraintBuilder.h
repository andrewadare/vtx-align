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
veci EdgeLadderLabels(SvxTGeo *geo, string arm, string coord, string t_or_b);
void AddConstraint(veci &labels, vecd &coords, ofstream &fs,
                   string comment = "", double sumto = 0.0);
void AddConstraint(veci &labels, SvxTGeo *geo,
                   void(*fl)(veci &, SvxTGeo *, vecd &),
                   ofstream &fs, string comment = "", double sumto = 0.0);
void AddConstraints(veci &wlabels, veci &elabels, SvxTGeo *geo,
                    void(*fl)(veci &, SvxTGeo *, vecd &),
                    ofstream &fs, string comment, double sumto = 0.0);

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
void      XPos(veci &labels, SvxTGeo *geo, vecd &x);
void      YPos(veci &labels, SvxTGeo *geo, vecd &x);

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
EdgeLadderLabels(SvxTGeo *geo, string arm, string coord, string t_or_b)
{
  // Boundary "wedge" units - sequences of 4 edge ladders at top and bottom.
  // Select "top" or "bottom" using t_or_b parameter.
  int wt[4] = {4,9,7,11};
  int wb[4] = {0,0,0,0};
  int et[4] = {5,10,8,12};
  int eb[4] = {9,19,15,23};
  veci v;

  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
  {
    if ((arm == "W" || arm == "w") && (t_or_b == "top"))
      v.push_back(Label(lyr, wt[lyr], coord));
    if ((arm == "W" || arm == "w") && (t_or_b == "bottom"))
      v.push_back(Label(lyr, wb[lyr], coord));
    if ((arm == "E" || arm == "e") && (t_or_b == "top"))
      v.push_back(Label(lyr, et[lyr], coord));
    if ((arm == "E" || arm == "e") && (t_or_b == "bottom"))
      v.push_back(Label(lyr, eb[lyr], coord));
  }

  if (v.size() == 0)
    Printf("WARNING from EdgeLadderLabels(): No labels returned.");

  return v;
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
XPos(veci &labels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    double xyz[3] = {0};
    geo->GetSensorXYZ(lyr,ldr,0,xyz);
    x[i] = xyz[0];
  }
}

void
YPos(veci &labels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    double xyz[3] = {0};
    geo->GetSensorXYZ(lyr,ldr,0,xyz);
    x[i] = xyz[1];
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

void
AddConstraints(veci &wlabels, veci &elabels, SvxTGeo *geo,
               void(*fl)(veci &, SvxTGeo *, vecd &),
               ofstream &fs, string comment, double sumto)
{
  vecd wx(wlabels.size(), 0.0);
  vecd ex(elabels.size(), 0.0);

  fl(wlabels, geo, wx);
  AddConstraint(wlabels, wx, fs, string("West " + comment), sumto);

  fl(elabels, geo, ex);
  AddConstraint(elabels, ex, fs, string("East " + comment), sumto);

  return;
}

#endif
