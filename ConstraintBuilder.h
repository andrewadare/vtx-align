#ifndef __CONSTRAINTBUILDER_H__
#define __CONSTRAINTBUILDER_H__

#include "SvxTGeo.h"
#include "ParameterDefs.h"

#include <fstream>
#include <string>

typedef vector<double> vecd;
typedef vector<int>    veci;

bool In(int val, veci v);
void Ones(veci &labels, veci &xlabels, vecd &x);
void Radii(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x);
void PhiAngles(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x);
void RPhi(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x);
void AddConstraint(veci &labels, vecd &coords, ofstream &fs,
                   string comment = "", double sumto=0.0);

bool In(int val, veci v)
{
  for (unsigned int j=0; j<v.size(); j++)
    if (val==v[j])
      return true;
  return false;
}

void
Ones(veci &labels, veci &xlabels, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  for (unsigned int i=0; i<labels.size(); i++)
    x[i] = In(labels[i], xlabels) ? 0.0 : 1.0;
}

void
Radii(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    x[i] = In(labels[i], xlabels) ? 0.0 : geo->SensorRadius(lyr, ldr, 0);
  }
}

void
PhiAngles(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x)
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
    x[i] = In(labels[i], xlabels) ? 0.0 : phi;
  }
}

void
RPhi(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x)
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
    x[i] = In(labels[i], xlabels) ? 0.0 : rphi;
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
#endif
