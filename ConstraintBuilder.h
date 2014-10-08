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
template<typename T> int Idx(T val, vector<T> &v);
veci LadderLabels(SvxTGeo *geo, string arm, string coord);
veci LadderLabels(SvxTGeo *geo, int arm, int layer, string coord);
veci EdgeLadderLabels(SvxTGeo *geo, string arm, string coord, string t_or_b);
veci HalfLayerLabels(SvxTGeo *geo, string arm, string coord);

void AddConstraint(veci &labels, vecd &coords, ofstream &fs,
                   string comment = "", double sumto = 0.0);
void AddConstraint(veci &labels, SvxTGeo *geo,
                   void(*fl)(veci &, SvxTGeo *, vecd &),
                   ofstream &fs, string comment = "", double sumto = 0.0);
void AddConstraints(veci &wlabels, veci &elabels, SvxTGeo *geo,
                    void(*fl)(veci &, SvxTGeo *, vecd &),
                    ofstream &fs, string comment, double sumto = 0.0);
void WriteLadderConstraints(const char *filename,
                            vecs &sgpars, vecs &zgpars,
                            vecd &sgpre, vecd &zgpre,
                            SvxTGeo *geo, TString opt = "");
void WriteHLConstraints(const char *filename,
                        vecs &sgpars, vecs &zgpars,
                        vecd &sgpre, vecd &zgpre,
                        SvxTGeo *geo, TString opt = "");
// void ApplyPreSigma(ofstream &fs, SvxTGeo *geo, vecs &gpars, double presigma);
void ApplyPreSigma(ofstream &fs, veci &labels, double presigma);

// Position querying callback functions: store current positions into x for the
// provided labels.
// These are used to provide the f_l factors in the global parameter
// sum constraint equations of the form sum_l(f_l*p_l) = const.
// See sec. 6 of the Millepede manual for some documentation.
// They should all return void and take the same parameter list.
void      Ones(veci &labels, SvxTGeo *geo, vecd &x);
void     Radii(veci &labels, SvxTGeo *geo, vecd &x);
void   HLRadii(veci &labels, SvxTGeo *geo, vecd &x);
void PhiAngles(veci &labels, SvxTGeo *geo, vecd &x);
void      RPhi(veci &labels, SvxTGeo *geo, vecd &x);
void      XPos(veci &labels, SvxTGeo *geo, vecd &x);
void      YPos(veci &labels, SvxTGeo *geo, vecd &x);

template<typename T>
bool In(T val, vector<T> &v)
{
  // Check for membership in a vector.
  if (find(v.begin(), v.end(), val) != v.end())
    return true;
  return false;
}

template<typename T>
int Idx(T val, vector<T> &v)
{
  // Return position of element in vector.
  // If element is not found, index will be v.size().
  // vector<T>::iterator it = find(begin(v), end(v), val);
  return int(find(begin(v), end(v), val) - begin(v));
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
        if (!Dead(lyr,ldr))
          l.push_back(Label(lyr, ldr, coord));
  }
  if (arm == "E" || arm == "e" || arm == "all")    // East arm
  {
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      for (int ldr=geo->GetNLadders(lyr)/2; ldr<geo->GetNLadders(lyr); ldr++)
        if (!Dead(lyr,ldr))
          l.push_back(Label(lyr, ldr, coord));
  }
  return l;
}

veci
LadderLabels(SvxTGeo *geo, int arm, int layer, string coord)
{
  veci l;
  int first = -1, last  = -1;
  geo->LadderRange(layer, arm, first, last);
  for (int ldr = first; ldr <= last; ldr++)
    if (!Dead(layer,ldr))
      l.push_back(Label(layer, ldr, coord));

  return l;
}

veci
HalfLayerLabels(SvxTGeo *geo, string arm, string coord)
{
  veci l;
  if (arm == "E" || arm == "e" || arm == "all")    // East arm (0)
  {
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      l.push_back(HalfLayerLabel(lyr, 0, coord));
  }
  if (arm == "W" || arm == "w" || arm == "all")    // West arm (1)
  {
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      l.push_back(HalfLayerLabel(lyr, 1, coord));
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
HLRadii(veci &labels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, arm;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    HalfLayerParInfo(labels[i], lyr, arm, s);

    // Get mean radius over ladders in halflayer
    int first = -1, last  = -1;
    geo->LadderRange(lyr, arm, first, last);
    int nladders = 0;
    for (int ldr = first; ldr <= last; ldr++)
    {
      x[i] += geo->SensorRadius(lyr, ldr, 0);
      ++nladders;
    }
    x[i] /= nladders;
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

void
WriteHLConstraints(const char *filename,
                   vecs &sgpars, vecs &zgpars,
                   vecd &sgpre, vecd &zgpre,
                   SvxTGeo *geo, TString opt)
{
  cout << "Writing " << filename << "..." << flush;
  ofstream fs(filename);

  if (opt.Contains("empty"))
  {
    fs.close();
    return;
  }

  bool xdof = In(string("x"), sgpars) || In(string("x"), zgpars);
  bool ydof = In(string("y"), sgpars) || In(string("y"), zgpars);
  bool zdof = In(string("z"), zgpars);
  bool sdof = In(string("s"), sgpars);

  fs << "! * Half-layer constraints file for input to pede *" << endl;

  veci hlx = HalfLayerLabels(geo, "all", "x");
  veci hly = HalfLayerLabels(geo, "all", "y");
  veci hlz = HalfLayerLabels(geo, "all", "z");
  veci hls = HalfLayerLabels(geo, "all", "s");

  veci ex = HalfLayerLabels(geo, "e", "x");
  veci ey = HalfLayerLabels(geo, "e", "y");
  veci ez = HalfLayerLabels(geo, "e", "z");
  veci es = HalfLayerLabels(geo, "e", "s");

  veci wx = HalfLayerLabels(geo, "w", "x");
  veci wy = HalfLayerLabels(geo, "w", "y");
  veci wz = HalfLayerLabels(geo, "w", "z");
  veci ws = HalfLayerLabels(geo, "w", "s");

  for (unsigned int ic=0; ic<sgpars.size(); ++ic)
  {
    veci ll = HalfLayerLabels(geo, "all", sgpars[ic]);
    ApplyPreSigma(fs, ll, sgpre.at(ic));
  }
  for (unsigned int ic=0; ic<zgpars.size(); ++ic)
  {
    veci ll = HalfLayerLabels(geo, "all", zgpars[ic]);
    ApplyPreSigma(fs, ll, zgpre.at(ic));
  }

  if (xdof)
  {
    // AddConstraint(hlx, geo, Ones, fs, "Total x translation");
    // AddConstraints(wx, ex, geo, Ones,  fs, "x translation");
    AddConstraints(wx, ex, geo, HLRadii, fs, "x r shear");
  }
  if (ydof)
  {
    // AddConstraint(hly, geo, Ones, fs, "Total y translation");
    // AddConstraints(wy, ey, geo, Ones,  fs, "y translation");
    AddConstraints(wy, ey, geo, HLRadii, fs, "y r shear");
  }
  if (zdof)
  {
    // AddConstraint(hlz, geo, Ones, fs, "Total z translation");
    // AddConstraints(wz, ez, geo, Ones,  fs, "z translation");
    AddConstraints(wz, ez, geo, HLRadii, fs, "z r shear");
  }
  if (sdof)
  {
    // AddConstraint(hls, geo, Ones, fs, "Total s translation");
    // AddConstraints(ws, es, geo, Ones,  fs, "s translation");
    AddConstraints(ws, es, geo, HLRadii, fs, "s r shear");
  }

  fs.close();
  Printf("done.");
  return;
}

void
ApplyPreSigma(ofstream &fs, veci &labels, double presigma)
{
  // Fix or regulate global parameters in the provided label set.
  fs << "! Presigma for global parameter regularization:" << endl;
  fs << "!  = 0: free." << endl;
  fs << "!  > 0: regularized. Smaller --> stronger regularization." << endl;
  fs << "!  < 0: fixed." << endl;
  fs << "Parameter ! Columns: label, initial value, presigma " << endl;

  // Currently leaving initial global par. at default value of zero.
  for (unsigned int i=0; i<labels.size(); ++i)
    fs << Form("%d 0.0 %g", labels[i], presigma) << endl;
  fs << endl;

  return;
}

void
WriteLadderConstraints(const char *filename,
                       vecs &sgpars, vecs &zgpars,
                       vecd &sgpre, vecd &zgpre,
                       SvxTGeo *geo, TString opt)
{
  cout << "Writing " << filename << "..." << flush;
  ofstream fs(filename);

  if (opt.Contains("empty"))
  {
    fs.close();
    return;
  }

  bool sdof = In(string("s"), sgpars);
  bool zdof = In(string("z"), zgpars);
  bool xdof = In(string("x"), sgpars) || In(string("x"), zgpars);
  bool ydof = In(string("y"), sgpars) || In(string("y"), zgpars);
  bool rdof = In(string("r"), sgpars) || In(string("r"), zgpars);

  fs << "! * Ladder constraints file for input to pede *" << endl;
  fs << "! Write linear constraints such that sum_l (f_l * p_l) = c" << endl;
  fs << "! where the f_l factors provide coord. dependence to constrain shear." << endl;
  fs << "! The format is as follows:" << endl;
  fs << "! Constraint c" << endl;
  fs << "! 123 1.0   ! Add a p_123 * f_123 term" << endl;
  fs << "! 234 2.345 ! Add a p_234 * f_234 term" << endl;
  fs << "! ..." << endl;
  fs << endl;

  // Global parameter labels for z and s coordinates in east and west arms
  veci wx = LadderLabels(geo, "w", "x");
  veci wy = LadderLabels(geo, "w", "y");
  veci wz = LadderLabels(geo, "w", "z");
  veci ws = LadderLabels(geo, "w", "s");
  veci wr = LadderLabels(geo, "w", "r");
  veci ex = LadderLabels(geo, "e", "x");
  veci ey = LadderLabels(geo, "e", "y");
  veci ez = LadderLabels(geo, "e", "z");
  veci es = LadderLabels(geo, "e", "s");
  veci er = LadderLabels(geo, "e", "r");

  veci wtz = EdgeLadderLabels(geo, "w", "z", "top");
  veci wts = EdgeLadderLabels(geo, "w", "s", "top");
  veci wbz = EdgeLadderLabels(geo, "w", "z", "bottom");
  veci wbs = EdgeLadderLabels(geo, "w", "s", "bottom");
  veci etz = EdgeLadderLabels(geo, "e", "z", "top");
  veci ets = EdgeLadderLabels(geo, "e", "s", "top");
  veci ebz = EdgeLadderLabels(geo, "e", "z", "bottom");
  veci ebs = EdgeLadderLabels(geo, "e", "s", "bottom");

  // ApplyPreSigma(fs, slabels_bad, 1e-4);
  // ApplyPreSigma(fs, zlabels_bad, 1e-4);

  if (!opt.Contains("sim"))
  {
    // Regularize a hand picked set of poorly-constrained ladders.
    // (low-stats, low-eff, faulty electronics, behind a dead ladder, etc).
    veci slabels_bad;
    veci zlabels_bad;
    RegularizedLadders(geo, sgpars, slabels_bad);
    RegularizedLadders(geo, zgpars, zlabels_bad);
    ApplyPreSigma(fs, slabels_bad, 1e-4);
    ApplyPreSigma(fs, zlabels_bad, 1e-4);
  }

  for (unsigned int ic=0; ic<sgpars.size(); ++ic)
  {
    veci ll = LadderLabels(geo, "all", sgpars[ic]);
    ApplyPreSigma(fs, ll, sgpre.at(ic));
  }
  for (unsigned int ic=0; ic<zgpars.size(); ++ic)
  {
    veci ll = LadderLabels(geo, "all", zgpars[ic]);
    ApplyPreSigma(fs, ll, zgpre.at(ic));
  }

  // Prevent net displacements/distortions, but only if the coordinate
  // is being used as a global parameter (otherwise the system will be
  // rank-deficient)
  if (sdof)
  {
    AddConstraints(ws, es, geo, Ones,      fs, "s translation");
    // AddConstraints(ws, es, geo, PhiAngles, fs, "s phi shear");
    // AddConstraints(ws, es, geo, RPhi,      fs, "s r-phi shear");
    // AddConstraints(ws, es, geo, Radii,     fs, "s r shear");
    AddConstraints(wts, ets, geo, Radii,   fs, "top s r shear");
    AddConstraints(wbs, ebs, geo, Radii,   fs, "bottom s r shear");
  }
  if (xdof)
  {
    AddConstraints(wx, ex, geo, Ones,      fs, "x translation");
    AddConstraints(wx, ex, geo, PhiAngles, fs, "x phi shear");
    // AddConstraints(wx, ex, geo, RPhi,      fs, "x r-phi shear");
    // AddConstraints(wx, ex, geo, Radii,     fs, "x r shear");
    // for (int lyr=0; lyr<4; lyr++)
    // {
    //   veci e = LadderLabels(geo, 0, lyr, "x");
    //   veci w = LadderLabels(geo, 1, lyr, "x");
    //   AddConstraints(w,e,geo,Ones,fs,Form("x translation (layer %d)", lyr));
    // }
  }
  if (ydof)
  {
    AddConstraints(wy, ey, geo, Ones,      fs, "y translation");
    AddConstraints(wy, ey, geo, PhiAngles, fs, "y phi shear");
    // AddConstraints(wy, ey, geo, RPhi,      fs, "y r-phi shear");
    // AddConstraints(wy, ey, geo, Radii,     fs, "y r shear");
    // for (int lyr=0; lyr<4; lyr++)
    // {
    //   veci e = LadderLabels(geo, 0, lyr, "y");
    //   veci w = LadderLabels(geo, 1, lyr, "y");
    //   AddConstraints(w,e,geo,Ones,fs,Form("y translation (layer %d)", lyr));
    // }
  }
  if (zdof)
  {
    AddConstraints(wz, ez, geo, Ones,      fs, "z translation");
    AddConstraints(wz, ez, geo, PhiAngles, fs, "z phi shear");
    // AddConstraints(wz, ez, geo, RPhi,      fs, "z r-phi shear");
    AddConstraints(wz, ez, geo, Radii,     fs, "z r shear");
    AddConstraints(wtz, etz, geo, Radii,   fs, "top z r shear");
    AddConstraints(wbz, ebz, geo, Radii,   fs, "bottom z r shear");
  }
  if (rdof)
  {
    AddConstraints(wr, er, geo, Ones,      fs, "r expansion/contraction");
  }

  fs.close();
  Printf("done.");
  return;
}

#endif
