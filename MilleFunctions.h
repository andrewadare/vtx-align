#ifndef __MILLEFUNCTIONS_H__
#define __MILLEFUNCTIONS_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "ParameterDefs.h"

using namespace std;

void MilleVtx(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars,
              TGraphErrors *bc = 0, TString opt = "");
void MilleCnt(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars,
              TGraphErrors *bc = 0, TString opt = "");
float GlobalDerivative(SvxGeoTrack &trk, int ihit, string res, string par,
                       TGraphErrors *bc = 0);
void
MilleVtx(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars,
         TGraphErrors *bc, TString opt)
{
  // Mille::mille() is called once for each residual (so, twice per hit:
  // one for the s coordinate, one for z).
  // These calls fill a binary file (standalone.bin) with info on
  // measured and expected residual sizes, and how the residuals change with
  // respect to local (track) and global (detector position) parameters.
  // The actual fitting is performed afterward by the pede program.

  // Get the number of global parameters for s and for z. This sets the number
  // of derivatives ds/d* and dz/d* to be computed.
  int nsgp = sgpars.size();
  int nzgp = zgpars.size();
  int arm = (trk.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.

  for (int j=0; j<trk.nhits; j++)
  {
    SvxGeoHit hit = trk.GetHit(j);

    // Fill vectors of labels and global derivatives for s and z residuals
    veci slabels;
    veci zlabels;
    vector<float> sdergl;
    vector<float> zdergl;
    for (int k=0; k<nsgp; k++)
    {
      if (opt.Contains("halflayer"))
        slabels.push_back(HalfLayerLabel(hit.layer, arm, sgpars[k]));
      else
        slabels.push_back(Label(hit.layer, hit.ladder, sgpars[k]));
      sdergl.push_back(GlobalDerivative(trk, j, "s", sgpars[k], bc));
    }
    for (int k=0; k<nzgp; k++)
    {
      if (opt.Contains("halflayer"))
        zlabels.push_back(HalfLayerLabel(hit.layer, arm, zgpars[k]));
      else
        zlabels.push_back(Label(hit.layer, hit.ladder, zgpars[k]));
      zdergl.push_back(GlobalDerivative(trk, j, "z", zgpars[k], bc));
    }

    // Assign sigmas for denominator of chi square function.
    // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
    // If millepede complains that chi^2/ndf is away from 1.0,
    // this is a good place to make adjustments.
    float sigs = hit.xsigma * ClusterXResolution(hit.layer);
    float sigz = hit.zsigma * ClusterZResolution(hit.layer);

    if (false)
      Printf("hit.ds %.3g, sigs %.3g, hit.dz %.3g, sigz %.3g",
             hit.ds, sigs, hit.dz, sigz);

    // Local derivatives
    // Split into two independent fits per track:
    // 1. in the r-z plane: z = z0 + slope*r
    // 2. in the x-y plane: y' = y0' + slope*r (where y' \perp r)
    float r = hit.x*hit.x + hit.y*hit.y;
    float sderlc[4] = {1.0,   r, 0.0, 0.0}; // dy(r)/dy0, dy(r)/dslope
    float zderlc[4] = {0.0, 0.0, 1.0,   r}; // dz(r)/dz0, dz(r)/dslope

    m.mille(4, sderlc, nsgp, &sdergl[0], &slabels[0], hit.ds, sigs);
    m.mille(4, zderlc, nzgp, &zdergl[0], &zlabels[0], hit.dz, sigz);
  }

  // Write residuals for this track to file & reset for next track.
  m.end();

  return;
}

void
MilleCnt(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars,
         TGraphErrors *bc, TString opt)
{
  // Get the number of global parameters for s and for z. This sets the number
  // of derivatives ds/d* and dz/d* to be computed.
  int nsgp = sgpars.size();
  int nzgp = zgpars.size();
  int arm = (trk.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.

  // Local derivatives
  float sderlc[1] = {1.0};
  float zderlc[1] = {1.0};

  for (int j=0; j<trk.nhits; j++)
  {
    SvxGeoHit hit = trk.GetHit(j);

    // Fill vectors of labels and global derivatives for s and z residuals
    veci slabels;
    veci zlabels;
    vector<float> sdergl;
    vector<float> zdergl;

    for (int k=0; k<nsgp; k++)
    {
      if (opt.Contains("halflayer"))
        slabels.push_back(HalfLayerLabel(hit.layer, arm, sgpars[k]));
      else
        slabels.push_back(Label(hit.layer, hit.ladder, sgpars[k]));
      sdergl.push_back(GlobalDerivative(trk, j, "s", sgpars[k], bc));
    }
    for (int k=0; k<nzgp; k++)
    {
      if (opt.Contains("halflayer"))
        zlabels.push_back(HalfLayerLabel(hit.layer, arm, zgpars[k]));
      else
        zlabels.push_back(Label(hit.layer, hit.ladder, zgpars[k]));
      zdergl.push_back(GlobalDerivative(trk, j, "z", zgpars[k], bc));
    }

    // Assign sigmas for denominator of chi square function.
    // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
    // If millepede complains that chi^2/ndf is away from 1.0,
    // this is a good place to make adjustments.
    float sigs = 1.6*hit.xsigma * ClusterXResolution(hit.layer);
    float sigz = 1.6*hit.zsigma * ClusterZResolution(hit.layer);

    // Here, fit clusters individually. Each cluster is treated as
    // one "local fit object".

    m.mille(0, sderlc, nsgp, &sdergl[0], &slabels[0], hit.ds, sigs);
    m.mille(1, sderlc,    0, &sdergl[0], &slabels[0],      0, sigs);
    m.end();

    m.mille(0, zderlc, nzgp, &zdergl[0], &zlabels[0], hit.dz, sigz);
    m.mille(1, zderlc,    0, &zdergl[0], &zlabels[0],      0, sigz);
    m.end();
  }

  return;
}

float
GlobalDerivative(SvxGeoTrack &trk, int ihit, string res, string par,
                 TGraphErrors *bc)
{
  // Return first derivative of residual with respect to global parameter.
  double bcx = 0., bcy = 0.;
  if (bc) // Include beam position in track fit. Rotate about (bcx,bcy).
  {
    int arm = (trk.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.
    bcx = bc->GetX()[arm];
    bcy = bc->GetY()[arm];
  }

  SvxGeoHit hit = trk.GetHit(ihit);
  float r = TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);
  if (res == "s")
  {
    // d(Delta_s)/dr
    if (par == "r")
      return hit.ds/r;

    // d(Delta_s)/ds -- 0.974 accounts for 13 degree tilt in pixel layers
    if (par == "s")
      return hit.layer < 2 ? 0.9741 : 1.0;

    // d(Delta_s)/dx
    if (par == "x")
      return hit.layer < 2 ? -0.9741*hit.y/r : -hit.y/r;

    // d(Delta_s)/dy
    if (par == "y")
      return hit.layer < 2 ? 0.9741*hit.x/r : hit.x/r;

    // d(Delta_s)/dz
    if (par == "z")
      return 0.0;
  }

  if (res == "z")
  {
    // d(Delta_z)/dr
    if (par == "r")
      return 1./TMath::Tan(trk.the0);

    // d(Delta_z)/ds
    if (par == "s")
      return 0.0;

    // d(Delta_z)/dx
    if (par == "x")
      return hit.x/r/TMath::Tan(trk.the0);

    // d(Delta_z)/dy
    if (par == "y")
      return hit.y/r/TMath::Tan(trk.the0);

    // d(Delta_z)/dz
    if (par == "z")
      return 1.0;
  }

  Printf("WARNING from GlobalDerivative() in MilleFunctions.h: "
         "No derivative d(Delta_%s)/d%s", res.c_str(), par.c_str());
  return 0.0;
}

#endif
