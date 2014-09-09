#ifndef __MILLEFUNCTIONS_H__
#define __MILLEFUNCTIONS_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "ParameterDefs.h"

using namespace std;

void MilleVtx(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars,
              TGraphErrors *bc = 0);
void MilleCnt(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars,
              TGraphErrors *bc = 0);
float GlobalDerivative(SvxGeoTrack &trk, int ihit, string res, string par,
                       TGraphErrors *bc = 0);

void
MilleVtx(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars, TGraphErrors *bc)
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
      slabels.push_back(Label(hit.layer, hit.ladder, sgpars[k]));
      sdergl.push_back(GlobalDerivative(trk, j, "s", sgpars[k], bc));
    }
    for (int k=0; k<nzgp; k++)
    {
      zlabels.push_back(Label(hit.layer, hit.ladder, zgpars[k]));
      zdergl.push_back(GlobalDerivative(trk, j, "z", zgpars[k], bc));
    }

    // Assign sigmas for denominator of chi square function.
    // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
    // If millepede complains that chi^2/ndf is away from 1.0,
    // this is a good place to make adjustments.
    float sigs = 4. * hit.xsigma * ClusterXResolution(hit.layer);
    float sigz = 4. * hit.zsigma * ClusterZResolution(hit.layer);

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
// void
// MilleVtxStandalone(Mille &m, SvxGeoTrack &trk)
// {
//   // Mille::mille() is called once for each residual (so, twice per hit:
//   // one for the s coordinate, one for z).
//   // These calls fill a binary file (standalone.bin) with info on
//   // measured and expected residual sizes, and how the residuals change with
//   // respect to local (track) and global (detector position) parameters.
//   // The actual fitting is performed afterward by the pede program.

//   for (int j=0; j<trk.nhits; j++)
//   {
//     SvxGeoHit hit = trk.GetHit(j);
//     int ls = Label(hit.layer, hit.ladder, "s");
//     int lz = Label(hit.layer, hit.ladder, "z");
//     int lr = Label(hit.layer, hit.ladder, "r");
//     int slabels[2] = {ls, lr};
//     int zlabels[2] = {lz, lr};
//     float r = hit.x*hit.x + hit.y*hit.y;
//     //float theta = TMath::ATan2(r, hit.z);

//     // Assign sigmas for denominator of chi square function.
//     // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
//     // If millepede complains that chi^2/ndf is away from 1.0,
//     // this is a good place to make adjustments.
//     float sigs = 4. * hit.xsigma * ClusterXResolution(hit.layer);
//     float sigz = 4. * hit.zsigma * ClusterZResolution(hit.layer);

//     if (false)
//       Printf("hit.ds %.3g, sigs %.3g, hit.dz %.3g, sigz %.3g",
//              hit.ds, sigs, hit.dz, sigz);

//     // Local derivatives
//     // Split into two independent fits per track:
//     // 1. in the r-z plane: z = z0 + slope*r
//     // 2. in the x-y plane: y' = y0' + slope*r (where y' \perp r)
//     float sderlc[4] = {1.0,   r, 0.0, 0.0}; // dy(r)/dy0, dy(r)/dslope
//     float zderlc[4] = {0.0, 0.0, 1.0,   r}; // dz(r)/dz0, dz(r)/dslope

//     // Approximate global derivatives
//     float dsdr = hit.ds/r;
//     float dzdr = TMath::Tan(trk.the0 - TMath::PiOver2());

//     // Use ngd to control whether derivatives of residuals
//     // w.r.t. radial coordinate are included in global fit.
//     int ngd = allowRShifts ? 2 : 1;

//     float sdergl[2] = {1., dsdr}; // ds/ds, ds/dr
//     float zdergl[2] = {1., dzdr}; // dz/dz, dz/dr
//     m.mille(4, sderlc, ngd, sdergl, slabels, hit.ds, sigs);
//     m.mille(4, zderlc, ngd, zdergl, zlabels, hit.dz, sigz);
//   }

//   // Write residuals for this track to file & reset for next track.
//   m.end();

//   return;
// }

void
MilleCnt(Mille &m, SvxGeoTrack &trk, vecs &sgpars, vecs &zgpars, TGraphErrors *bc)
{
  // Get the number of global parameters for s and for z. This sets the number
  // of derivatives ds/d* and dz/d* to be computed.
  int nsgp = sgpars.size();
  int nzgp = zgpars.size();

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
      string res = "s";
      string par = sgpars[k];
      slabels.push_back(Label(hit.layer, hit.ladder, sgpars[k]));
      sdergl.push_back(GlobalDerivative(trk, j, res, par, bc));
    }
    for (int k=0; k<nzgp; k++)
    {
      string res = "z";
      string par = zgpars[k];
      zlabels.push_back(Label(hit.layer, hit.ladder, zgpars[k]));
      zdergl.push_back(GlobalDerivative(trk, j, res, par, bc));
    }

    // Assign sigmas for denominator of chi square function.
    // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
    // If millepede complains that chi^2/ndf is away from 1.0,
    // this is a good place to make adjustments.
    float sigs = 4. * hit.xsigma * ClusterXResolution(hit.layer);
    float sigz = 4. * hit.zsigma * ClusterZResolution(hit.layer);

    // Here, fit clusters individually. Each cluster is treated as
    // one "local fit object".

    // No global fit to s residual
    m.mille(1, sderlc, 0, &sdergl[0], &slabels[0], 0, .000001);
    m.mille(1, zderlc, nzgp, &zdergl[0], &zlabels[0], hit.dz, sigz);
    m.end();

    // No global fit to z residual
    m.mille(1, sderlc, nsgp, &sdergl[0], &slabels[0], hit.ds, sigs);
    m.mille(1, zderlc, 0, &zdergl[0], &zlabels[0], 0, .000001);
    m.end();
  }

  return;
}

float
GlobalDerivative(SvxGeoTrack &trk, int ihit, string res, string par,
                 TGraphErrors *bc)
{
  // Return first derivative of residual with respect to global parameter.
  float deriv = 0.0;
  double bcx = 0., bcy = 0.;
  if (bc) // Include beam position in track fit. Rotate about (bcx,bcy).
  {
    int arm = (trk.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.
    bcx = bc->GetX()[arm];
    bcy = bc->GetY()[arm];
  }

  SvxGeoHit hit = trk.GetHit(ihit);
  float r = hit.x*hit.x + hit.y*hit.y;
  float xp = bcx + r*TMath::Cos(trk.phi0);
  float yp = bcy + r*TMath::Sin(trk.phi0) + trk.vy; // vy is y-intercept from fit not vertex
  double dx = xp - hit.x;
  double dy = yp - hit.y;
  double dxy = TMath::Sqrt(dx*dx + dy*dy);

  if (res == "s")
  {
    // d(Delta_s)/dr
    if (par == "r")
      return hit.ds/r;

    // d(Delta_s)/ds
    if (par == "s")
      return 1.0;

    // d(Delta_s)/dx
    if (par == "x")
      return dx/dxy;

    // d(Delta_s)/dy
    if (par == "y")
      return dy/dxy;

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
  return deriv;
}

#endif
