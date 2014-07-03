#ifndef __BEAMCENTERFUNCS_H__
#define __BEAMCENTERFUNCS_H__

#include "GLSFitter.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"

typedef vector<SvxGeoTrack> geoTracks;
TVectorD BeamCenter(geoTracks &tracks, int ntrk, TString arm);
TVectorD IPVec(TVectorD &a, TVectorD &n, TVectorD &p);
TVectorD IPVec(SvxGeoTrack &t, TVectorD &p);

TVectorD
BeamCenter(geoTracks &tracks, int ntrk, TString arm)
{
  // Compute least-squares beam center from ntrk tracks.
  // Warning: TrackFitSResid() overwrites tracks[i].ds as a side effect.

  int n = TMath::Min(ntrk, (int)tracks.size());
  TMatrixD M(n,2);   // "Design matrix" containing track slopes
  TVectorD y0(n);    // Vector of track y-intercept values
  TMatrixD L(n,n);   // Covariance matrix for track fit parameters y0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dy0 from track fits
  TMatrixD cov(2,2); // Assigned in SolveGLS()

  for (int i=0; i<n; i++)
  {
    // Perform straight-line fit --> residuals, track parameters
    double pars[2] = {0}; /* y0, phi */
    TrackFitSResid(tracks[i], pars);
    double yint = pars[0], phi = pars[1];
    bool east = (phi > 0.5*TMath::Pi() && phi < 1.5*TMath::Pi());
    if ((arm=="east" && east) || (arm=="west" && !east))
    {
      M(i, 0) = -TMath::Tan(phi);
      M(i, 1) = 1;
      y0(i) = yint;
      i++;
    }
  }

  TVectorD bc = SolveGLS(M, y0, L, cov);
  Printf("%s x,y (%.3f +- %.3f, %.3f +- %.3f)",
         arm.Data(),
         bc(0), TMath::Sqrt(cov(0,0)),
         bc(1), TMath::Sqrt(cov(1,1)));
  cov.Print();
  return bc;
}

TVectorD
IPVec(TVectorD &a, TVectorD &n, TVectorD &p)
{
  // Compute impact parameter vector from point p to line x = a + tn
  // where
  // - x is a straight-line track trajectory
  // - a is a point on vector x (e.g. y-intercept point (0, y0))
  // - n is a unit vector directing x (e.g. (cos(phi), sin(phi)).
  return a - p - ((a - p)*n)*n;
}

TVectorD
IPVec(SvxGeoTrack &t, TVectorD &p)
{
  TVectorD a(2); a(1) = t.vy;
  TVectorD n(2);
  n(0) = TMath::Cos(t.phi0);
  n(1) = TMath::Sin(t.phi0);
  return IPVec(a,n,p);
}
#endif
