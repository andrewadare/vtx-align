#ifndef __BEAMCENTERFUNCS_H__
#define __BEAMCENTERFUNCS_H__

#include "GLSFitter.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"

typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

TVectorD BeamCenter(geoTracks &tracks, int ntrk, TString arm, TString opt="");
TVectorD BeamCenter(geoEvents &events, int ntrk, TString arm, TString opt="");
TVectorD IPVec(TVectorD &a, TVectorD &n, TVectorD &p);
TVectorD IPVec(SvxGeoTrack &t, TVectorD &p);
bool East(double phi);
void FillSystem(geoEvents &events, TMatrixD &M, TVectorD &y, TString arm);

TVectorD
Vertex(geoTracks &event, TString arm)
{
  // Compute least-squares vertex from tracks.
  int n = (int)event.size();
  TMatrixD M(n,3);   // "Design matrix" containing track slopes
  TVectorD y0(n);    // Vector of track y-intercept values
  TMatrixD L(n,n);   // Covariance matrix for track fit parameters y0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dy0 from track fits
  TMatrixD cov(3,3); // Assigned in SolveGLS()

  for (int i=0; i<n; i++)
  {
    SvxGeoTrack track = event.at(i);
    bool east = East(track.phi0);
    if (arm=="" || (arm=="east" && east) || (arm=="west" && !east))
    {
      // TODO: code this
      // M(i, 0) = -TMath::Tan(phi);
      // M(i, 1) = 1;
      // y0(i) = yint;
      // i++;
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
BeamCenter(geoTracks &tracks, int ntrk, TString arm, TString opt)
{
  // Compute least-squares beam center from ntrk tracks.
  if (0) Printf("%s", opt.Data());

  assert(ntrk <= (int)tracks.size());

  int n = ntrk; 
  TMatrixD M(n,2);   // "Design matrix" containing track slopes
  TVectorD y0(n);    // Vector of track y-intercept values
  TMatrixD L(n,n);   // Covariance matrix for track fit parameters y0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dy0 from track fits
  TMatrixD cov(2,2); // Assigned in SolveGLS()

  int row=0;
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    if (row==n)
      break;
    double yint = tracks[i].vy;
    double phi = tracks[i].phi0;
    bool east = East(phi);
    if ((arm=="east" && east) || (arm=="west" && !east))
    {
      M(row, 0) = -TMath::Tan(phi);
      M(row, 1) = 1;
      y0(row) = yint;
      row++;
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
BeamCenter(geoEvents &events, int ntrk, TString arm, TString opt)
{
  // Compute least-squares beam center from ntrk tracks.
  if (0) Printf("%s", opt.Data());

  int n = ntrk;
  TMatrixD M(n,2);   // "Design matrix" containing track slopes
  TVectorD y0(n);    // Vector of track y-intercept values
  TMatrixD L(n,n);   // Covariance matrix for track fit parameters y0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dy0 from track fits
  TMatrixD cov(2,2); // Assigned in SolveGLS()

  FillSystem(events, M, y0, arm);

  TVectorD bc = SolveGLS(M, y0, L, cov);
  Printf("%s x,y (%.3f +- %.3f, %.3f +- %.3f)",
         arm.Data(),
         bc(0), TMath::Sqrt(cov(0,0)),
         bc(1), TMath::Sqrt(cov(1,1)));
  cov.Print();
  return bc;
}

void
FillSystem(geoEvents &events, TMatrixD &M, TVectorD &y0, TString arm)
{
  int row = 0;
  int nrows = M.GetNrows();
  assert(y0.GetNrows()==nrows);

  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      if (row==nrows)
        return;
      double yint = events[ev][t].vy;
      double phi = events[ev][t].phi0;
      bool east = East(phi);
      if ((arm=="east" && east) || (arm=="west" && !east))
      {
        M(row, 0) = -TMath::Tan(phi);
        M(row, 1) = 1;
        y0(row) = yint;
        row++;
      }
    }

  // Function shouldn't reach this point....
  if (row < nrows)
    Printf("Warning! Linear system dimensions too large; not enough tracks.");

  return;
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

bool
East(double phi)
{
  return (phi > TMath::PiOver2() && phi < 3*TMath::PiOver2());
}

#endif
