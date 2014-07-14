#ifndef __BEAMCENTERFUNCS_H__
#define __BEAMCENTERFUNCS_H__

#include "GLSFitter.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "DcaFunctions.h"

typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

TVectorD BeamCenter(geoTracks &tracks, int ntrk, TString arm, TString opt="");
TVectorD BeamCenter(geoEvents &events, int ntrk, TString arm, TString opt="");
void FillSystem(geoEvents &events, TMatrixD &M, TVectorD &y, TString arm);

TVectorD
BeamCenter(geoTracks &tracks, int ntrk, TString arm, TString opt)
{
  return XYCenter(tracks, arm, ntrk, opt);
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
  // TODO: Add logic to allow filling from both arms (arm=="")
  //       See GLSFitter.h XYCenter().
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

  // Function shouldn't reach this point.
  // If it does, either decrease ntrk or provide larger tracks vector.
  // Or, (TODO) Consider providing option to on-the-fly resize M and y0,
  //       e.g. M.ResizeTo(row, M.GetNcols());
  //            y0.ResizeTo(row);
  if (row < nrows)
    Printf("BeamCenterFuncs.h FillSystem(): Warning!! "
           "Not enough tracks: %d requested, %d available.", nrows, row);

  return;
}

#endif
