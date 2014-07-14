#ifndef __VERTEXFINDER_H__
#define __VERTEXFINDER_H__

#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "GLSFitter.h"

typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;
TVectorD Vertex(geoTracks &event, TString arm);
double ZVertex(geoTracks &event, TString arm);

// Code snippet for median calculation:
// root [0] double x[] = {1,2,3,4,5,6,7,8,9,10};
// root [1] double probs[] = {0.25, 0.5, 0.75};
// root [2] double quantiles[3] = {0};
// root [3] TMath::Quantiles(10, 3, x, quantiles, probs, false)
// root [4] quantiles
// (double [3]) { 3.250000e+00, 5.500000e+00, 7.750000e+00 }

double
ZVertex(geoTracks &tracks, TString arm)
{
  int nz = 0;
  const int maxnz = tracks.size();
  assert(maxnz>0);
  double zs[maxnz];
  double probs[1] = {0.5}; // For median = 50% quantile
  double quantiles[1] = {0};

  for (unsigned int i=0; i<tracks.size(); i++)
  {
    bool east = East(tracks[i].phi0);
    if (arm=="")
      zs[nz++] = tracks[i].vz;
    else if ((arm=="east" && east) || (arm=="west" && !east))
      zs[nz++] = tracks[i].vz;
    // else
    //   Printf("VertexFinder.h ZVertex(): Unknown option %s", arm.Data());
  }
  TMath::Quantiles(nz, 1, zs, quantiles, probs, false);

  return quantiles[0];
}

TVectorD
Vertex(geoTracks &event, TString arm)
{
  // Compute least-squares vertex from tracks.

  int n = event.size();
  assert(n>0);
  TVectorD xy = XYCenter(event, arm);

  TVectorD xyz(3);
  xyz(0) = xy(0);
  xyz(1) = xy(1);
  xyz(2) = ZVertex(event, arm);
  return xyz;
}


#endif
