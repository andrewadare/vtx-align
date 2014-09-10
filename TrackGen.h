#ifndef __TRACKGEN_H__
#define __TRACKGEN_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"

void GenerateTrack(SvxTGeo *geo, TVectorD &vertex, SvxGeoTrack &t, double B=0);
void GenerateEvent(SvxTGeo *geo, int ntracks, TVectorD &vertex, geoTracks &tracks);
void AddHitNoise(SvxGeoTrack &t, TVectorD &xsigma, TVectorD &zsigma);
bool TrackOk(SvxGeoTrack &t);

void
GenerateEvent(SvxTGeo *geo, int ntracks, TVectorD &vertex, geoTracks &tracks)
{

  while ((int)tracks.size() < ntracks)
  {
    SvxGeoTrack trk;
    GenerateTrack(geo, vertex, trk);

    // Inefficient...could move AddHitNoise() to its own loop later
    TVectorD xsig(4);
    TVectorD zsig(4);
    for (int j=0; j<4; j++)
    {
      xsig(j) = ClusterXResolution(j);
      zsig(j) = ClusterZResolution(j);
    }
    AddHitNoise(trk, xsig, zsig);

    // trk.Print();

    if (TrackOk(trk))
      tracks.push_back(trk);
  }

  // Printf("%lu tracks generated.", tracks.size());

  return;
}

void
GenerateTrack(SvxTGeo *geo, TVectorD &vertex, SvxGeoTrack &t, double B)
{
  static TRandom3 ran3;
  SvxProj proj;
  proj.SetVerbosity(0);
  static long trackid = 0;

  t.nhits  = 0;
  t.charge = 0;
  t.vx     = vertex(0);
  t.vy     = vertex(1);
  t.vz     = vertex(2);
  t.mom    = ran3.Uniform(0.5, 5.0);
  t.phi0   = ran3.Uniform(0,TMath::TwoPi());
  t.the0   = ran3.Uniform(TMath::PiOver4(), 3*TMath::PiOver4());
  t.bfield = B;

  if (false)
    Printf("Initialized track at (vx,vy,vz) %.2f,%.2f,%.2f with "
           "(phi0,the0) %.2f, %.2f charge %d, magField %.2f, mom %.2f",
           t.vx, t.vy, t.vz, t.phi0, t.the0, t.charge, t.bfield, t.mom);

  // Add hits to track
  proj.FindHitsFromVertex(t, geo);

  for (int j=0; j<t.nhits; j++)
    t.hits[j].trkid = trackid;
  trackid++;

  return;
}

void
AddHitNoise(SvxGeoTrack &t, TVectorD &xsigma, TVectorD &zsigma)
{
  static TRandom3 ran3;
  for (int j=0; j<t.nhits; j++)
  {
    SvxGeoHit hit = t.hits[j];
    // Store cluster size as an integer (channel count?), as done now in data.
    // Set to 1 or 2 as a rough approximation of the data.
    t.hits[j].xsigma = (ran3.Uniform() < 0.67) ? 1 : 2;
    t.hits[j].zsigma = (ran3.Uniform() < 0.67) ? 1 : 2;
    double dx = t.hits[j].xsigma*xsigma(hit.layer)*ran3.Gaus();
    double dz = t.hits[j].xsigma*zsigma(hit.layer)*ran3.Gaus();

    // t.hits[j].xsigma = xsigma(hit.layer);
    // t.hits[j].zsigma = zsigma(hit.layer);
    // double dx = t.hits[j].xsigma*ran3.Gaus();
    // double dz = t.hits[j].zsigma*ran3.Gaus();

    t.hits[j].xs += dx;
    t.hits[j].zs += dz;

    t.hits[j].ds = dx;
    t.hits[j].dz = dz;

    double lxyz[3] = {hit.xs, hit.ys, hit.zs};
    double gxyz[3] = {0};
    TGeoNode *s = hit.node;
    if (s)
      s->GetMatrix()->LocalToMaster(lxyz, gxyz);

    t.hits[j].x = gxyz[0];
    t.hits[j].y = gxyz[1];
    t.hits[j].z = gxyz[2];
  }
  return;
}

bool
TrackOk(SvxGeoTrack &t)
{
  if (t.nhits == 4)
  {
    if (t.hits[0].layer==0 &&
        t.hits[1].layer==1 &&
        t.hits[2].layer==2 &&
        t.hits[3].layer==3)
    {
      return true;
    }
  }
  return false;
}

#endif
