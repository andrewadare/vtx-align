#ifndef __VTXIO_H__
#define __VTXIO_H__

#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include <TNtuple.h>
#include <TLeaf.h>
#include <vector>

typedef vector<SvxGeoTrack> geoTracks;

void GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks);
void FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple);

void
GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks)
{
  // This function reads hit ntuple variables of the form
  // "layer:ladder:sensor:xs:ys:zs:x:y:z:xsigma:zsigma:dz:ds:trkid"
  // into the tracks vector.

  assert(t);
  assert(geo);

  long nentries = t->GetEntries();
  int previd = -1;

  Printf("Reading tracks from NTuple (%d entries)...", (int)nentries);
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int id = t->GetLeaf("trkid")->GetValue();
    if (id != previd) // New track
    {
      SvxGeoTrack t;
      // Track info is not stored in tree. If it were, it would go here.
      tracks.push_back(t);
    }

    SvxGeoHit hit;
    hit.layer  = t->GetLeaf("layer")->GetValue();
    hit.ladder = t->GetLeaf("ladder")->GetValue();
    hit.sensor = t->GetLeaf("sensor")->GetValue();
    hit.xs     = t->GetLeaf("lx")->GetValue();
    hit.ys     = t->GetLeaf("ly")->GetValue();
    hit.zs     = t->GetLeaf("lz")->GetValue();
    hit.x      = t->GetLeaf("gx")->GetValue();
    hit.y      = t->GetLeaf("gy")->GetValue();
    hit.z      = t->GetLeaf("gz")->GetValue();
    hit.xsigma = t->GetLeaf("x_size")->GetValue();
    hit.zsigma = t->GetLeaf("z_size")->GetValue();
    hit.dz     = t->GetLeaf("res_z")->GetValue();
    hit.ds     = t->GetLeaf("res_s")->GetValue();
    hit.trkid  = t->GetLeaf("trkid")->GetValue();
    hit.node   = geo->SensorNode(hit.layer, hit.ladder, hit.sensor);
    tracks.back().nhits++;
    tracks.back().hits.push_back(hit);
    previd = id;
  }

  Printf("%d tracks imported.", (int)tracks.size());
  return;
}

void
FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple)
{
  // TNtuple *ntuple = new TNtuple("t", "SvxGeoHit variables",
  // "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid"


  assert(ntuple);

  if (false)
    Printf("ds (%8.4f,%8.4f,%8.4f,%8.4f),  "
           "dz (%8.4f,%8.4f,%8.4f,%8.4f)",
           gt.hits[0].ds, gt.hits[1].ds, gt.hits[2].ds, gt.hits[3].ds,
           gt.hits[0].dz, gt.hits[1].dz, gt.hits[2].dz, gt.hits[3].dz);

  for (int ihit=0; ihit<gt.nhits; ihit++)
  {
    int nj = 14;
    std::vector<float> vars(nj, 0.);
    int j = 0;
    SvxGeoHit hit = gt.GetHit(ihit);

    vars[j++] = hit.layer  ;
    vars[j++] = hit.ladder ;
    vars[j++] = hit.sensor ;
    vars[j++] = hit.xs     ;
    vars[j++] = hit.ys     ;
    vars[j++] = hit.zs     ;
    vars[j++] = hit.x      ;
    vars[j++] = hit.y      ;
    vars[j++] = hit.z      ;
    vars[j++] = hit.xsigma ;
    vars[j++] = hit.zsigma ;
    vars[j++] = hit.dz     ;
    vars[j++] = hit.ds     ;
    vars[j++] = hit.trkid  ;
    ntuple->Fill(&vars[0]);
  }

  return;
}
#endif