#ifndef __VTXIO_H__
#define __VTXIO_H__

#include "VtxAlignBase.h"
#include <TLeaf.h>
#include <utility>
#include <map>

using namespace std;

void GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &trks, int nmax=-1);
void GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &trks, multimap<int, int> &eventsync, int nmax=-1,
		       TString opt = "");
void GetEventsFromTree(TNtuple *t, SvxTGeo *geo, geoEvents &evts, int nmax=-1,
                       TString opt = "");
void FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple, int event = -1);
void FillNTuple(geoEvents &events, TNtuple *ntuple);
int GetCorrections(const char *resFile, std::map<int, double> &mpc);
void WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfiles,
                    double regFactor = 1.0, double preSigma = 0.01);

void
GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks, int nmax)
{
  // This function reads hit ntuple variables of the form
  // "layer:ladder:sensor:xs:ys:zs:x:y:z:xsigma:zsigma:dz:ds:trkid"
  // into the tracks vector.

  assert(t);
  assert(geo);

  long nentries = t->GetEntries();
  int previd = -1;
  int nhits  = 0;

  cout << "Reading tracks from " << t->GetName()
       << " (" << nentries << " available hits)..."
       << flush;

  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int id = t->GetLeaf("trkid")->GetValue();
    if (id != previd) // New track
    {
      if (nmax > 0 && (int)tracks.size() == nmax)
      {
        Printf("%d tracks imported (%d hits).",
               (int)tracks.size(), nhits);
        return;
      }
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

    // Overwrite x,y,z with local-->global transformation on xs,ys,zs
    hit.Update();

    tracks.back().nhits++;
    tracks.back().hits.push_back(hit);
    previd = id;
    nhits++;
  }

  Printf("%d tracks imported (%d hits).",
         (int)tracks.size(), nhits);

  return;
}

void
GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks, std::multimap<int, int> &eventsync , int nmax, TString opt)
{
  // This function reads hit ntuple variables of the form
  // "layer:ladder:sensor:xs:ys:zs:x:y:z:xsigma:zsigma:dz:ds:trkid"
  // into the tracks vector.

  assert(t);
  assert(geo);

  long nentries = t->GetEntries();
  int previd = -1;
  int nhits  = 0;
  int itrack = 0;

  cout << "Reading tracks from " << t->GetName()
       << " (" << nentries << " available hits)..."
       << flush;

  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int id = t->GetLeaf("trkid")->GetValue();
    if (id != previd) // New track
    {
      if (nmax > 0 && (int)tracks.size() == nmax)
      {
        Printf("%d tracks imported (%d hits).",
               (int)tracks.size(), nhits);
        return;
      }
      int event = t->GetLeaf("event")->GetValue();
      eventsync.insert(pair<int,int>(itrack,event));
      itrack++;

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

    // Overwrite x,y,z with local-->global transformation on xs,ys,zs
    hit.Update();

    tracks.back().nhits++;
    tracks.back().hits.push_back(hit);
    previd = id;
    nhits++;
  }

  Printf("%d tracks imported (%d hits).",
         (int)tracks.size(), nhits);

  return;
}


void
GetEventsFromTree(TNtuple *t, SvxTGeo *geo, geoEvents &events, int nmax,
                  TString opt)
{
  // This function reads hit ntuple variables of the form
  // "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid:event"
  // into the tracks vector.
  assert(t);
  assert(geo);

  long nentries = t->GetEntries();
  int previd  = -1;
  int prevev  = -1;
  int ntracks = 0;
  int nhits   = 0;

  Printf("Reading events from NTuple (%d available hits)...", (int)nentries);
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);

    int ev = t->GetLeaf("event")->GetValue();
    int id = t->GetLeaf("trkid")->GetValue();

    bool newevent = (ev != prevev); // New event
    bool newtrack = (id != previd); // New track

    if (newevent && !newtrack)
      Error("GetEventsFromTree() in VtxIO.h",
            "Event/track ID problem !!! event %d track %d", ev, id);

    if (ev != prevev) // New event
    {
      if (nmax > 0 && (int)events.size() == nmax)
      {
        Printf("%d events imported (%d tracks, %d hits).",
               (int)events.size(), ntracks, nhits);
        return;
      }

      geoTracks tracks;
      events.push_back(tracks);
    }
    if (id != previd) // New track
    {
      SvxGeoTrack trk;
      if (opt.Contains("cnt"))
      {
        trk.phi0 = t->GetLeaf("trkphi")->GetValue();
        trk.the0 = t->GetLeaf("trktheta")->GetValue();
      }
      events.back().push_back(trk);
      ntracks++;
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
    hit.xsigma = t->GetLeaf("x_size")->GetValue(); // N.B. this is not sigma!
    hit.zsigma = t->GetLeaf("z_size")->GetValue(); // N.B. this is not sigma!
    hit.dz     = t->GetLeaf("res_z")->GetValue();
    hit.ds     = t->GetLeaf("res_s")->GetValue();
    hit.trkid  = t->GetLeaf("trkid")->GetValue();
    hit.node   = geo->SensorNode(hit.layer, hit.ladder, hit.sensor);

    // Overwrite x,y,z with local-->global transformation on xs,ys,zs
    hit.Update();

    events.back().back().nhits++;
    events.back().back().hits.push_back(hit);
    previd = id;
    prevev = ev;
    nhits++;
  }

  Printf("%d events imported (%d tracks, %d hits).",
         (int)events.size(), ntracks, nhits);

  return;
}

void
FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple, int event)
{
  // TNtuple *ntuple = new TNtuple("t", "SvxGeoHit variables",
  // "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:"
  // "res_z:res_s:trkid:event"

  assert(ntuple);


  if (false)
    Printf("ds (%8.4f,%8.4f,%8.4f,%8.4f),  "
           "dz (%8.4f,%8.4f,%8.4f,%8.4f)",
           gt.hits[0].ds, gt.hits[1].ds, gt.hits[2].ds, gt.hits[3].ds,
           gt.hits[0].dz, gt.hits[1].dz, gt.hits[2].dz, gt.hits[3].dz);

  for (int ihit=0; ihit<gt.nhits; ihit++)
  {
    int nj = 15;
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
    vars[j++] = event      ;
    ntuple->Fill(&vars[0]);
  }

  return;
}

void
FillNTuple(geoEvents &events, TNtuple *ntuple)
{
  // TNtuple *ntuple = new TNtuple("t", "SvxGeoHit variables",
  // "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid:event"

  assert(ntuple);
  int nj = 15;
  int j = 0;

  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      for (int ihit=0; ihit<trk.nhits; ihit++)
      {
        std::vector<float> vars(nj, 0.);
        SvxGeoHit hit = trk.GetHit(ihit);

        j = 0;
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
        vars[j++] = ev         ;
        ntuple->Fill(&vars[0]);
      }
    }

  return;
}

int
GetCorrections(const char *resFile, std::map<int, double> &mpc)
{
  // Load global parameter (i.e. position) corrections from pede output
  // (called millepede.res by default) into the mpc map.

  std::ifstream filein(resFile);
  int label;
  double p, col3, col4, col5;

  if (!filein)
  {
    Error("GetCorrections() in VtxAlign.C", "Problem opening %s", resFile);
    return -1;
  }
  else
    for (std::string line; std::getline(filein, line);)
    {
      if (filein >> label >> p >> col3 >> col4 >> col5)
        mpc[label] = p;
    }

  return (int)mpc.size();
}

void
WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfiles,
               double regFactor, double preSigma)
{
  cout << "Writing " << filename << "..." << flush;

  ofstream fs(filename);

  fs << Form("* This is %s", filename) << endl;
  fs << Form("* Pass this file to pede: pede %s", filename) << endl;
  fs << endl;

  fs << "Fortranfiles  ! Fortran/text inputs listed here:" << endl;
  for (unsigned int i=0; i<constfiles.size(); i++)
    fs << constfiles[i] << " ! constraints file" << endl;
  fs << endl;

  fs << "Cfiles  ! c/c++ binary input files listed here:" << endl;
  for (unsigned int i=0; i<binfiles.size(); i++)
    fs << binfiles[i] << " ! binary data file" << endl;
  fs << endl;

  fs << Form("regularisation %.3g %.3g ", regFactor, preSigma);
  fs << "! regularisation factor, pre-sigma" << endl;
  fs << endl;

  fs << "method inversion 5 0.0001  ! Gauss. elim., #iterations, tol." << endl;
  fs << "end" << endl;

  fs.close();

  Printf("done.");

  return;
}


#endif
