#ifndef __VTXIO_H__
#define __VTXIO_H__

#include "VtxAlignBase.h"
#include <TLeaf.h>
#include <TTree.h>

#include <utility>
#include <map>

using namespace std;

void GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &trks, int nmax = -1);
void GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &trks, multimap<int, int> &eventsync, int nmax = -1,
                       TString opt = "");
void GetEventsFromClusTree(TNtuple *t, SvxTGeo *geo, geoEvents &evts, map<int, int> &eventsync, int nmax = -1,
                           TString opt = "");
void GetSyncedEventsFromClusTree(TNtuple *tseg, TNtuple *tcnt,//<-- input
                                 geoEvents &segevts, geoEvents &cntevts,//<-- output
                                 SvxTGeo, int nmax = -1);
TTree *CreateTree(const char *name);
//void FillTree(SvxGeoTrack &gt, TTree *t, int evt = -1);
//void FillTree(geoTracks &trks, TTree *t, int evt = -1);
void FillTree(geoEvents &events, TTree *t);
void GetEventsFromTree(TTree *t, SvxTGeo *geo, geoEvents &evts, int nmax = -1,
                       TString opt = "");
void FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple, int event = -1);
void FillNTuple(geoEvents &events, TNtuple *ntuple);
int GetCorrections(const char *resFile, std::map<int, double> &mpc);
void WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfiles,
                    double regFactor = 1.0, double preSigma = 0.01);
vecs SplitLine(const string &line, const char *delim = " ");

int GetOffsetsFromConfig(const char *configFile, float e2w[], float v2c[]);
string GetGeomFromConfig(const char *configFile);

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

  for (int i = 0; i < nentries; i++)
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
GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks,
                  std::multimap<int, int> &eventsync, int nmax, TString /*opt*/)
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

  for (int i = 0; i < nentries; i++)
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
      eventsync.insert(pair<int, int>(itrack, event));
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
GetEventsFromClusTree(TNtuple *t, SvxTGeo *geo, geoEvents &events, map<int, int> &eventsync,
                      int nmax, TString opt)
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
  for (int i = 0; i < nentries; i++)
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
      //add entry to multimap
      //-> Key = original event index
      //-> Element = geoEvents vector index
      if (eventsync.find(ev) != eventsync.end())
      {
        cout << "ERROR!! event index=" << ev
             << " is already in map. Aborting"
             << endl;
        return;
      }
      else
      {
        eventsync.insert(pair<int, int>(ev, events.size() - 1));
      }
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


void GetSyncedEventsFromClusTree(TNtuple *tseg, TNtuple *tcnt,//<-- input
                                 geoEvents &segevts, geoEvents &cntevts,//<-- output
                                 SvxTGeo *geo, int nmax)
{
  //This option reads in entries from segment and svxcnt cluster
  //ntuples, matches event indecies between them
  //and output's them to geoEvents object, where the same
  //index in both geoEvents corresponds to the same
  //event


  //need a pair of temporary geoEvent objects for reading complete entries
  geoEvents tmpsegevts;
  geoEvents tmpcntevts;

  //maps for syncing events
  map<int, int> segmap;
  map<int, int> cntmap;

  //read the trees
  GetEventsFromClusTree(tseg, geo, tmpsegevts, segmap, nmax);
  GetEventsFromClusTree(tcnt, geo, tmpcntevts, cntmap, nmax, "cnt");

  //Loop over cnt events (more restrictive) and find corresponding
  //seg events, write results to output geoEvents
  map<int, int>::iterator cntpos;
  map<int, int>::iterator segpos;
  for (cntpos = cntmap.begin(); cntpos != cntmap.end(); ++cntpos)
  {
    //find the corresponding seg index
    segpos = segmap.find(cntpos->first);
    if (segpos != segmap.end())
    {
      cntevts.push_back(tmpcntevts[cntpos->second]);
      segevts.push_back(tmpsegevts[segpos->second]);
    }
  }

  return;

}


TTree *
CreateTree(const char *name)
{
  // Create A TTree to hold track information.
  // For use within alignment code.
  // Each entry is track based, with arrays for
  //  each cluster.

  // Declare the tree branch variables
  // is there a way to get around this,
  // since they aren't being used here?
  int event;
  int NTracks;
  int trkid;
  float trkphi;
  float trktheta;
  int NClus;
  int layer[8];
  int ladder[8];
  int sensor[8];
  float lx[8];
  float ly[8];
  float lz[8];
  float gx[8];
  float gy[8];
  float gz[8];
  float x_size[8];
  float z_size[8];
  float res_z[8];
  float res_s[8];

  TTree *t = new TTree(name, "Track based TTree for use in VTX Alignment");
  t->Branch("event", &event, "event/I");
  t->Branch("NTracks", &NTracks, "NTracks/I");
  t->Branch("trkid", &trkid, "trkid/I");
  t->Branch("trkphi", &trkphi, "trkphi/F");
  t->Branch("trktheta", &trktheta, "trktheta/F");
  t->Branch("NClus", &NClus, "NClus/I");
  t->Branch("layer", &layer, "layer[NClus]/I");
  t->Branch("ladder", &ladder, "ladder[NClus]/I");
  t->Branch("sensor", &sensor, "sensor[NClus]/I");
  t->Branch("lx", &lx, "lx[NClus]/F");
  t->Branch("ly", &ly, "ly[NClus]/F");
  t->Branch("lz", &lz, "lz[NClus]/F");
  t->Branch("gx", &gx, "gx[NClus]/F");
  t->Branch("gy", &gy, "gy[NClus]/F");
  t->Branch("gz", &gz, "gz[NClus]/F");
  t->Branch("x_size", &x_size, "x_size[NClus]/F");
  t->Branch("z_size", &z_size, "z_size[NClus]/F");
  t->Branch("res_z", &res_z, "res_z[NClus]/F");
  t->Branch("res_s", &res_s, "res_s[NClus]/F");


  return t;

}
/*
void
FillTree(SvxGeoTrack &gt, TTree *t, int evt)
{
  // Function to fill the TTree with
  // information from SvxGeoTrack object.
  // See CreateTree() for tree structure.
  assert(t);

  // Declare the tree branch variables
  int event;
  int NTracks;
  int trkid;
  float trkphi;
  float trktheta;
  int NClus;
  int layer[8];
  int ladder[8];
  int sensor[8];
  float lx[8];
  float ly[8];
  float lz[8];
  float gx[8];
  float gy[8];
  float gz[8];
  float x_size[8];
  float z_size[8];
  float res_z[8];
  float res_s[8];

  // Set the branch addresses
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("NTracks", &NTracks);
  t->SetBranchAddress("trkid", &trkid);
  t->SetBranchAddress("trkphi", &trkphi);
  t->SetBranchAddress("trktheta", &trktheta);
  t->SetBranchAddress("NClus", &NClus);
  t->SetBranchAddress("layer", &layer);
  t->SetBranchAddress("ladder", &ladder);
  t->SetBranchAddress("sensor", &sensor);
  t->SetBranchAddress("lx", &lx);
  t->SetBranchAddress("ly", &ly);
  t->SetBranchAddress("lz", &lz);
  t->SetBranchAddress("gx", &gx);
  t->SetBranchAddress("gy", &gy);
  t->SetBranchAddress("gz", &gz);
  t->SetBranchAddress("x_size", &x_size);
  t->SetBranchAddress("z_size", &z_size);
  t->SetBranchAddress("res_z", &res_z);
  t->SetBranchAddress("res_s", &res_s);

  event = evt;
  NTracks = -1; //if filling w/ this method, NTracks is unknown
  trkphi = gt.phi0;
  trktheta = gt.the0;
  NClus = gt.nhits;
  for (int ihit = 0; ihit < gt.nhits; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    layer[ihit] = hit.layer;
    ladder[ihit] = hit.ladder;
    sensor[ihit] = hit.sensor;
    lx[ihit] = hit.xs;
    ly[ihit] = hit.ys;
    lz[ihit] = hit.zs;
    gx[ihit] = hit.x;
    gy[ihit] = hit.y;
    gz[ihit] = hit.z;
    x_size[ihit] = hit.xsigma;
    z_size[ihit] = hit.zsigma;
    res_z[ihit] = hit.dz;
    res_s[ihit] = hit.ds;
    trkid = hit.trkid;
  }

  t->Fill();

  return;
}

void
FillTree(geoTracks &trks, TTree *t, int evt)
{
  // Function to fill the TTree with
  // information from SvxGeoTrack object.
  // See CreateTree() for tree structure.
  assert(t);

  // Declare the tree branch variables
  int event;
  int NTracks;
  int trkid;
  float trkphi;
  float trktheta;
  int NClus;
  int layer[8];
  int ladder[8];
  int sensor[8];
  float lx[8];
  float ly[8];
  float lz[8];
  float gx[8];
  float gy[8];
  float gz[8];
  float x_size[8];
  float z_size[8];
  float res_z[8];
  float res_s[8];

  // Set the branch addresses
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("NTracks", &NTracks);
  t->SetBranchAddress("trkid", &trkid);
  t->SetBranchAddress("trkphi", &trkphi);
  t->SetBranchAddress("trktheta", &trktheta);
  t->SetBranchAddress("NClus", &NClus);
  t->SetBranchAddress("layer", &layer);
  t->SetBranchAddress("ladder", &ladder);
  t->SetBranchAddress("sensor", &sensor);
  t->SetBranchAddress("lx", &lx);
  t->SetBranchAddress("ly", &ly);
  t->SetBranchAddress("lz", &lz);
  t->SetBranchAddress("gx", &gx);
  t->SetBranchAddress("gy", &gy);
  t->SetBranchAddress("gz", &gz);
  t->SetBranchAddress("x_size", &x_size);
  t->SetBranchAddress("z_size", &z_size);
  t->SetBranchAddress("res_z", &res_z);
  t->SetBranchAddress("res_s", &res_s);

  for (unsigned int itrk = 0; itrk < trks.size(); itrk++)
  {
    SvxGeoTrack gt = trks.at(itrk);
    event = evt;
    NTracks = trks.size(); //if filling w/ this method, NTracks is unknown
    trkphi = gt.phi0;
    trktheta = gt.the0;
    NClus = gt.nhits;
    for (int ihit = 0; ihit < gt.nhits; ihit++)
    {
      SvxGeoHit hit = gt.GetHit(ihit);
      layer[ihit] = hit.layer;
      ladder[ihit] = hit.ladder;
      sensor[ihit] = hit.sensor;
      lx[ihit] = hit.xs;
      ly[ihit] = hit.ys;
      lz[ihit] = hit.zs;
      gx[ihit] = hit.x;
      gy[ihit] = hit.y;
      gz[ihit] = hit.z;
      x_size[ihit] = hit.xsigma;
      z_size[ihit] = hit.zsigma;
      res_z[ihit] = hit.dz;
      res_s[ihit] = hit.ds;
      trkid = hit.trkid;
    }
    t->Fill();
  }

  return;
}
*/
void
FillTree(geoEvents &events, TTree *t)
{
  // Function to fill the TTree with
  // information from SvxGeoTrack object.
  // See CreateTree() for tree structure.
  assert(t);

  // Declare the tree branch variables
  int event;
  int NTracks;
  int trkid;
  float trkphi;
  float trktheta;
  int NClus;
  int layer[8];
  int ladder[8];
  int sensor[8];
  float lx[8];
  float ly[8];
  float lz[8];
  float gx[8];
  float gy[8];
  float gz[8];
  float x_size[8];
  float z_size[8];
  float res_z[8];
  float res_s[8];

  // Set the branch addresses
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("NTracks", &NTracks);
  t->SetBranchAddress("trkid", &trkid);
  t->SetBranchAddress("trkphi", &trkphi);
  t->SetBranchAddress("trktheta", &trktheta);
  t->SetBranchAddress("NClus", &NClus);
  t->SetBranchAddress("layer", &layer);
  t->SetBranchAddress("ladder", &ladder);
  t->SetBranchAddress("sensor", &sensor);
  t->SetBranchAddress("lx", &lx);
  t->SetBranchAddress("ly", &ly);
  t->SetBranchAddress("lz", &lz);
  t->SetBranchAddress("gx", &gx);
  t->SetBranchAddress("gy", &gy);
  t->SetBranchAddress("gz", &gz);
  t->SetBranchAddress("x_size", &x_size);
  t->SetBranchAddress("z_size", &z_size);
  t->SetBranchAddress("res_z", &res_z);
  t->SetBranchAddress("res_s", &res_s);

  for (unsigned int ev = 0; ev < events.size(); ev++)
  {
    for (unsigned int itrk = 0; itrk < events[ev].size(); itrk++)
    {

      SvxGeoTrack trk = events[ev][itrk];

      event = (int)ev;
      NTracks = (int)events[ev].size();
      trkid = (int)itrk;
      trkphi = trk.phi0;
      trktheta = trk.the0;
      NClus = trk.nhits;

      for (int ihit = 0; ihit < trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);

        layer[ihit] = hit.layer;
        ladder[ihit] = hit.ladder;
        sensor[ihit] = hit.sensor;
        lx[ihit] = hit.xs;
        ly[ihit] = hit.ys;
        lz[ihit] = hit.zs;
        gx[ihit] = hit.x;
        gy[ihit] = hit.y;
        gz[ihit] = hit.z;
        x_size[ihit] = hit.xsigma;
        z_size[ihit] = hit.zsigma;
        res_z[ihit] = hit.dz;
        res_s[ihit] = hit.ds;
      }

      t->Fill();
    }
  }

  return;
}

void
GetEventsFromTree(TTree *t, SvxTGeo *geo, geoEvents &evts, int nmax,
                  TString /*opt*/)
{
  // Get events from TTree and put information in geoEvents
  assert(t);
  assert(geo);

  // Declare the tree branch variables
  int event;
  int NTracks;
  int trkid;
  float trkphi;
  float trktheta;
  int NClus;
  int layer[8];
  int ladder[8];
  int sensor[8];
  float lx[8];
  float ly[8];
  float lz[8];
  float gx[8];
  float gy[8];
  float gz[8];
  float x_size[8];
  float z_size[8];
  float res_z[8];
  float res_s[8];

  // Set the branch addresses
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("NTracks", &NTracks);
  t->SetBranchAddress("trkid", &trkid);
  t->SetBranchAddress("trkphi", &trkphi);
  t->SetBranchAddress("trktheta", &trktheta);
  t->SetBranchAddress("NClus", &NClus);
  t->SetBranchAddress("layer", &layer);
  t->SetBranchAddress("ladder", &ladder);
  t->SetBranchAddress("sensor", &sensor);
  t->SetBranchAddress("lx", &lx);
  t->SetBranchAddress("ly", &ly);
  t->SetBranchAddress("lz", &lz);
  t->SetBranchAddress("gx", &gx);
  t->SetBranchAddress("gy", &gy);
  t->SetBranchAddress("gz", &gz);
  t->SetBranchAddress("x_size", &x_size);
  t->SetBranchAddress("z_size", &z_size);
  t->SetBranchAddress("res_z", &res_z);
  t->SetBranchAddress("res_s", &res_s);


  long nentries = t->GetEntries();
  int prevev  = -1;
  int ntracks = 0;
  int nhits   = 0;


  Printf("Reading events from Tree (%d available hits)...", (int)nentries);
  for (int i = 0; i < nentries; i++)
  {
    t->GetEntry(i);

    if (event != prevev) // New event
    {
      if (nmax > 0 && (int)evts.size() == nmax)
      {
        Printf("%d events imported (%d tracks, %d hits).",
               (int)evts.size(), ntracks, nhits);
        return;
      }

      geoTracks tracks;
      //if we know the number of tracks in the event
      //reserve the space to avoid reallocating the vector
      if (NTracks >= 0)
        tracks.reserve(NTracks);

      evts.push_back(tracks);
    }

    SvxGeoTrack trk;
    trk.phi0 = trkphi;
    trk.the0 = trktheta;
    trk.nhits = NClus;

    //save some reallocation by knowing the number of clusters
    trk.hits.reserve(NClus);

    for (int j = 0; j < NClus; j++)
    {
      SvxGeoHit hit;
      hit.layer  = layer[j];
      hit.ladder = ladder[j];
      hit.sensor = sensor[j];
      hit.xs     = lx[j];
      hit.ys     = ly[j];
      hit.zs     = lz[j];
      hit.x      = gx[j];
      hit.y      = gy[j];
      hit.z      = gz[j];
      hit.xsigma = x_size[j]; // N.B. this is not sigma!
      hit.zsigma = z_size[j]; // N.B. this is not sigma!
      hit.dz     = res_z[j];
      hit.ds     = res_s[j];
      hit.trkid  = trkid;
      hit.node   = geo->SensorNode(hit.layer, hit.ladder, hit.sensor);

      // Overwrite x,y,z with local-->global transformation on xs,ys,zs
      hit.Update();

      trk.hits.push_back(hit);
      nhits++;
    }

    evts.back().push_back(trk);
    ntracks++;

    prevev = event;
  }

  Printf("%d events imported (%d tracks, %d hits).",
         (int)evts.size(), ntracks, nhits);

  return;
}




void
FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple, int event)
{
  // TNtuple *ntuple = new TNtuple("t", "SvxGeoHit variables",
  // "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:"
  // "res_z:res_s:trkid:event"

  assert(ntuple);

  std::cout << std::endl;
  std::cout << "WARNING!! This function is deprecated and "
            << "should no longer be used! "
            << "Please update to FillTree()."
            << std::endl;

  if (false)
    Printf("ds (%8.4f,%8.4f,%8.4f,%8.4f),  "
           "dz (%8.4f,%8.4f,%8.4f,%8.4f)",
           gt.hits[0].ds, gt.hits[1].ds, gt.hits[2].ds, gt.hits[3].ds,
           gt.hits[0].dz, gt.hits[1].dz, gt.hits[2].dz, gt.hits[3].dz);

  for (int ihit = 0; ihit < gt.nhits; ihit++)
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

  std::cout << std::endl;
  std::cout << "WARNING!! This function is deprecated and "
            << "should no longer be used! "
            << "Please update to FillTree()."
            << std::endl;

  for (unsigned int ev = 0; ev < events.size(); ev++)
    for (unsigned int t = 0; t < events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      for (int ihit = 0; ihit < trk.nhits; ihit++)
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
    Error("GetCorrections() in VtxIO.h", "Problem opening %s", resFile);
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
  for (unsigned int i = 0; i < constfiles.size(); i++)
    fs << constfiles[i] << " ! constraints file" << endl;
  fs << endl;

  fs << "Cfiles  ! c/c++ binary input files listed here:" << endl;
  for (unsigned int i = 0; i < binfiles.size(); i++)
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

int
GetOffsetsFromConfig(const char *configFile, float e2w[], float v2c[])
{
  // Read east-to-west and vtx-to-cnt offsets from configFile
  assert(configFile);
  ifstream fin;
  fin.open(configFile);

  if (!fin)
  {
    Error("GetOffsetsFromConfig() in VtxIO.h", "Problem opening %s",
          configFile);
    return -1;
  }

  Info("", "Reading offsets from %s...", configFile);

  for (string line; getline(fin, line);)
  {
    if (line.find("east-to-west") != string::npos)
    {
      vecs words = SplitLine(line);
      assert(words.size() == 4);
      for (int i = 0; i < 3; ++i)
        e2w[i] = atof(words[i + 1].c_str());

      Printf("E-W offset: %f %f %f", e2w[0], e2w[1], e2w[2]);
    }
    if (line.find("vtx-to-cnt") != string::npos)
    {
      vecs words = SplitLine(line);
      assert(words.size() == 4);
      for (int i = 0; i < 3; ++i)
        v2c[i] = atof(words[i + 1].c_str());

      Printf("VTX-CNT offset: %f %f %f", v2c[0], v2c[1], v2c[2]);
    }
  }

  return 0;
}

string
GetGeomFromConfig(const char *configFile)
{
  // Get the geometry file name from a config file
  assert(configFile);
  ifstream fin;
  fin.open(configFile);

  if (!fin)
  {
    Error("GetGeomFromConfig() in VtxIO.h", "Problem opening %s",
          configFile);
    return "";
  }

  Info("", "Reading geom file from %s...", configFile);

  for (string line; getline(fin, line);)
  {
    if (line.find("geomfile") != string::npos)
    {
      vecs words = SplitLine(line);
      assert(words.size() == 2);
      return words[1];
    }
  }

  Error("GetGeomFromConfig() in VtxIO.h",
        "Unable to find geomfile in %s", configFile);

  return "";
}

vecs
SplitLine(const string &line, const char *delim)
{
  vecs result;

  size_t prev = 0;
  size_t next = 0;

  while ((next = line.find_first_of(delim, prev)) != string::npos)
  {
    if (next - prev != 0)
      result.push_back(line.substr(prev, next - prev));
    prev = next + 1;
  }

  if (prev < line.size())
    result.push_back(line.substr(prev));

  return result;
}

#endif
