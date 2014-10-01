#ifndef __VTXIO_H__
#define __VTXIO_H__

#include "VtxAlignBase.h"
#include <TLeaf.h>
#include <TTree.h>

#include <utility>
#include <map>

using namespace std;

void GetEventsFromClusTree(TNtuple *t, SvxTGeo *geo, geoEvents &evts,
                           map<int, int> &eventsync, int nmax = -1,
                           TString opt = "");
void GetSyncedEventsFromClusTree(TNtuple *tseg, TNtuple *tcnt,//<-- input
                                 geoEvents &segevts, geoEvents &cntevts,//<-- output
                                 SvxTGeo, int nmax = -1);
TTree *CreateTree(const char *name);
void FillTree(geoEvents &events, TTree *t);
void GetEventsFromTree(TTree *t, SvxTGeo *geo, geoEvents &evts, int nmax = -1,
                       TString opt = "");
int GetCorrections(const char *resFile, std::map<int, double> &mpc);
void WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfiles,
                    double regFactor = 1.0, double preSigma = 0.01);
vecs SplitLine(const string &line, const char *delim = " ");

int GetOffsetsFromConfig(const char *configFile, float e2w[], float v2c[]);
string GetGeomFromConfig(const char *configFile);









void
GetEventsFromClusTree(TNtuple *t, SvxTGeo *geo, geoEvents &events,
                      map<int, int> &eventsync, int nmax, TString opt)
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

struct VtxTrackVars
{
  VtxTrackVars();
  ~VtxTrackVars() {};

  Int_t   event;
  Int_t   NTracks;
  Int_t   trkid;
  Int_t   NClus;
  Int_t   layer[8];
  Int_t   ladder[8];
  Int_t   sensor[8];
  Float_t trkphi;
  Float_t trktheta;
  Float_t primvtx[3];
  Float_t xydca;
  Float_t zdca;
  Float_t phirot;
  Float_t yp0;
  Float_t z0;
  Float_t lx[8];
  Float_t ly[8];
  Float_t lz[8];
  Float_t gx[8];
  Float_t gy[8];
  Float_t gz[8];
  Float_t x_size[8];
  Float_t z_size[8];
  Float_t res_z[8];
  Float_t res_s[8];
};

VtxTrackVars::VtxTrackVars() :
  event(0),
  NTracks(0),
  trkid(0),
  NClus(0),
  trkphi(0.),
  trktheta(0.),
  xydca(0.),
  zdca(0.),
  phirot(0.),
  yp0(0.),
  z0(0.)
{

  for (int i = 0; i < 3; i++)
  {
    primvtx[i] = 0;
  }

  for (int i = 0; i < 8; ++i)
  {
    layer[i] = 0;
    ladder[i] = 0;
    sensor[i] = 0;
    lx[i] = 0;
    ly[i] = 0;
    lz[i] = 0;
    gx[i] = 0;
    gy[i] = 0;
    gz[i] = 0;
    x_size[i] = 0;
    z_size[i] = 0;
    res_z[i] = 0;
    res_s[i] = 0;
  }
}

TTree *
CreateTree(const char *name)
{
  // Create A TTree to hold track information.
  // For use within alignment code.
  // Each entry is track based, with arrays for
  //  each cluster.

  TTree *t = new TTree(name, "Track based TTree for use in VTX Alignment");
  VtxTrackVars vtv;

  t->Branch("event"    , &vtv.event   , "event/I");
  t->Branch("NTracks"  , &vtv.NTracks , "NTracks/I");
  t->Branch("trkid"    , &vtv.trkid   , "trkid/I");
  t->Branch("trkphi"   , &vtv.trkphi  , "trkphi/F");
  t->Branch("trktheta" , &vtv.trktheta, "trktheta/F");
  t->Branch("primvtx"  , &vtv.primvtx , "primvtx[3]/F");
  t->Branch("xydca"    , &vtv.xydca   , "xydca/F");
  t->Branch("zdca"     , &vtv.zdca    , "zdca/F");
  t->Branch("phirot"   , &vtv.phirot  , "phirot/F");
  t->Branch("yp0"      , &vtv.yp0     , "yp0/F");
  t->Branch("z0"       , &vtv.z0      , "z0/F");
  t->Branch("NClus"    , &vtv.NClus   , "NClus/I");
  t->Branch("layer"    , &vtv.layer   , "layer[NClus]/I");
  t->Branch("ladder"   , &vtv.ladder  , "ladder[NClus]/I");
  t->Branch("sensor"   , &vtv.sensor  , "sensor[NClus]/I");
  t->Branch("lx"       , &vtv.lx      , "lx[NClus]/F");
  t->Branch("ly"       , &vtv.ly      , "ly[NClus]/F");
  t->Branch("lz"       , &vtv.lz      , "lz[NClus]/F");
  t->Branch("gx"       , &vtv.gx      , "gx[NClus]/F");
  t->Branch("gy"       , &vtv.gy      , "gy[NClus]/F");
  t->Branch("gz"       , &vtv.gz      , "gz[NClus]/F");
  t->Branch("x_size"   , &vtv.x_size  , "x_size[NClus]/F");
  t->Branch("z_size"   , &vtv.z_size  , "z_size[NClus]/F");
  t->Branch("res_z"    , &vtv.res_z   , "res_z[NClus]/F");
  t->Branch("res_s"    , &vtv.res_s   , "res_s[NClus]/F");

  return t;
}

void
FillTree(geoEvents &events, TTree *t)
{
  // Function to fill the TTree with
  // information from SvxGeoTrack object.
  // See CreateTree() for tree structure.
  assert(t);

  VtxTrackVars                    vtv;
  t->SetBranchAddress("event",    &vtv.event);
  t->SetBranchAddress("NTracks",  &vtv.NTracks);
  t->SetBranchAddress("trkid",    &vtv.trkid);
  t->SetBranchAddress("trkphi",   &vtv.trkphi);
  t->SetBranchAddress("trktheta", &vtv.trktheta);
  t->SetBranchAddress("primvtx",  &vtv.primvtx);
  t->SetBranchAddress("xydca",    &vtv.xydca);
  t->SetBranchAddress("zdca",     &vtv.zdca);
  t->SetBranchAddress("phirot",   &vtv.phirot);
  t->SetBranchAddress("yp0",      &vtv.yp0);
  t->SetBranchAddress("z0",       &vtv.z0);
  t->SetBranchAddress("NClus",    &vtv.NClus);
  t->SetBranchAddress("layer",    &vtv.layer);
  t->SetBranchAddress("ladder",   &vtv.ladder);
  t->SetBranchAddress("sensor",   &vtv.sensor);
  t->SetBranchAddress("lx",       &vtv.lx);
  t->SetBranchAddress("ly",       &vtv.ly);
  t->SetBranchAddress("lz",       &vtv.lz);
  t->SetBranchAddress("gx",       &vtv.gx);
  t->SetBranchAddress("gy",       &vtv.gy);
  t->SetBranchAddress("gz",       &vtv.gz);
  t->SetBranchAddress("x_size",   &vtv.x_size);
  t->SetBranchAddress("z_size",   &vtv.z_size);
  t->SetBranchAddress("res_z",    &vtv.res_z);
  t->SetBranchAddress("res_s",    &vtv.res_s);

  for (unsigned int ev = 0; ev < events.size(); ev++)
  {
    for (unsigned int itrk = 0; itrk < events[ev].size(); itrk++)
    {

      SvxGeoTrack trk = events[ev][itrk];

      vtv.event      = (int)ev;
      vtv.NTracks    = (int)events[ev].size();
      vtv.trkid      = (int)itrk;
      vtv.trkphi     = trk.phi0;
      vtv.trktheta   = trk.the0;
      vtv.primvtx[0] = trk.vx;
      vtv.primvtx[1] = trk.vy;
      vtv.primvtx[2] = trk.vz;
      vtv.xydca      = trk.xydca;
      vtv.zdca       = trk.zdca;
      vtv.phirot     = trk.phirot;
      vtv.yp0        = trk.yp0;
      vtv.z0         = trk.z0;
      vtv.NClus      = trk.nhits;

      for (int ihit = 0; ihit < trk.nhits; ihit++)
      {
        SvxGeoHit hit = trk.GetHit(ihit);

        vtv.layer[ihit] = hit.layer;
        vtv.ladder[ihit] = hit.ladder;
        vtv.sensor[ihit] = hit.sensor;
        vtv.lx[ihit] = hit.xs;
        vtv.ly[ihit] = hit.ys;
        vtv.lz[ihit] = hit.zs;
        vtv.gx[ihit] = hit.x;
        vtv.gy[ihit] = hit.y;
        vtv.gz[ihit] = hit.z;
        vtv.x_size[ihit] = hit.xsigma;
        vtv.z_size[ihit] = hit.zsigma;
        vtv.res_z[ihit] = hit.dz;
        vtv.res_s[ihit] = hit.ds;
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

  VtxTrackVars                    vtv;
  t->SetBranchAddress("event",    &vtv.event);
  t->SetBranchAddress("NTracks",  &vtv.NTracks);
  t->SetBranchAddress("trkid",    &vtv.trkid);
  t->SetBranchAddress("trkphi",   &vtv.trkphi);
  t->SetBranchAddress("trktheta", &vtv.trktheta);
  t->SetBranchAddress("primvtx",  &vtv.primvtx);
  t->SetBranchAddress("xydca",    &vtv.xydca);
  t->SetBranchAddress("zdca",     &vtv.zdca);
  t->SetBranchAddress("phirot",   &vtv.phirot);
  t->SetBranchAddress("yp0",      &vtv.yp0);
  t->SetBranchAddress("z0",       &vtv.z0);
  t->SetBranchAddress("NClus",    &vtv.NClus);
  t->SetBranchAddress("layer",    &vtv.layer);
  t->SetBranchAddress("ladder",   &vtv.ladder);
  t->SetBranchAddress("sensor",   &vtv.sensor);
  t->SetBranchAddress("lx",       &vtv.lx);
  t->SetBranchAddress("ly",       &vtv.ly);
  t->SetBranchAddress("lz",       &vtv.lz);
  t->SetBranchAddress("gx",       &vtv.gx);
  t->SetBranchAddress("gy",       &vtv.gy);
  t->SetBranchAddress("gz",       &vtv.gz);
  t->SetBranchAddress("x_size",   &vtv.x_size);
  t->SetBranchAddress("z_size",   &vtv.z_size);
  t->SetBranchAddress("res_z",    &vtv.res_z);
  t->SetBranchAddress("res_s",    &vtv.res_s);

  long nentries = t->GetEntries();
  int prevev  = -1;
  int ntracks = 0;
  int nhits   = 0;


  Printf("Reading events from Tree (%d available hits)...", (int)nentries);
  for (int i = 0; i < nentries; i++)
  {
    t->GetEntry(i);

    if (vtv.event != prevev) // New event
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
      if (vtv.NTracks >= 0)
        tracks.reserve(vtv.NTracks);

      evts.push_back(tracks);
    }

    SvxGeoTrack trk;
    trk.phi0   = vtv.trkphi;
    trk.the0   = vtv.trktheta;
    trk.vx     = vtv.primvtx[0];
    trk.vy     = vtv.primvtx[1];
    trk.vz     = vtv.primvtx[2];
    trk.xydca  = vtv.xydca;
    trk.zdca   = vtv.zdca;
    trk.phirot = vtv.phirot;
    trk.yp0    = vtv.yp0;
    trk.z0     = vtv.z0;
    trk.nhits  = vtv.NClus;

    //save some reallocation by knowing the number of clusters
    trk.hits.reserve(vtv.NClus);

    for (int j = 0; j < vtv.NClus; j++)
    {
      SvxGeoHit hit;
      hit.layer  = vtv.layer[j];
      hit.ladder = vtv.ladder[j];
      hit.sensor = vtv.sensor[j];
      hit.xs     = vtv.lx[j];
      hit.ys     = vtv.ly[j];
      hit.zs     = vtv.lz[j];
      hit.x      = vtv.gx[j];
      hit.y      = vtv.gy[j];
      hit.z      = vtv.gz[j];
      hit.xsigma = vtv.x_size[j]; // N.B. this is not sigma!
      hit.zsigma = vtv.z_size[j]; // N.B. this is not sigma!
      hit.dz     = vtv.res_z[j];
      hit.ds     = vtv.res_s[j];
      hit.trkid  = vtv.trkid;
      hit.node   = geo->SensorNode(hit.layer, hit.ladder, hit.sensor);

      // Overwrite x,y,z with local-->global transformation on xs,ys,zs
      hit.Update();

      trk.hits.push_back(hit);
      nhits++;
    }

    evts.back().push_back(trk);
    ntracks++;

    prevev = vtv.event;
  }

  Printf("%d events imported (%d tracks, %d hits).",
         (int)evts.size(), ntracks, nhits);

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
