#include "DataRejector.h"
#include <TFile.h>

typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

void FilterData(const char *infilename = "rootfiles/411768_july17_ideal.root",
                const char *outfilename = "rootfiles/test.root")
{
  TFile *inFile = new TFile(infilename, "read");
  TFile *outFile = new TFile(outfilename, "recreate");

  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  if (!svxseg)
  {
    svxseg = (TNtuple *)inFile->Get("ht2");
    svxseg->SetName("seg_clusntuple");
  }

  assert(svxseg);

  const char *hitvars = "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:"
                        "x_size:z_size:res_z:res_s:trkid:event";
  TNtuple *ht1 = new TNtuple("ht1", "SvxGeoHit variables", hitvars);

  SvxTGeo *tgeo = new SvxTGeo;
  tgeo->ReadParFile("geom/svxPISA-ideal.par");
  tgeo->MakeTopVolume(100, 100, 100);
  tgeo->AddSensors();
  TGeoManager *mgr = tgeo->GeoManager();
  mgr->CloseGeometry();

  geoEvents vtxevents;
  GetEventsFromTree(svxseg, tgeo, vtxevents, -1);
  FitTracks(vtxevents);
  FilterData(vtxevents, 0.01, 0.99, 0.2, 0.1, 0.1, ht1);

  Printf("%.1f%% of initial data rejected.", 
         100*(1.0 - (float)ht1->GetEntries()/svxseg->GetEntries()));

  Printf("Writing output to %s", outfilename);
  outFile->cd();
  ht1->Write("seg_clusntuple");

  return;
}