#include "DataRejector.h"
#include <TFile.h>

typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

void FilterData(const char *infilename = "rootfiles/411768_july17_ideal.root",
                const char *outfilename = "rootfiles/411768_july17_ideal_filt.root",
                double vertexprobmin = 0.01,
                double vertexprobmax = 0.99,
                double maxdca = 0.5,
                double maxres_s = 0.1,
                double maxres_z = 0.1)
{
  TFile *inFile = new TFile(infilename, "read");

  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  assert(svxseg);
    
  TFile *outFile = new TFile(outfilename, "recreate");
  
  const char *hitvars = "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:"
                        "x_size:z_size:res_z:res_s:trkid:event";
  TNtuple *ht0 = new TNtuple("ht0", "SvxGeoHit variables", hitvars);

  SvxTGeo *tgeo = new SvxTGeo;
  tgeo->ReadParFile("geom/svxPISA-ideal.par");
  tgeo->MakeTopVolume(100, 100, 100);
  tgeo->AddSensors();
  TGeoManager *mgr = tgeo->GeoManager();
  mgr->CloseGeometry();

  geoEvents vtxevents;
  GetEventsFromTree(svxseg, tgeo, vtxevents, -1);
  FitTracks(vtxevents);
  FilterData(vtxevents,
             vertexprobmin,
             vertexprobmax,
             maxdca,
             maxres_s,
             maxres_z,
             ht0);

  Printf("%.1f%% of initial data rejected.",
         100*(1.0 - (float)ht0->GetEntries()/svxseg->GetEntries()));

  Printf("Writing output to %s", outfilename);
  outFile->cd();
  ht0->Write(0, TObject::kOverwrite);

  return;
}