#include "VtxAlignBase.h"
#include "DataRejector.h"


// typedef vector<SvxGeoTrack> geoTracks;
// typedef vector<geoTracks> geoEvents;

void FilterData(const char *infilename = "rootfiles/411768_july17_ideal.root",
                const char *outfilename = "rootfiles/411768-0-0.root",
                const char *pisafilename = "geom/svxPISA-ideal.par",
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
  TNtuple *ht = new TNtuple("vtxhits", "VTX hit variables", HITVARS);

  SvxTGeo *tgeo = VTXModel(pisafilename);

  geoEvents vtxevents;
  GetEventsFromTree(svxseg, tgeo, vtxevents, -1);
  FitTracks(vtxevents);
  FilterData(vtxevents,
             vertexprobmin,
             vertexprobmax,
             maxdca,
             maxres_s,
             maxres_z,
             ht);

  Printf("%.1f%% of initial data rejected.",
         100*(1.0 - (float)ht->GetEntries()/svxseg->GetEntries()));

  Printf("Writing output to %s", outfilename);
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);

  return;
}