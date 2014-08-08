#include "VtxAlignBase.h"
#include "DataRejector.h"

void FilterData(const char *infilename = "rootfiles/anavtxcluster_406541-pro1.root",
                const char *outfilename = "rootfiles/406541-1-99.root",
                const char *pisafilename = "geom/411768-7-2.par",
                double vertexprobmin = 0.02,
                double vertexprobmax = 0.98,
                double maxdca = 0.5,
                double maxres_s = 0.1,
                double maxres_z = 0.1)
{
  TFile *inFile = new TFile(infilename, "read");

  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  assert(svxseg);

  TFile *outFile = new TFile(outfilename, "recreate");
  TNtuple *vtxhits = new TNtuple("vtxhits", "VTX hit variables", HITVARS);
  // TNtuple *cnthits = new TNtuple("cnthits", "CNT hit variables", HITVARS);

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
             vtxhits);

  Printf("%.1f%% of initial data rejected.",
         100*(1.0 - (float)vtxhits->GetEntries()/svxseg->GetEntries()));

  Printf("Writing output to %s", outfilename);
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);

  return;
}
