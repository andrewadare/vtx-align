#include "VtxAlignBase.h"
#include "DataRejector.h"

// This "preprocessing" script does some initial outlier rejection.
// It also renames the ntuple(s).
// It is conventional for outfilename and pisafilename to have the same base.
// For record-keeping transparency, symlink the needed geometry file to 
// run-prod-iter.par so the base matches the ROOT file.
// The alignment script VtxAlign.C requires a .par file with such a name.

void FilterData(const char *infilename = "rootfiles/anavtxcluster_411768-pro8.root",
                const char *outfilename = "rootfiles/411768-8-0.root",
                const char *pisafilename = "geom/411768-8-0.par",
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
