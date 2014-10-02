#include "VtxAlignBase.h"
#include "DataRejector.h"

#include <map>

// This "preprocessing" script does some initial outlier rejection.
// It also renames the ntuple(s).
// It is conventional for outfilename and pisafilename to have the same base.
// For record-keeping transparency, symlink the needed geometry file to
// run-prod-iter.par so the base matches the ROOT file.
// The alignment script VtxAlign.C requires a .par file with such a name.

void FilterData(const char *infilename = "rootfiles/anavtxcluster_411768-pro0-no-vtx2cnt.root",
                const char *outfilename = "rootfiles/411768-0-0.root",
                const char *pisafilename = "geom/411768-0-0.par",
                double vertexprobmin = 0.02,
                double vertexprobmax = 0.98,
                double maxdca = 0.5,
                double maxres_s = 0.1,
                double maxres_z = 0.1,
                int nhitsmin = 3,
                int nevents = -1, // -1 = everything
                float frac4hit = -1, // -1 = no filter
                TString opt = "cnt") // "": vtx only. "cnt": vtxtrks & cnttrks.
{
  std::cout << "-- Opening " << infilename << " --" << std::endl;
  TFile *inFile = new TFile(infilename, "read");

  std::cout << "-- Reading seg_clusntuple --" << std::endl;
  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  TNtuple *svxcnt;
  assert(svxseg);

  std::cout << "-- Initializing SvxTGeo --" << std::endl;
  SvxTGeo *tgeo = VTXModel(pisafilename);

  std::cout << "-- Creating output file " << outfilename << " --" << std::endl;
  TFile *outFile = new TFile(outfilename, "recreate");
  // TNtuple *vtxhits = new TNtuple("vtxhits", "VTX hit variables", HITVARS);
  // TNtuple *cnthits = new TNtuple("cnthits", "CNT hit variables", HITVARS);

  std::cout << "-- Creating vtxtrks Tree --" << std::endl;
  TTree *vtxtrks = CreateTree("vtxtrks");
  TTree *cnttrks;

  geoEvents vtxevents;
  geoEvents cntevents;

  if (opt.Contains("cnt"))
  {
    std::cout << "-- Reading cnt_clusntuple --" << std::endl;
    svxcnt = (TNtuple *) inFile->Get("cnt_clusntuple");
    assert(svxcnt);

    std::cout << "-- Creating vtxtrks Tree --" << std::endl;
    cnttrks = CreateTree("cnttrks");

    GetSyncedEventsFromClusTree(svxseg, svxcnt,
                                vtxevents, cntevents,
                                tgeo, nevents);

    //need to fit segments to get phi, theta
    std::cout << "-- Fitting tracks --" << std::endl;
    FitTracks(vtxevents);

    std::cout << "-- Filtering data --" << std::endl;
    FilterData(vtxevents, cntevents,
               vtxtrks, cnttrks,
               vertexprobmin, vertexprobmax,
               maxdca,
               maxres_s, maxres_z, 
               10, 
               10, 
               nhitsmin,
               frac4hit);
  }
  else
  {
    std::cout << "-- Reading events from ntuple --" << std::endl;
    map <int, int> tmpmap;
    GetEventsFromClusTree(svxseg, tgeo, vtxevents, tmpmap, nevents);

    std::cout << "-- Fitting tracks --" << std::endl;
    FitTracks(vtxevents);

    std::cout << "-- Filtering data --" << std::endl;
    FilterData(vtxevents,
               vtxtrks,
               vertexprobmin,
               vertexprobmax,
               maxdca,
               maxres_s,
               maxres_z,
               10,
               10,
               nhitsmin,
               frac4hit);
  }

  Printf("-- Writing output to %s --", outfilename);
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);

  // For now, just write the cnt ntuple to file without filtering.
  if (opt.Contains("cnt"))
  {
    cnttrks->Write();
  }

  return;
}
