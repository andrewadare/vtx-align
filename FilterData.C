#include "VtxAlignBase.h"
#include "DataRejector.h"

#include <map>

// This "preprocessing" script does some initial outlier rejection.
// It also renames the ntuple(s).
// This script uses the PISA geometry file listed in the provided config file. 
// The opt string may contain:
//      "": vtx only.
//      "cnt": vtxtrks & cnttrks.
//      "no-refit": don't refit the tracks after filtering
//      "fixed-bc": use beamcenter from rootfiles/bc-411768.root when refitting
//                  tracks after filtering

void FilterData(const char *infilename = "rootfiles/anavtxcluster_411768-0-3.root",
                const char *outfilename = "rootfiles/411768-1-0.root",
                const char *configfilename = "production/config/config-411768-1-0.txt",
                double vertexprobmin = 0.02,
                double vertexprobmax = 0.98,
                double maxdca = 0.5,
                double maxres_s = 0.1,
                double maxres_z = 0.1,
                int nhitsmin = 3,
                int nevents = -1, // -1 = everything
                float frac4hit = -1, // -1 = no filter
                TString opt = "fixed-bc") // see above
{
  std::cout << "-- Opening " << infilename << " --" << std::endl;
  TFile *inFile = new TFile(infilename, "read");

  std::cout << "\n-- Reading seg_clusntuple --" << std::endl;
  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  TNtuple *svxcnt;
  assert(svxseg);

  std::cout << "\n-- Initializing SvxTGeo --" << std::endl;
  float bc[2] = {0};
  float e2w[3] = {0};
  float v2c[3] = {0};
  string parFileName;
  GetParamFromConfig(configfilename, bc, e2w, v2c, parFileName);
  SvxTGeo *tgeo = VTXModel(parFileName.c_str());

  std::cout << "\n-- Creating output file " << outfilename << " --" << std::endl;
  TFile *outFile = new TFile(outfilename, "recreate");

  std::cout << "\n-- Creating vtxtrks Tree --" << std::endl;
  TTree *vtxtrks = CreateTree("vtxtrks");
  TTree *cnttrks;

  geoEvents vtxevents;
  geoEvents cntevents;

  if (opt.Contains("cnt"))
  {
    std::cout << "\n-- Reading cnt_clusntuple --" << std::endl;
    svxcnt = (TNtuple *) inFile->Get("cnt_clusntuple");
    assert(svxcnt);

    std::cout << "\n-- Creating cnttrks Tree --" << std::endl;
    cnttrks = CreateTree("cnttrks");

    GetSyncedEventsFromClusTree(svxseg, svxcnt,
                                vtxevents, cntevents,
                                tgeo, nevents);
  }
  else
  {
    std::cout << "\n-- Reading events from ntuple --" << std::endl;
    map <int, int> tmpmap;
    GetEventsFromClusTree(svxseg, tgeo, vtxevents, tmpmap, nevents);
  }

  //need to fit segments to get phi, theta
  //NOTE: For the purposes of rejecting outlider, the
  //      first fit should always find the primary vertex
  //      and calculate dca with respect to that.
  std::cout << "\n-- Fitting tracks --" << std::endl;
  FitTracks(vtxevents, 0, "find_vertex, calc_dca");

  std::cout << "\n-- Filtering data --" << std::endl;
  geoEvents accvtxevents;
  geoEvents acccntevents;
  FilterData(vtxevents, cntevents,
             accvtxevents, acccntevents,
             vertexprobmin, vertexprobmax,
             maxdca,
             maxres_s, maxres_z,
             10,
             10,
             nhitsmin,
             frac4hit);

  //refit tracks with fixed bc if desired
  TGraphErrors *gbc = NULL;
  if (opt.Contains("no-refit"))
  {
    //don't refit
  }
  else if (opt.Contains("fixed-bc"))
  {
    //use fixed bc from designated file
    TString bcFileIn   = "rootfiles/bc-411768.root";
    TFile *bcf         = new TFile(bcFileIn.Data(), "read");
    gbc  = (TGraphErrors *) bcf->Get("gbc");
  }
  else
  {
    //use bc from config file
    gbc = new TGraphErrors();
    gbc->SetTitle("Beamcenter from east arm (point 0) and west arm (point 1)");
    gbc->SetPoint(0, bc[0] + e2w[0], bc[1] + e2w[1]);
    gbc->SetPoint(1, bc[0], bc[1]);
    gbc->SetPointError(0, 0.015, 0.010);
    gbc->SetPointError(1, 0.015, 0.010);
  }

  if (gbc)
  {
    Printf("\n-- Refitting tracks with desired bc --");
    Info("", " E: (%.3f, %.3f)", gbc->GetX()[0], gbc->GetY()[0]);
    Info("", " W: (%.3f, %.3f)", gbc->GetX()[1], gbc->GetY()[1]);

    // Use gbc in track fits. Find primary vertex.
    FitTracks(accvtxevents, gbc, "fit_to_bc, find_vertex");

    // Refit using z vertex from previous fit.
    // Finally, compute DCA wrt gbc and z vertex.
    FitTracks(accvtxevents, gbc, "fit_to_bc, fit_to_z_vertex, calc_dca");
  }

  Printf("\n-- Writing output to %s --", outfilename);
  FillTree(accvtxevents, vtxtrks);
  if (cntevents.size() > 0)
    FillTree(acccntevents, cnttrks);

  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);

  // For now, just write the cnt ntuple to file without filtering.
  if (opt.Contains("cnt"))
  {
    cnttrks->Write();
  }

  return;
}
