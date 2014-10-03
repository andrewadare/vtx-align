#include "VtxAlignBase.h"
#include "DataRejector.h"

#include <map>

// This "preprocessing" script does some initial outlier rejection.
// It also renames the ntuple(s).
// It is conventional for outfilename and pisafilename to have the same base.
// For record-keeping transparency, symlink the needed geometry file to
// run-prod-iter.par so the base matches the ROOT file.
// The alignment script VtxAlign.C requires a .par file with such a name.
// opt may contain:
//      "": vtx only. 
//      "cnt": vtxtrks & cnttrks. 
//      "fixed-bc" use beamcenter from rootfiles/bc-411768.root
void FilterData(const char *infilename = "rootfiles/anavtxcluster_411768-pro0-no-vtx2cnt.root",
                const char *outfilename = "rootfiles/411768-0-0.root",
                const char *configfilename = "production/config/config-zf-411768-0-0.txt",
                double vertexprobmin = 0.02,
                double vertexprobmax = 0.98,
                double maxdca = 0.5,
                double maxres_s = 0.1,
                double maxres_z = 0.1,
                int nhitsmin = 3,
                int nevents = 10000, // -1 = everything
                float frac4hit = -1, // -1 = no filter
                TString opt = "cnt, fixed-bc") // see above
{
  std::cout << "-- Opening " << infilename << " --" << std::endl;
  TFile *inFile = new TFile(infilename, "read");

  std::cout << "-- Reading seg_clusntuple --" << std::endl;
  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  TNtuple *svxcnt;
  assert(svxseg);

  std::cout << "-- Initializing SvxTGeo --" << std::endl;
  float bc[2] = {0};
  float e2w[3] = {0};
  float v2c[3] = {0};
  string geo;
  GetParamFromConfig(configfilename, bc, e2w, v2c, geo);
  SvxTGeo *tgeo = VTXModel(geo.c_str());

  TGraphErrors *gbc;
  if (opt.Contains("fixed-bc"))
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


  std::cout << "-- Creating output file " << outfilename << " --" << std::endl;
  TFile *outFile = new TFile(outfilename, "recreate");


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
    if (gbc)
    {
      Info("", "Using beamcenter in fit:");
      Info("", " E: (%.3f, %.3f)", gbc->GetX()[0], gbc->GetY()[0]);
      Info("", " W: (%.3f, %.3f)", gbc->GetX()[1], gbc->GetY()[1]);
    }
    FitTracks(vtxevents, gbc);

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
    if (gbc)
    {
      Info("", "Using beamcenter from %s", bcf->GetName());
      Info("", " E: (%.3f, %.3f)", gbc->GetX()[0], gbc->GetY()[0]);
      Info("", " W: (%.3f, %.3f)", gbc->GetX()[1], gbc->GetY()[1]);
    }
    FitTracks(vtxevents, gbc);

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
