#include "VtxAlignBase.h"
#include "DataRejector.h"

// This "preprocessing" script does some initial outlier rejection.
// It also renames the ntuple(s).
// It is conventional for outfilename and pisafilename to have the same base.
// For record-keeping transparency, symlink the needed geometry file to
// run-prod-iter.par so the base matches the ROOT file.
// The alignment script VtxAlign.C requires a .par file with such a name.

void FilterData(const char *infilename = "rootfiles/zf-411768-0-0_3.root",
                const char *outfilename = "rootfiles/411768-0-0.root",
                const char *pisafilename = "geom/411768-0-0.par",
                double vertexprobmin = 0.02,
                double vertexprobmax = 0.98,
                double maxdca = 0.5,
                double maxres_s = 0.1,
                double maxres_z = 0.1,
                int nevents = -1, // -1 = everything
                TString opt = "cnt")
{
    std::cout << "-- Opening " << infilename << " --" << std::endl;
    TFile *inFile = new TFile(infilename, "read");

    std::cout << "-- Reading seg_clusntuple --" << std::endl;
    TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
    assert(svxseg);

    TNtuple *svxcnt;
    if (opt.Contains("cnt"))
    {
        std::cout << "-- Reading cnt_clusntuple --" << std::endl;
        svxcnt = (TNtuple *) inFile->Get("cnt_clusntuple");
        assert(svxcnt);
    }

    std::cout << "-- Creating output file " << outfilename << " --" << std::endl;
    TFile *outFile = new TFile(outfilename, "recreate");
    // TNtuple *vtxhits = new TNtuple("vtxhits", "VTX hit variables", HITVARS);
    // TNtuple *cnthits = new TNtuple("cnthits", "CNT hit variables", HITVARS);


    std::cout << "-- Creating vtxtrks --" << std::endl;
    TTree *vtxtrks = CreateTree("vtxtrks");

    TTree *cnttrks;
    if (opt.Contains("cnt"))
    {
        cnttrks = CreateTree("cnttrks");
    } 

    std::cout << "-- Initializing SvxTGeo --" << std::endl;
    SvxTGeo *tgeo = VTXModel(pisafilename);


    std::cout << "-- Reading events from ntuple --" << std::endl;
    geoEvents vtxevents;
    GetEventsFromClusTree(svxseg, tgeo, vtxevents, nevents);

    geoEvents cntevents;
    if (opt.Contains("cnt"))
    {
        GetEventsFromClusTree(svxcnt, tgeo, cntevents, nevents, "cnt");
        FillTree(cntevents, cnttrks);
    }

    std::cout << "-- Fitting tracks --" << std::endl;
    FitTracks(vtxevents);


    std::cout << "-- Filtering data --" << std::endl;
    FilterData(vtxevents,
               vertexprobmin,
               vertexprobmax,
               maxdca,
               maxres_s,
               maxres_z,
               vtxtrks);



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
