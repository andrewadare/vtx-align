#include "VtxAlignBase.h"
#include "VtxIO.h"
#include "DataRejector.h"

#include <utility>
#include <map>
#include <iomanip>

using namespace std;

void test_eventsync(const char *inFile = "rootfiles/anavtxcluster_0000411768-0000.root",
                    const char *geoFile = "geom/411768-0-0.par",
                    int nevents = 10000)
{

    //-- Get ntuples from file
    TFile *fin = TFile::Open(inFile);
    assert(fin);

    TNtuple *segclus = (TNtuple *) fin->Get("seg_clusntuple");
    TNtuple *cntclus = (TNtuple *) fin->Get("cnt_clusntuple");
    assert(segclus);
    assert(cntclus);

    SvxTGeo *tgeo = VTXModel(geoFile);


    //-- Create output file
    TFile *fout = new TFile("test_sync.root","RECREATE");

    //-- Read events from tree
    cout << endl;
    cout << "--> Getting unsynced events from ntuple" << endl;
    geoEvents segevt;
    map<int, int> segmap;
    GetEventsFromClusTree(segclus, tgeo, segevt, segmap, nevents);
    cout << "  n seg events = " << segevt.size() << endl;

    geoEvents cntevt;
    map<int, int> cntmap;
    GetEventsFromClusTree(cntclus, tgeo, cntevt, cntmap, nevents, "cnt");
    cout << "  n cnt events = " << cntevt.size() << endl;

    //-- Test original FilterData
    cout << endl;
    cout << "--> Filtering Unsynced seg data" << endl;
    //TTree *segunsynctree = CreateTree("segunsynctree");
    //FilterData(segevt, 0.01, 0.99, 0.2, 0.1, 0.1, segunsynctree);

    //-- Print a table of the number of tracks / event
    /*
    cout << endl;
    cout << "--> # of tracks per event" << endl;
    cout << setw(6) << "event"
         << setw(7) << "index"
         << setw(6) << "Nseg"
         << setw(6) << "Ncnt"
         << endl;
    map<int, int>::iterator cntpos;
    map<int, int>::iterator segpos;
    for (segpos = segmap.begin(); segpos != segmap.end(); ++segpos)
    {
        cout << setw(6) << segpos->first
             << setw(7) << segpos->second
             << setw(6) << segevt[segpos->second].size();

        cntpos = cntmap.find(segpos->first);
        if (cntpos != cntmap.end())
        {
            cout << setw(6) << cntevt[cntpos->second].size();
        }

        cout << endl;


    }
    */

    //-- Get Synced events from tree
    cout << "--> Getting synced events from ntuples" << endl;
    geoEvents segsyncevt;
    geoEvents cntsyncevt;
    GetSyncedEventsFromClusTree(segclus, cntclus, segsyncevt, cntsyncevt, tgeo, nevents);

    cout << "  n seg events = " << segsyncevt.size() << endl;
    cout << "  n cnt events = " << cntsyncevt.size() << endl;

    /*
    cout << endl;
    cout << "--> # of tracks per synced event" << endl;
    cout << setw(6) << "event"
         << setw(7) << "index"
         << setw(6) << "Nseg"
         << setw(6) << "Ncnt"
         << endl;
    for (int i = 0; i < segsyncevt.size(); i++)
    {
        cout << setw(6) << i
             << setw(7) << i
             << setw(6) << segsyncevt.at(i).size()
             << setw(6) << cntsyncevt.at(i).size()
             << endl;
    }
    */

    //-- Test data rejection
    cout << endl;
    cout << "--> Filterering Synced Data (" << segsyncevt.size() << " events)" << endl;

    TTree *segtree = CreateTree("segtree");
    cout << "  created tree " << segtree->GetName() << endl;

    TTree *cnttree = CreateTree("cnttree");
    cout << "  created tree " << cnttree->GetName() << endl;

    FitTracks(segsyncevt);
    FilterData(segsyncevt, cntsyncevt, segtree, cnttree);
    //FilterData(segsyncevt, 0.01, 0.99, 0.2, 0.1, 0.1, segtree);




}
