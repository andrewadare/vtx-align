#include "VtxAlignBase.h"
#include "VtxIO.h"
#include <utility>
#include <map>
#include <iomanip>

using namespace std;

void test_eventsync(const char *inFile = "rootfiles/anavtxcluster_0000411768-0000.root",
                    const char *geoFile = "geom/411768-0-0.par",
                    int nevents = 8)
{

    //-- Get ntuples from file
    TFile *fin = TFile::Open(inFile);
    assert(fin);

    TNtuple *segclus = (TNtuple *) fin->Get("seg_clusntuple");
    TNtuple *cntclus = (TNtuple *) fin->Get("cnt_clusntuple");
    assert(segclus);
    assert(cntclus);

    SvxTGeo *tgeo = VTXModel(geoFile);


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

    //-- Print a table of the number of tracks / event
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

    //-- Get Synced events from tree
    cout << "--> Getting synced events from ntuples" << endl;
    geoEvents segsyncevt;
    geoEvents cntsyncevt;
    GetSyncedEventsFromClusTree(segclus, cntclus, segsyncevt, cntsyncevt, tgeo, nevents);

    cout << "  n seg events = " << segsyncevt.size() << endl;
    cout << "  n cnt events = " << cntsyncevt.size() << endl;

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




}
