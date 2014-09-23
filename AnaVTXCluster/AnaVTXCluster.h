#ifndef __ANAVTXCLUSTER_H__
#define __ANAVTXCLUSTER_H__


#include "SubsysReco.h"

class PHCompositeNode;
class TFile;
class TH1F;
class TH2F;
class TNtuple;
class TTree;

#include <string>
#include <map>

class AnaVTXCluster: public SubsysReco
{

public:

    AnaVTXCluster();
    AnaVTXCluster(std::string filename);
    virtual ~AnaVTXCluster() {}

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void SetEventOffset(int off){event_offset = off;}

protected:

    void reset_variables();

    TH1F *standalone_tracks_per_event;
    TH1F *cnt_tracks_per_event;
    TNtuple *clusntuple; //ntuple containing cluster information from SvxCentralTrack

    TNtuple *seg_clusntuple;
    long int seg_trkid;
    TNtuple *cnt_clusntuple;
    long int cnt_trkid;

    //Tree containing event information
    TTree * ntp_event;
    int run;
    int event_offset;
    int event;
    float vtx[3];
    float vtxE[3];
    float vtxW[3];
    std::string which_vtx;
    float centrality;
    float bbcq[2];
    float bbcz;
    bool tickCut;


    bool m_print_pevent;
    int m_runnumber;
    int nEvent;

    std::string OutFileName;

    TFile *OutputFile;

};



#endif //__ANAVTXCLUSTER_H__
