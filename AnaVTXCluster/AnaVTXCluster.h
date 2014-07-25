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


protected:


    TH1F *standalone_tracks_per_event;
    TH1F *cnt_tracks_per_event;
    TNtuple *clusntuple; //ntuple containing cluster information from SvxCentralTrack

    TNtuple *seg_clusntuple;
    long int seg_trkid;
    TNtuple *cnt_clusntuple;
    long int cnt_trkid;
    //-- Track & cluster information for SvxCentralTracks --//
    TTree *svxcntclus;
    int run;
    float res_z;
    float res_s;
    int event;
    int svxindex;
    int dchindex;
    int dchquality;
    float pT;
    float phi;
    float theta;
    float charge;
    float chisq;
    float ndf;
    float score;
    float vtx[3];
    int Nclus;
    char hitPattern;
    int Nclus_layer[4];
    int clus_layer[7];
    int clus_ladder[7];
    int clus_sensor[7];
    float clus_lx[7];
    float clus_ly[7];
    float clus_lz[7];
    float clus_gx[7];
    float clus_gy[7];
    float clus_gz[7];
    float clus_xsize[7];
    float clus_zsize[7];
    //-- Track & cluster information for SvxSegment --//
    // shares many of the above variables
    TTree* svxsegclus;

    TTree* cntntuple;
    float gl_itrack;
    float gl_phi0;
    float gl_tr_mom;
    float gl_the0;
    float gl_pt;
    float gl_charge;
    float gl_emcdz;
    float gl_emcdphi;
    float gl_pc3dz;
    float gl_pc3dphi;
    float gl_pc2dz;
    float gl_pc2dphi;
    float gl_trk_quality;
    float gl_n0;
    float gl_ecore;
    float gl_alpha;
    float gl_zed;
    float centrality;
    bool found_phcnttrk;

    int m_runnumber;
    int nEvent;

    std::string OutFileName;

    TFile *OutputFile;

};



#endif //__ANAVTXCLUSTER_H__
