#ifndef __AnaAlignmentProd_H__
#define __AnaAlignmentProd_H__


#include "SubsysReco.h"

class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class PHCentralTrack;
class SvxCentralTrackList;
class PHGlobal;
class VtxOut;
class PreviousEvent;
class SvxClusterList;
class PHSnglCentralTrack;
class SvxCentralTrack;
class TH1F;
class TH2F;
class SvxSegmentList;
class TProfile;
class ConversionVeto;
class EventHeader;

#include <string>
#include <map>

class AnaAlignmentProd: public SubsysReco
{

public:

    AnaAlignmentProd();
    AnaAlignmentProd(std::string filename);
    virtual ~AnaAlignmentProd() {}

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);
    void SetZeroField(bool zf){m_zerofield = zf;};


protected:

    //--nodes
    PHCentralTrack *d_phcnttrk;
    SvxCentralTrackList *d_svxcnttrklist;
    VtxOut *d_vtxout;
    PHGlobal *d_global;
    PreviousEvent *d_pevent;
    SvxClusterList *d_svxcluslist;
    SvxSegmentList *d_svxseglist;
    EventHeader *d_eventhead;


    //--functions
    void reset_variables();
    void print_variables();
    bool pass_conversionVeto(SvxCentralTrack *svxcnt, float charge);

    //--variables
    int m_fieldPolarity;
    int m_runNumber;
    int m_nEvent;
    int m_EventSeqNumber;
    int m_nWrite;
    bool m_zerofield;
    ConversionVeto *m_convveto;

    //--histograms
    TH2F *hcluster_phiz[4];
    TH2F *hcluster_phiz_svxcnt[4];
    TH2F *htrack_phized;
    TH1F *hevent;
    TProfile *pDC_alpha_phi[2];
    TH2F *hDC_alpha_phi[2];


    std::string OutFileName;

    TFile *OutputFile;

    //TTree containing information for each SvxCentralTrack
    TTree *ntp_SVXCNT;
    int run;
    int event;
    int eventseq;
    int svxindex;
    int dchindex;
    float dchquality;
    float mom;
    float pT;
    float E;
    int   charge;
    float the0;
    float phi0;
    int dcarm;
    int emcsect;
    int sector;
    float centrality;
    float zed;
    float phi;
    int n0;
    int sn0;
    float chi2npe0;
    float schi2snpe0;
    float disp;
    float sdisp;
    float prob;
    float dep;
    float emcdphi;
    float emcdz;
    float emcsdphi;
    float emcsdz;
    float emcsdphi_e;
    float emcsdz_e;
    float chisq;
    float ndf;
    float quality;
    float score;
    float dca2D;
    float dcaz;
    float vtx[3];
    bool convVeto;
    bool tickCut;
    int Nclus;
    char hitPattern;
    int Nclus_layer[4];
    int layer[7];
    int ladder[7];
    float res_s[7];
    float res_z[7];
    float clus_global_x[7];
    float clus_global_y[7];
    float clus_global_z[7];


    //TTree containing PHCentralTrack information
    TTree *ntp_CNT;
    float alpha;
    float ppc1[3]; //track projection at pc1 [x,y,z]
    float ppc2[3]; //track projection at pc2 [x,y,z]
    float ppc3[3]; //track projection at pc3 [x,y,z]
    float pc2dphi;
    float pc2dz;
    float pc3dphi;
    float pc3dz;

    //Tree containing SvxSegment information
    TTree *ntp_SEG;
    float p[3];

    //Tree containing event information
    TTree *ntp_event;
    //int run;
    //int event;
    //float vtx[3];
    std::string which_vtx;
    //float centrality;
    //bool tickCut;
    float vtxE[3];
    float vtxW[3];
    float bbcq[2];
    float bbcz;
    int nsvxcnt;
    int nseg;
    int ncnt;
    int nclus;


    bool m_print_pevent;

};



#endif //__AnaAlignmentProd_H__
