///////////////////////////////////////////////////////////////////
//
// This Analysis Module is designed to provide information from
// SvxCentralTracks, SvxSegments reconstructed from compact format
//
///////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 10-16-2013
//
///////////////////////////////////////////////////////////////////
//
// Insiration from:
//  offline/analysis/svxanalysis/SvxAnalysis.cc
//  offline/analysis/svx_cent_ana/svxana.cc
//
///////////////////////////////////////////////////////////////////


#include <iostream>
#include <map>

#include "AnaAlignmentProd.h"

//phenix libraries
#include "getClass.h"

#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include "TriggerHelper.h"

#include <VtxOut.h>
#include <PHPoint.h>
#include "RunHeader.h"
#include "PreviousEvent.h"
#include "EventHeader.h"

#include "PHGlobal.h"
#include "PHCentralTrack.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"
#include "SvxClusterInfo.h"
#include "SvxSegmentList.h"
#include "SvxSegment.h"
#include "SvxClusterList.h"
#include "SvxCluster.h"

#include "ConversionVeto.h"

//root libraries
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TProfile.h"

static SvxCluster *searchSvxHit(SvxClusterList *d_svxhit, int clusterID);
static int findCloseSvxHits(SvxClusterList *, int, float, float, float[], float[], int[]);
static int layerNumber(float);
static const int MAXHIT = 50;

//===========================================================================
AnaAlignmentProd::AnaAlignmentProd() :
  d_phcnttrk(NULL),
  d_svxcnttrklist(NULL),
  d_vtxout(NULL),
  d_global(NULL),
  d_pevent(NULL),
  d_svxcluslist(NULL),
  d_svxseglist(NULL),
  m_fieldPolarity(-9999),
  m_runNumber(-9999),
  m_nEvent(0),
  m_nWrite(0),
  m_zerofield(false),
  m_convveto(NULL),
  OutFileName("AnaAlignmentProd_output.root"),
  OutputFile(NULL),
  m_print_pevent(true)
{
  ThisName = "AnaAlignmentProd";

  for (int i = 0; i < 4; i++)
    hcluster_phiz[i] = NULL;
  htrack_phized = NULL;
  hevent        = NULL;
  ntp_SVXCNT    = NULL;
  ntp_CNT       = NULL;
  ntp_SEG       = NULL;
  ntp_event     = NULL;
  //OutFileName = "AnaAlignmentProd_output.root";

  reset_variables();
}

//===========================================================================
AnaAlignmentProd::AnaAlignmentProd(std::string filename) :
  d_phcnttrk(NULL),
  d_svxcnttrklist(NULL),
  d_vtxout(NULL),
  d_global(NULL),
  d_pevent(NULL),
  d_svxcluslist(NULL),
  d_svxseglist(NULL),
  m_fieldPolarity(-9999),
  m_runNumber(-9999),
  m_nEvent(0),
  m_nWrite(0),
  m_zerofield(false),
  OutFileName(filename),
  OutputFile(NULL),
  m_print_pevent(true)
{
  ThisName = "AnaAlignmentProd";

  for (int i = 0; i < 4; i++)
    hcluster_phiz[i] = NULL;
  htrack_phized = NULL;
  hevent        = NULL;
  ntp_SVXCNT    = NULL;
  ntp_CNT       = NULL;
  ntp_SEG       = NULL;
  ntp_event     = NULL;

  reset_variables();
}


//===========================================================================
int AnaAlignmentProd::Init(PHCompositeNode *topNode)
{
  std::cout << "AnaAlignmentProd::Init()" << std::endl;

  OutputFile = new TFile(OutFileName.c_str(), "RECREATE");
  std::cout << "--> Opened file - " << OutFileName << std::endl;

  //-- Register Histograms --//

  for (int ilayer = 0; ilayer < 4; ilayer++)
  {
    hcluster_phiz[ilayer] = new TH2F(Form("hcluster_phized_B%i", ilayer),
                                     Form("B%i; z [cm]; phi", ilayer),
                                     400, -20, 20,
                                     800, -4, 4);

    hcluster_phiz_svxcnt[ilayer] =
      new TH2F(Form("hcluster_phized_svxcnt_B%i", ilayer),
               Form("B%i; z [cm]; phi", ilayer),
               400, -20, 20,
               800, -4, 4);
  }

  htrack_phized = new TH2F("htrack_phized",
                           "SvxCentralTracks; zed [cm]; phi",
                           500, -90, 90,
                           500, -1, 4);

  hevent = new TH1F("hevent", ";N_{events};", 1, 0.5, 1.5);
  //hm->registerHisto(hevent);

  for (int idc = 0; idc < 2; idc++)
  {
    pDC_alpha_phi[idc] = new TProfile(Form("pDC_alpha_phi_%i", idc), Form("DC Arm=%i;phi0;#alpha", idc), 500, -1, 4);
    hDC_alpha_phi[idc] = new TH2F(Form("hDC_alpha_phi_%i", idc), Form("DC Arm=%i;phi0;#alpha", idc), 500, -1, 4, 200, -0.02, 0.02);
  }

  //-- Create Ntuples --//


  //make the tree for SvxCentralTracks
  ntp_SVXCNT = new TTree("ntp_SVXCNT", "Tree containing SvxCentralTrack information");
  ntp_SVXCNT->Branch("run", &run, "run/I");
  ntp_SVXCNT->Branch("event", &event, "event/I");
  ntp_SVXCNT->Branch("eventseq", &eventseq, "eventseq/I");
  ntp_SVXCNT->Branch("svxindex", &svxindex, "svxindex/I");
  ntp_SVXCNT->Branch("dchindex", &dchindex, "dchindex/I");
  ntp_SVXCNT->Branch("dchquality", &dchquality, "dchquality/F");
  ntp_SVXCNT->Branch("mom", &mom, "mom/F");
  ntp_SVXCNT->Branch("pT", &pT, "pT/F");
  ntp_SVXCNT->Branch("E", &E, "E/F");
  ntp_SVXCNT->Branch("charge", &charge, "charge/I");
  ntp_SVXCNT->Branch("the0", &the0, "the0/F");
  ntp_SVXCNT->Branch("phi0", &phi0, "phi0/F");
  ntp_SVXCNT->Branch("dcarm", &dcarm, "dcarm/I");
  ntp_SVXCNT->Branch("emcsect", &emcsect, "emcsect/I");
  ntp_SVXCNT->Branch("sector", &sector, "sector/I");
  ntp_SVXCNT->Branch("centrality", &centrality, "centrality/F");
  ntp_SVXCNT->Branch("zed", &zed, "zed/F");
  ntp_SVXCNT->Branch("phi", &phi, "phi/F");
  ntp_SVXCNT->Branch("n0", &n0, "n0/I");
  ntp_SVXCNT->Branch("sn0", &sn0, "sn0/I");
  ntp_SVXCNT->Branch("chi2npe0", &chi2npe0, "chi2npe0/F");
  ntp_SVXCNT->Branch("schi2snpe0", &schi2snpe0, "schi2snpe0/F");
  ntp_SVXCNT->Branch("disp", &disp, "disp/F");
  ntp_SVXCNT->Branch("sdisp", &sdisp, "sdisp/F");
  ntp_SVXCNT->Branch("prob", &prob, "prob/F");
  ntp_SVXCNT->Branch("dep", &dep, "dep/F");
  ntp_SVXCNT->Branch("emcdphi", &emcdphi, "emcdphi/F");
  ntp_SVXCNT->Branch("emcdz", &emcdz, "emcdz/F");
  ntp_SVXCNT->Branch("emcsdphi", &emcsdphi, "emcsdphi/F");
  ntp_SVXCNT->Branch("emcsdz", &emcsdz, "emcsdz/F");
  ntp_SVXCNT->Branch("emcsdphi_e", &emcsdphi_e, "emcsdphi_e/F");
  ntp_SVXCNT->Branch("emcsdz_e", &emcsdz_e, "emcsdz_e/F");
  ntp_SVXCNT->Branch("chisq", &chisq, "chisq/F");
  ntp_SVXCNT->Branch("ndf", &ndf, "ndf/F");
  ntp_SVXCNT->Branch("quality", &quality, "quality/F");
  ntp_SVXCNT->Branch("score", &score, "score/F");
  ntp_SVXCNT->Branch("dca2D", &dca2D, "dca2D/F");
  ntp_SVXCNT->Branch("dcaz", &dcaz, "dcaz/F");
  ntp_SVXCNT->Branch("convVeto", &convVeto, "convVeto/O");
  ntp_SVXCNT->Branch("tickCut", &tickCut, "tickCut/O");
  ntp_SVXCNT->Branch("vtx", &vtx, "vtx[3]/F");
  ntp_SVXCNT->Branch("which_vtx", &which_vtx);
  ntp_SVXCNT->Branch("Nclus", &Nclus, "Nclus/I");
  ntp_SVXCNT->Branch("Nclus_layer", &Nclus_layer, "Nclus_layer[4]/I");
  ntp_SVXCNT->Branch("layer", &layer, "layer[Nclus]/I");
  ntp_SVXCNT->Branch("ladder", &ladder, "ladder[Nclus]/I");
  ntp_SVXCNT->Branch("res_s", &res_s, "res_s[Nclus]/F");
  ntp_SVXCNT->Branch("res_z", &res_z, "res_z[Nclus]/F");
  ntp_SVXCNT->Branch("clus_global_x", &clus_global_x, "clus_global_x[Nclus]/F");
  ntp_SVXCNT->Branch("clus_global_y", &clus_global_y, "clus_global_y[Nclus]/F");
  ntp_SVXCNT->Branch("clus_global_z", &clus_global_z, "clus_global_z[Nclus]/F");
  ntp_SVXCNT->Branch("hitPattern", &hitPattern, "hitPattern/B");
  ntp_SVXCNT->Branch("bbcz", &bbcz, "bbcz/F");
  ntp_SVXCNT->Branch("bbcq", &bbcq, "bbcq[2]/F");



  //make tree for PHCentralTracks
  ntp_CNT = new TTree("ntp_CNT", "Tree containing PHCentralTrack information");
  ntp_CNT->Branch("run", &run, "run/I");
  ntp_CNT->Branch("event", &event, "event/I");
  ntp_CNT->Branch("eventseq", &eventseq, "eventseq/I");
  ntp_CNT->Branch("centrality", &centrality, "centrality/F");
  ntp_CNT->Branch("vtx", &vtx, "vtx[3]/F");
  ntp_CNT->Branch("which_vtx", &which_vtx);
  ntp_CNT->Branch("bbcz", &bbcz, "bbcz/F");
  ntp_CNT->Branch("bbcq", &bbcq, "bbcq[2]/F");
  ntp_CNT->Branch("dchindex", &dchindex, "dchindex/I");
  ntp_CNT->Branch("dchquality", &dchquality, "dchquality/F");
  ntp_CNT->Branch("mom", &mom, "mom/F");
  ntp_CNT->Branch("pT", &pT, "pT/F");
  ntp_CNT->Branch("E", &E, "E/F");
  ntp_CNT->Branch("charge", &charge, "charge/I");
  ntp_CNT->Branch("the0", &the0, "the0/F");
  ntp_CNT->Branch("phi0", &phi0, "phi0/F");
  ntp_CNT->Branch("dcarm", &dcarm, "dcarm/I");
  ntp_CNT->Branch("emcsect", &emcsect, "emcsect/I");
  ntp_CNT->Branch("sector", &sector, "sector/I");
  ntp_CNT->Branch("zed", &zed, "zed/F");
  ntp_CNT->Branch("phi", &phi, "phi/F");
  ntp_CNT->Branch("n0", &n0, "n0/I");
  ntp_CNT->Branch("sn0", &sn0, "sn0/I");
  ntp_CNT->Branch("chi2npe0", &chi2npe0, "chi2npe0/F");
  ntp_CNT->Branch("schi2snpe0", &schi2snpe0, "schi2snpe0/F");
  ntp_CNT->Branch("disp", &disp, "disp/F");
  ntp_CNT->Branch("sdisp", &sdisp, "sdisp/F");
  ntp_CNT->Branch("prob", &prob, "prob/F");
  ntp_CNT->Branch("dep", &dep, "dep/F");
  ntp_CNT->Branch("emcdphi", &emcdphi, "emcdphi/F");
  ntp_CNT->Branch("emcdz", &emcdz, "emcdz/F");
  ntp_CNT->Branch("emcsdphi", &emcsdphi, "emcsdphi/F");
  ntp_CNT->Branch("emcsdz", &emcsdz, "emcsdz/F");
  ntp_CNT->Branch("emcsdphi_e", &emcsdphi_e, "emcsdphi_e/F");
  ntp_CNT->Branch("emcsdz_e", &emcsdz_e, "emcsdz_e/F");
  ntp_CNT->Branch("alpha", &alpha, "alpha/F");
  ntp_CNT->Branch("ppc1", &ppc1, "ppc1[3]/F");
  ntp_CNT->Branch("ppc2", &ppc2, "ppc2[3]/F");
  ntp_CNT->Branch("ppc3", &ppc3, "ppc3[3]/F");
  ntp_CNT->Branch("pc2dphi", &pc2dphi, "pc2dphi/F");
  ntp_CNT->Branch("pc2dz", &pc2dz, "pc2dz/F");
  ntp_CNT->Branch("pc3dphi", &pc3dphi, "pc3dphi/F");
  ntp_CNT->Branch("pc3dz", &pc3dz, "pc3dz/F");



  //make the tree for SvxSegments
  ntp_SEG = new TTree("ntp_SEG", "Tree containing SvxCentralTrack information");
  ntp_SEG->Branch("run", &run, "run/I");
  ntp_SEG->Branch("event", &event, "event/I");
  ntp_SEG->Branch("eventseq", &eventseq, "eventseq/I");
  ntp_SEG->Branch("svxindex", &svxindex, "svxindex/I");
  ntp_SEG->Branch("mom", &mom, "mom/F");
  ntp_SEG->Branch("pT", &pT, "pT/F");
  ntp_SEG->Branch("p", &p, "p[3]/F");
  ntp_SEG->Branch("charge", &charge, "charge/I");
  ntp_SEG->Branch("centrality", &centrality, "centrality/F");
  ntp_SEG->Branch("chisq", &chisq, "chisq/F");
  ntp_SEG->Branch("ndf", &ndf, "ndf/F");
  ntp_SEG->Branch("quality", &quality, "quality/F");
  ntp_SEG->Branch("score", &score, "score/F");
  ntp_SEG->Branch("dca2D", &dca2D, "dca2D/F");
  ntp_SEG->Branch("Nclus", &Nclus, "Nclus/I");
  ntp_SEG->Branch("Nclus_layer", &Nclus_layer, "Nclus_layer[4]/I");
  ntp_SEG->Branch("tickCut", &tickCut, "tickCut/O");
  ntp_SEG->Branch("vtx", &vtx, "vtx[3]/F");
  ntp_SEG->Branch("which_vtx", &which_vtx);
  ntp_SEG->Branch("bbcq", &bbcq, "bbcq[2]/F");
  ntp_SEG->Branch("bbcz", &bbcz, "bbcz/F");



  //make the tree for event information
  ntp_event = new TTree("ntp_event", "Tree containing event information");
  ntp_event->Branch("run", &run, "run/I");
  ntp_event->Branch("event", &event, "event/I");
  ntp_event->Branch("eventseq", &eventseq, "eventseq/I");
  ntp_event->Branch("vtx", &vtx, "vtx[3]/F");
  ntp_event->Branch("vtxE", &vtxE, "vtxE[3]/F");
  ntp_event->Branch("vtxW", &vtxW, "vtxW[3]/F");
  ntp_event->Branch("which_vtx", &which_vtx);
  ntp_event->Branch("centrality", &centrality, "centrality/F");
  ntp_event->Branch("bbcq", &bbcq, "bbcq[2]/F");
  ntp_event->Branch("bbcz", &bbcz, "bbcz/F");
  ntp_event->Branch("tickCut", &tickCut, "tickCut/O");
  ntp_event->Branch("nsvxcnt", &nsvxcnt, "nsvxcnt/I");
  ntp_event->Branch("nseg", &nseg, "nseg/I");
  ntp_event->Branch("ncnt", &ncnt, "ncnt/I");
  ntp_event->Branch("nclus", &nclus, "nclus/I");
  ntp_event->Branch("Nclus_layer", &Nclus_layer, "Nclus_layer[4]/I");


  //conversion veto function
  m_convveto = new ConversionVeto();


  //reset event counter
  m_runNumber = -1;
  m_nEvent    = 0;
  m_nWrite    = 0;

  return 0;
}

//===========================================================================
int AnaAlignmentProd::End(PHCompositeNode *topNode)
{

  std::cout << "AnaAlignmentProd::End()" << std::endl;
  std::cout << "--> Saved info for " << m_nWrite << " tracks." << std::endl;
  std::cout << "--> Writing results to " << OutFileName << std::endl;
  OutputFile->Write();

  OutputFile->Close();
  delete OutputFile;

  delete m_convveto;

  return 0;
}

//===========================================================================
int AnaAlignmentProd::InitRun(PHCompositeNode *topNode)
{

  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!runheader)
  {
    std::cout << "AnaAlignmentProd::InitRun() : No RunHeader, do nothing and return!" << std::endl;
    return 1;
  }


  //-- field polarity --//
  if (runheader->get_currentCentral() > 0)
  {
    m_fieldPolarity = 1.0;
  }
  else
  {
    m_fieldPolarity = -1.0;
  }
  std::cout << "AnaAlignmentProd::InitRun  fieldScale=" << m_fieldPolarity << std::endl;

  //-- run number --//
  m_runNumber = runheader->get_RunNumber();

  //reset event counter
  m_nEvent = 0;

  return 0;

}

//===========================================================================
int AnaAlignmentProd::process_event(PHCompositeNode *topNode)
{
  m_nEvent++;

  //-------------------EventHeader----------------------------------
  d_eventhead = NULL;
  d_eventhead = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!d_eventhead)
  {
      std::cout<< "ERROR!! EventHeader node not found." << std::endl;
      return -1;
  }
  //------------------------------------------------------------------

  //----------------------- VtxOut -----------------------------------//
  d_vtxout = NULL;
  d_vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if (!d_vtxout)
  {
    std::cout << "ERROR!! Can't find VtxOut!" << std::endl;
    return -1;
  }
  //------------------------------------------------------------------//

  //---------------------- PHGlobal ----------------------------------//
  d_global = NULL;
  d_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if (!d_global)
  {
    std::cout << "ERROR!! Can't find PHGlobal" << std::endl;
    return -1;
  }
  //------------------------------------------------------------------//


  //---------------------- SvxClusterList ----------------------------//
  d_svxcluslist = NULL;
  d_svxcluslist = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
  if (!d_svxcluslist)
  {
    std::cout << "ERROR!! Can't find SvxClusterList!" << std::endl;
    return -1;
  }
  //------------------------------------------------------------------//

  //---------------------- SvxSegmentList ----------------------------//
  d_svxseglist = NULL;
  d_svxseglist = findNode::getClass<SvxSegmentList>(topNode, "SvxSegmentList");
  if (!d_svxseglist)
  {
    std::cout << "ERROR!! Can't find SvxSegmentList!" << std::endl;
    return -1;
  }
  //------------------------------------------------------------------//

  //------------------ PreviousEvent ---------------------------------//
  d_pevent = NULL;
  d_pevent    = findNode::getClass<PreviousEvent>(topNode, "PreviousEvent");
  if (!d_pevent && m_print_pevent)
  {
    std::cout << "WARNING!! Can't find PreviousEvent! (suppressing further warnings)" << std::endl;
    m_print_pevent = false;
    //return -1;
  }
  //------------------------------------------------------------------//

  //---------------------- PHCentralTrack ----------------------------//
  d_phcnttrk = NULL;
  d_phcnttrk = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  // if (!d_phcnttrk)
  // {
  //   std::cout << "ERROR!! Can't find PHCentralTrack." << std::endl;
  //   return -1;
  // }
  //------------------------------------------------------------------//

  //---------------------- SvxCentralTrackList -----------------------//
  d_svxcnttrklist = NULL;
  d_svxcnttrklist = findNode::getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
  // if (!d_svxcnttrklist)
  // {
  //   std::cout << "ERROR!! Can't find SvxCentralTrackList." << std::endl;
  //   return -1;
  // }
  //------------------------------------------------------------------//


  //output some event info
  if (m_nEvent % 100000 == 0 )
  {
    std::cout << PHWHERE << " Event #" << m_nEvent << std::endl;
    //std::cout << PHWHERE << "        NSvxCentralTracks = " << d_svxcnttrklist->get_nCentralTracks() << std::endl;
  }

  m_EventSeqNumber = d_eventhead->get_EvtSequence();

  //---------------------------------------------//
  // CHECK TICK CUT
  // ABORT IF FALSE
  //---------------------------------------------//
  bool pass_tick = false;

  if (d_pevent)
  {
    int pticks[3] = {0};
    for ( int i = 0; i < 3; i++ )
      pticks[i] = d_pevent->get_clockticks(i);

    pass_tick =  !( ( 50 < pticks[0] && pticks[0] < 120) ||
                    (700 < pticks[1] && pticks[1] < 780) );
  }

  if (!pass_tick) return 0;


  //---------------------------------------------//
  // CHECK TRIGGERS
  // ABORT IF NOT FOUND
  //---------------------------------------------//
  TriggerHelper d_trghelp(topNode);

  //set up for Run 14 AuAu 200
  if (m_runNumber >= 405327 && m_runNumber <= 414988)
  {
    bool isnarrow = d_trghelp.didLevel1TriggerGetScaled("BBCLL1(>1 tubes) narrowvtx") ||
                    d_trghelp.didLevel1TriggerGetScaled("BBCLL1(>1 tubes) narrowvtx CopyA") ||
                    d_trghelp.didLevel1TriggerGetScaled("BBCLL1(>1 tubes) narrowvtx Copy B");

    if (!isnarrow) return 0;
  }



  //---------------------------------------------//
  // GET THE EVENT VERTEX AND CENTRALITY
  //---------------------------------------------//
  PHPoint vtxpos;
  vtxpos = d_vtxout->get_Vertex();

  //get the East vertex (if available, else fill with 0's)
  PHPoint vtxposE;
  vtxposE = d_vtxout->get_Vertex("SVX_PRECISEE");

  //get the West vertex (if available, else fill with 0's)
  PHPoint vtxposW;
  vtxposW = d_vtxout->get_Vertex("SVX_PRECISEW");



  double gcent = d_global->getCentrality();
  float bbc_qn = d_global->getBbcChargeN();
  float bbc_qs = d_global->getBbcChargeS();
  float bbc_z = (d_vtxout->get_Vertex("BBC")).getZ();

  //make sure a precise vertex was found
  std::string s_vtx = d_vtxout->which_Vtx();

  //make a BBC charge & z-vertex cut
  if (bbc_qn < 1 || bbc_qs < 1 || TMath::Abs(vtxpos.getZ()) > 10.)
    return 0;


  //---------------------------------------------//
  // FILL ALPHA V PHI TPROFILES
  //---------------------------------------------//
  if (d_phcnttrk)
  {
    for (unsigned int ipart = 0; ipart < d_phcnttrk->get_npart(); ipart++)
    {
      float dchquality = d_phcnttrk->get_quality(ipart);
      float zed        = d_phcnttrk->get_zed(ipart);
      float n0         = d_phcnttrk->get_n0(ipart);

      float alpha = d_phcnttrk->get_alpha(ipart);
      float phi0  = d_phcnttrk->get_phi0(ipart);
      int dcarm   = d_phcnttrk->get_dcarm(ipart);

      if (dcarm >= 0 && dcarm < 2 &&
          (dchquality == 31 || dchquality == 63) &&
          n0 < 1 &&
          fabs(zed) < 75)
      {
        pDC_alpha_phi[dcarm]->Fill(phi0, alpha);
        hDC_alpha_phi[dcarm]->Fill(phi0, alpha);
      }

    }
  }

  //---------------------------------------------//
  // CUT ON THE PRECISE VERTEX FOR REMAINING
  // INFORMATION
  //---------------------------------------------//
  /*
  if (s_vtx != "SVX_PRECISE")
      return 0;
  */

  //cound this event
  hevent->Fill(1);

  //---------------------------------------------//
  // FILL CLUSTER DISTRIBUTION HISTOGRAMS
  // AND COUNT THE NUMBER OF CLUSTERS IN EACH
  // LAYER
  //---------------------------------------------//
  int nclusB[4] = {0};
  for (int iclus = 0; iclus < d_svxcluslist->get_nClusters(); iclus++)
  {
    SvxCluster *clus = d_svxcluslist->get_Cluster(iclus);
    if (!clus) continue;

    int layer = clus->get_layer();
    float x = clus->get_xyz_global(0);
    float y = clus->get_xyz_global(1);
    float z = clus->get_xyz_global(2);
    float phi = TMath::ATan2(y, x);

    if (layer < 0 || layer > 3) continue;

    //fill histogram
    hcluster_phiz[layer]->Fill(z, phi);
    nclusB[layer]++;
  }


  //---------------------------------------------//
  // FILL EVENT TREE
  //---------------------------------------------//

  //fill the event ntuple
  reset_variables();
  run        = m_runNumber;
  event      = m_nEvent;
  eventseq   = m_EventSeqNumber;
  vtx[0]     = vtxpos.getX();
  vtx[1]     = vtxpos.getY();
  vtx[2]     = vtxpos.getZ();
  vtxE[0]    = vtxposE.getX();
  vtxE[1]    = vtxposE.getY();
  vtxE[2]    = vtxposE.getZ();
  vtxW[0]    = vtxposW.getX();
  vtxW[1]    = vtxposW.getY();
  vtxW[2]    = vtxposW.getZ();
  which_vtx  = s_vtx;
  centrality = gcent;
  //tick cut
  //false - passes the tick cut
  //true  - falis the tick cut
  tickCut = !pass_tick;
  bbcq[0] = bbc_qn;
  bbcq[1] = bbc_qs;
  bbcz    = bbc_z;
  if (d_svxcnttrklist)
    nsvxcnt = d_svxcnttrklist->get_nCentralTracks();
  else
    nsvxcnt = -1;
  nseg    = d_svxseglist->get_nSegments();
  if (d_phcnttrk)
    ncnt    = d_phcnttrk->get_npart();
  else
    ncnt = -1;
  nclus   = d_svxcluslist->get_nClusters();
  for (int i = 0; i < 4; i++)
    Nclus_layer[i] = nclusB[i];


  ntp_event->Fill();





  //---------------------------------------------//
  // Fill SVXCNT tree
  //---------------------------------------------//
  //std::cout << "--> Filling SvxCentralTrack details for " << svxcnttrklist->get_nCentralTracks() << " tracks." << std::endl;
  if (d_svxcnttrklist)
  {
    for (int itrk = 0; itrk < d_svxcnttrklist->get_nCentralTracks(); itrk++)
    {
      reset_variables();

      SvxCentralTrack *svxtrk = d_svxcnttrklist->getCentralTrack(itrk);
      if (!svxtrk) continue;

      //this index links the svx data to cnt traks
      unsigned int cntindex = svxtrk->getDchIndex();

      //double check that this central track exists
      if (cntindex >= d_phcnttrk->get_npart())
      {
        std:: cout << "WARNING!! Can not find cnt track for index=" << itrk << " (cntindex=" << cntindex << ")" << std::endl;
        continue;
      }

      //fill the ntuple
      run        = m_runNumber;
      event      = m_nEvent;
      eventseq   = m_EventSeqNumber;
      svxindex   = itrk;
      dchindex   = cntindex;
      dchquality = d_phcnttrk->get_quality(cntindex);
      mom        = d_phcnttrk->get_mom(cntindex);
      E          = d_phcnttrk->get_ecore(cntindex);

      float px                            = d_phcnttrk->get_mom(cntindex) * sin(d_phcnttrk->get_the0(cntindex)) * cos(d_phcnttrk->get_phi0(cntindex));
      float py                            = d_phcnttrk->get_mom(cntindex) * cos(d_phcnttrk->get_the0(cntindex)) * sin(d_phcnttrk->get_phi0(cntindex));
      pT                                  = sqrt(px * px + py * py);
      charge                              = d_phcnttrk->get_charge(cntindex);
      the0                                = d_phcnttrk->get_the0(cntindex);
      phi0                                = d_phcnttrk->get_phi0(cntindex);
      dcarm                               = d_phcnttrk->get_dcarm(cntindex);
      emcsect                             = d_phcnttrk->get_sect(cntindex);
      sector                              = dcarm * 4 + emcsect;
      centrality                          = gcent;
      zed                                 = d_phcnttrk->get_zed(cntindex);
      phi                                 = d_phcnttrk->get_phi(cntindex);
      n0                                  = d_phcnttrk->get_n0(cntindex);
      sn0                                 = d_phcnttrk->get_sn0(cntindex);
      if (d_phcnttrk->get_chi2(cntindex) >= 0 && d_phcnttrk->get_npe0(cntindex) >= 0)
        chi2npe0 = d_phcnttrk->get_chi2(cntindex) / d_phcnttrk->get_npe0(cntindex);
      if (d_phcnttrk->get_schi2(cntindex) >= 0 && d_phcnttrk->get_snpe0(cntindex) >= 0)
        schi2snpe0 = d_phcnttrk->get_schi2(cntindex) / d_phcnttrk->get_snpe0(cntindex);
      disp       = d_phcnttrk->get_disp(cntindex);
      sdisp      = d_phcnttrk->get_sdisp(cntindex);
      prob       = d_phcnttrk->get_prob(cntindex);
      dep        = d_phcnttrk->get_dep(cntindex);
      emcdphi    = d_phcnttrk->get_emcdphi(cntindex);
      emcdz      = d_phcnttrk->get_emcdz(cntindex);
      emcsdphi   = d_phcnttrk->get_emcsdphi(cntindex);
      emcsdz     = d_phcnttrk->get_emcsdz(cntindex);
      emcsdphi_e = d_phcnttrk->get_emcsdphi_e(cntindex);
      emcsdz_e   = d_phcnttrk->get_emcsdz_e(cntindex);
      chisq      = svxtrk->getChiSquare();
      ndf        = svxtrk->getNDF();
      quality    = svxtrk->getLinkQuality();
      score      = svxtrk->getLinkScore();
      dca2D      = svxtrk->getDCA2D();
      dcaz       = svxtrk->getDCAZ();
      vtx[0]     = vtxpos.getX();
      vtx[1]     = vtxpos.getY();
      vtx[2]     = vtxpos.getZ();
      which_vtx  = s_vtx;
      bbcq[0]    = bbc_qn;
      bbcq[1]    = bbc_qs;
      bbcz       = bbc_z;

      //conversion veto
      //false - not a conversion (passes pass_conversionVeto())
      //true  - conversion (fails pass_converstionVeto())
      //convVeto = !(pass_conversionVeto(svxtrk, charge));
      convVeto = !m_convveto->calculate(m_fieldPolarity,
                                        pT,
                                        charge,
                                        svxtrk,
                                        d_svxcluslist);

      //tick cut
      //false - passes the tick cut
      //true  - falis the tick cut
      tickCut = !pass_tick;


      //get cluster information
      Nclus = svxtrk->getNhits();
      Nclus_layer[0] = svxtrk->getNhits(0); //B0
      Nclus_layer[1] = svxtrk->getNhits(1); //B1

      //B2
      Nclus_layer[2] = 0;
      for (int i = 2; i < 5; i++)
        Nclus_layer[2] += svxtrk->getNhits(i);

      //B3
      Nclus_layer[3] = 0;
      for (int i = 5; i < 8; i++)
        Nclus_layer[3] += svxtrk->getNhits(i);

      hitPattern = svxtrk->getHitPattern();

      for (int ihit = 0; ihit < Nclus; ihit++)
      {
        SvxClusterInfo *svxclusinfo = svxtrk->getClusterInfo(ihit);

        if (!svxclusinfo)
        {
          std::cout << "ERROR: No ClusInfo Central Track" << std::endl;
          continue;
        }
        layer[ihit]  = svxclusinfo->getLayer();
        ladder[ihit] = svxclusinfo->getLadder();
        res_z[ihit]  = svxclusinfo->getdz();
        res_s[ihit]  = svxclusinfo->getdphi();
        clus_global_x[ihit] = svxclusinfo->getPosition(0);
        clus_global_y[ihit] = svxclusinfo->getPosition(1);
        clus_global_z[ihit] = svxclusinfo->getPosition(2);

        //fill cluster phi vs z histograms for SvxCentralTracks

        if ( ( (dchquality == 31 || dchquality == 63) &&
               chisq / ndf > 0 && chisq / ndf < 3 &&
               Nclus >= 3 &&
               Nclus_layer[0] > 0 && Nclus_layer[1] > 0 &&
               zed > -75 && zed < 75 &&
               (!m_zerofield && pT >= 1.0)
             ) )
        {
          float phi =
            TMath::ATan2(clus_global_y[ihit], clus_global_x[ihit]);
          hcluster_phiz_svxcnt[layer[ihit]]->Fill(clus_global_z[ihit], phi);

        }

      }


      //print_variables();
      //cuts regardless of field settings
      if ( ( (dchquality == 31 || dchquality == 63) &&
             chisq / ndf > 0 && chisq / ndf < 3 &&
             Nclus >= 3 &&
             Nclus_layer[0] > 0 && Nclus_layer[1] > 0 &&
             zed > -75 && zed < 75
           )
           //|| //swapped
           //( (dchquality == 31 || dchquality == 63) &&
           //  score > 60 &&
           //  chisq / ndf > 0 && chisq / ndf < 6 &&
           //  pT >= 1.0 &&
           //  dep>-2.5 &&
           //  sn0 >= 2 &&
           //  schi2snpe0 < 7 &&
           //  sdisp < 5)
         )
      {
        //if not zerofield, cut on pT
        if (!m_zerofield && pT >= 0.75)
        {
          ntp_SVXCNT->Fill();
          m_nWrite++;
          htrack_phized->Fill(zed, phi);
        }
        else if (m_zerofield)
        {
          ntp_SVXCNT->Fill();
          m_nWrite++;
          htrack_phized->Fill(zed, phi);
        }
      }
    }
  }

  //---------------------------------------------//
  // Fill PHCNT tree
  //---------------------------------------------//
  if (d_phcnttrk)
  {
    for (unsigned int itrk = 0; itrk < d_phcnttrk->get_npart(); itrk++)
    {
      reset_variables();

      //fill the ntuple
      run                             = m_runNumber;
      event                           = m_nEvent;
      eventseq                        = m_EventSeqNumber;
      vtx[0]                          = vtxpos.getX();
      vtx[1]                          = vtxpos.getY();
      vtx[2]                          = vtxpos.getZ();
      which_vtx                       = s_vtx;
      bbcq[0]                         = bbc_qn;
      bbcq[1]                         = bbc_qs;
      bbcz                            = bbc_z;
      centrality                      = gcent;
      dchindex                        = itrk;
      dchquality                      = d_phcnttrk->get_quality(itrk);
      mom                             = d_phcnttrk->get_mom(itrk);
      E                               = d_phcnttrk->get_ecore(itrk);
      float px                        = d_phcnttrk->get_mom(itrk) * sin(d_phcnttrk->get_the0(itrk)) * cos(d_phcnttrk->get_phi0(itrk));
      float py                        = d_phcnttrk->get_mom(itrk) * cos(d_phcnttrk->get_the0(itrk)) * sin(d_phcnttrk->get_phi0(itrk));
      pT                              = sqrt(px * px + py * py);
      charge                          = d_phcnttrk->get_charge(itrk);
      the0                            = d_phcnttrk->get_the0(itrk);
      phi0                            = d_phcnttrk->get_phi0(itrk);
      dcarm                           = d_phcnttrk->get_dcarm(itrk);
      emcsect                         = d_phcnttrk->get_sect(itrk);
      sector                          = dcarm * 4 + emcsect;
      zed                             = d_phcnttrk->get_zed(itrk);
      phi                             = d_phcnttrk->get_phi(itrk);
      n0                              = d_phcnttrk->get_n0(itrk);
      sn0                             = d_phcnttrk->get_sn0(itrk);
      if (d_phcnttrk->get_chi2(itrk) >= 0 && d_phcnttrk->get_npe0(itrk) >= 0)
        chi2npe0 = d_phcnttrk->get_chi2(itrk) / d_phcnttrk->get_npe0(itrk);
      if (d_phcnttrk->get_schi2(itrk) >= 0 && d_phcnttrk->get_snpe0(itrk) >= 0)
        schi2snpe0 = d_phcnttrk->get_schi2(itrk) / d_phcnttrk->get_snpe0(itrk);
      disp       = d_phcnttrk->get_disp(itrk);
      sdisp      = d_phcnttrk->get_sdisp(itrk);
      prob       = d_phcnttrk->get_prob(itrk);
      dep        = d_phcnttrk->get_dep(itrk);
      emcdphi    = d_phcnttrk->get_emcdphi(itrk);
      emcdz      = d_phcnttrk->get_emcdz(itrk);
      emcsdphi   = d_phcnttrk->get_emcsdphi(itrk);
      emcsdz     = d_phcnttrk->get_emcsdz(itrk);
      emcsdphi_e = d_phcnttrk->get_emcsdphi_e(itrk);
      emcsdz_e   = d_phcnttrk->get_emcsdz_e(itrk);
      alpha      = d_phcnttrk->get_alpha(itrk);
      ppc1[0]    = d_phcnttrk->get_ppc1x(itrk);
      ppc1[1]    = d_phcnttrk->get_ppc1y(itrk);
      ppc1[2]    = d_phcnttrk->get_ppc1z(itrk);
      ppc2[0]    = d_phcnttrk->get_ppc2x(itrk);
      ppc2[1]    = d_phcnttrk->get_ppc2y(itrk);
      ppc2[2]    = d_phcnttrk->get_ppc2z(itrk);
      ppc3[0]    = d_phcnttrk->get_ppc3x(itrk);
      ppc3[1]    = d_phcnttrk->get_ppc3y(itrk);
      ppc3[2]    = d_phcnttrk->get_ppc3z(itrk);
      pc2dphi    = d_phcnttrk->get_pc2dphi(itrk);
      pc2dz      = d_phcnttrk->get_pc2dz(itrk);
      pc3dphi    = d_phcnttrk->get_pc3dphi(itrk);
      pc3dz      = d_phcnttrk->get_pc3dz(itrk);

      //cuts regardless of field
      if ( (dchquality == 31 || dchquality == 63) &&
           zed > -75 && zed < 75 )
      {
        //cuts dependent on field
        if (!m_zerofield && pT >= 0.75)
        {
          ntp_CNT->Fill();
        }
        else if (m_zerofield)
        {
          ntp_CNT->Fill();
        }
      }
    }
  }

  //---------------------------------------------//
  // Fill SEG tree
  //---------------------------------------------//
  for (int itrk = 0; itrk < d_svxseglist->get_nSegments(); itrk++)
  {
    reset_variables();

    SvxSegment *svxtrk = d_svxseglist->get_segment(itrk);
    if (!svxtrk) continue;


    run        = m_runNumber;
    event      = m_nEvent;
    eventseq   = m_EventSeqNumber;
    svxindex   = itrk;
    centrality = gcent;
    vtx[0]     = vtxpos.getX();
    vtx[1]     = vtxpos.getY();
    vtx[2]     = vtxpos.getZ();
    which_vtx  = s_vtx;
    bbcq[0]    = bbc_qn;
    bbcq[1]    = bbc_qs;
    bbcz       = bbc_z;

    //tick cut
    //false - passes the tick cut
    //true  - falis the tick cut
    tickCut = !pass_tick;

    mom = svxtrk->getMomentum();

    float px       = svxtrk->get3Momentum(0);
    float py       = svxtrk->get3Momentum(1);
    pT             = sqrt(px * px + py * py);
    p[0]           = svxtrk->get3Momentum(0);
    p[1]           = svxtrk->get3Momentum(1);
    p[2]           = svxtrk->get3Momentum(2);
    charge         = svxtrk->IsPositive() ? 1 : -1;
    chisq          = svxtrk->getChiSq();
    ndf            = svxtrk->getNDF();
    quality        = svxtrk->getQuality();
    score          = svxtrk->getSegmentScore();
    dca2D          = svxtrk->getDCA2D();
    Nclus_layer[0] = svxtrk->getNhits(0);
    Nclus_layer[1] = svxtrk->getNhits(1);
    Nclus_layer[2] = svxtrk->getNhits(2);
    Nclus_layer[3] = svxtrk->getNhits(3);
    Nclus          = Nclus_layer[0] + Nclus_layer[1] + Nclus_layer[2] + Nclus_layer[3];

    if (!m_zerofield &&
        pT > 0.3 &&
        chisq / ndf > 0 && chisq / ndf < 6 &&
        Nclus >= 3 &&
        Nclus_layer[0] > 0 && Nclus_layer[1] > 0)
    {
      ntp_SEG->Fill();
    }
    else
      ntp_SEG->Fill();
  }



  return 0;
}

//===========================================================================
void AnaAlignmentProd::print_variables()
{
  std::cout << " run=" << run << std::endl;
  std::cout << " event=" << event << std::endl;
  std::cout << " svxindex=" << svxindex << std::endl;
  std::cout << " dchindex=" << dchindex << std::endl;
  std::cout << " dchquality=" << dchquality << std::endl;
  std::cout << " mom=" << mom << std::endl;
  std::cout << " pT=" << pT << std::endl;
  std::cout << " E=" << E << std::endl;
  std::cout << " charge=" << charge << std::endl;
  std::cout << " the0=" << the0 << std::endl;
  std::cout << " phi0=" << phi0 << std::endl;
  std::cout << " chisq=" << chisq << std::endl;
  std::cout << " ndf=" << ndf << std::endl;
  std::cout << " quality=" << quality << std::endl;
  std::cout << " score=" << score << std::endl;
  std::cout << " dca2D=" << dca2D << std::endl;
  std::cout << " dcaz=" << dcaz << std::endl;
  std::cout << " convVeto=" << convVeto << std::endl;
  std::cout << " tickCut=" << tickCut << std::endl;
  std::cout << " Nclus=" << Nclus << std::endl;
  for (int i = 0; i < 4; i++)
    std::cout << " Nclus_layer[" << i << "]=" << Nclus_layer[i] << std::endl;



}

//===========================================================================
void AnaAlignmentProd::reset_variables()
{
  run        = -9999;
  event      = -9999;
  eventseq   = -9999;
  svxindex   = -9999;
  dchindex   = -9999;
  dchquality = -9999;
  mom        = -9999;
  pT         = -9999;
  E          = -9999;
  charge     = -9999;
  the0       = -9999;
  phi0       = -9999;
  dcarm      = -9999;
  emcsect    = -9999;
  sector     = -9999;
  centrality = -9999.;
  zed        = -9999.;
  phi        = -9999.;
  n0         = -9999;
  sn0        = -9999;
  chi2npe0   = -9999.;
  schi2snpe0 = -9999.;
  disp       = -9999.;
  sdisp      = -9999.;
  prob       = -9999.;
  dep        = -9999.;
  emcdphi    = -9999.;
  emcdz      = -9999.;
  emcsdphi   = -9999.;
  emcsdz     = -9999.;
  emcsdphi_e = -9999.;
  emcsdz_e   = -9999.;
  chisq      = -9999;
  ndf        = -9999;
  quality    = -9999;
  score      = -9999;
  dca2D      = -9999;
  dcaz       = -9999;
  convVeto   = true;
  tickCut    = true;
  Nclus      = -9999;
  hitPattern = 0;
  bbcz       = -9999.;
  nsvxcnt    = -9999;
  nseg       = -9999;
  ncnt       = -9999;
  nclus      = -9999;
  alpha      = -9999.;
  pc2dphi    = -9999.;
  pc2dz      = -9999.;
  pc3dphi    = -9999.;
  pc3dz      = -9999.;

  for (int i = 0; i < 3; i++)
  {
    vtx[i]  = -9999.;
    vtxE[i] = -9999.;
    vtxW[i] = -9999.;
    p[i]    = -9999.;

    ppc1[i] = -9999.;
    ppc2[i] = -9999.;
    ppc3[i] = -9999.;
  }

  for (int i = 0; i < 4; i++)
  {
    Nclus_layer[i] = -9999;
  }

  which_vtx = "";

  for (int i = 0; i < 2; i++)
  {
    bbcq[i] = -9999.;
  }

  for (int i = 0; i < 7; i++)
  {
    layer[i]  = -9999;
    ladder[i] = -9999;
    res_s[i]  = -9999.;
    res_z[i]  = -9999.;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////
//
// ConversionVeto(): Taken, with some modification, from offline/AnalysisTrain/Run11VtxAna
//
//
// return: true if passes conversion veto cut, false otherwise
//
////////////////////////////////////////////////////////////////////////////////////////////
bool AnaAlignmentProd::pass_conversionVeto(SvxCentralTrack *svxcnt, float charge)
{

  // vars to save close hits
  int   vcid[MAXHIT];
  float vphi[MAXHIT], vz[MAXHIT];

  int B0p_near_hit = 0;
  int B0m_near_hit = 0;
  int B1p_near_hit = 0;
  int B1m_near_hit = 0;
  int B2_near_hit = 0;
  int B3p_near_hit = 0;
  int B3m_near_hit = 0;

  int nclust = svxcnt->getNhits();
  for (int ic = 0; ic < nclust; ic++)
  {
    SvxClusterInfo *info = svxcnt->getClusterInfo(ic);
    if (info == NULL) continue;

    int layer          = info->getLayer();
    int clusterID      = info->getClusterId();
    SvxCluster *svxhit = searchSvxHit(d_svxcluslist, clusterID);

    if (svxhit == NULL)
    {
      std::cout << "svxhit is not found" << std::endl;
    }
    else
    {
      float xhit = svxhit->get_xyz_global(0);
      float yhit = svxhit->get_xyz_global(1);
      float zhit = svxhit->get_xyz_global(2);

      //--float rhit = sqrt(xhit*xhit+yhit*yhit);
      float phihit  = atan2(yhit, xhit);
      if (phihit < - 1.5708) phihit += 6.283184;
      //--float phi  = atan2(yhit-y0,xhit-x0);
      //--if(phi < -1.5708) phi += 6.283184;

      //Look for other hits near this cluster
      int nfound = findCloseSvxHits(d_svxcluslist, layer, phihit, zhit, vphi, vz, vcid);

      for (int ifound = 0; ifound < nfound; ifound++)
      {
        //int   hitid = vcid[ifound];

        float dphi  = vphi[ifound] - phihit;
        float dz    = vz[ifound] - zhit;
        float cdphi = charge * dphi * (-m_fieldPolarity);

        //check there is a hit hear B0
        if (layer == 0 && fabs(dz) < 0.05 && 0.001 < cdphi && cdphi < 0.04)
          B0p_near_hit++;

        if (layer == 0 && fabs(dz) < 0.05 && -0.02 < cdphi && cdphi < -0.001)
          B0m_near_hit++;


        //check there is a hit hear B1
        if (layer == 1 && fabs(dz) < 0.05 && 0.001 < cdphi && cdphi < 0.06)
          B1p_near_hit++;

        if (layer == 1 && fabs(dz) < 0.05 && -0.02 < cdphi && cdphi < -0.001)
          B1m_near_hit++;


        //check there is a hit near B2
        if (layer == 2 && fabs(dz) < 0.1 &&
            ((0.002 < cdphi && cdphi < 0.08) || (-0.04 < cdphi && cdphi < -0.002)))
          B2_near_hit++;


        //check there is a hit hear B3
        if (layer == 3 && fabs(dz) < 0.1 && 0.001 < cdphi && cdphi < 0.08)
          B3p_near_hit++;

        if (layer == 3 && fabs(dz) < 0.1 && -0.02 < cdphi && cdphi < -0.001)
          B3m_near_hit++;


      }//for(ifound)
    }//if(svxhit)
  }//for(ic)


  return (B0m_near_hit == 0 && B0p_near_hit == 0 &&
          B1m_near_hit == 0 && B1p_near_hit == 0 &&
          B2_near_hit == 0  &&
          B3m_near_hit == 0 && B3p_near_hit == 0);

}

static SvxCluster *searchSvxHit(SvxClusterList *d_svxcluslist, int clusterID)
{
  if (d_svxcluslist == NULL)
  {
    std::cout << "searchSvxHit: d_svxcluslist is NULL!" << std::endl;
    return NULL;
  }
  else
  {
    int nhit = d_svxcluslist->get_nClusters();
    for (int i = 0; i < nhit; i++)
    {
      SvxCluster *svxhit = d_svxcluslist->get_Cluster(i);
      if (svxhit == NULL)
      {
        std::cout << "cluster NULL : " << i << std::endl;
        return NULL;
      }

      if (svxhit->get_hitID() == clusterID) return svxhit;
    }
    std::cout << "no cluster found in SvxHitMapEntry" << std::endl;
    return NULL;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////
//
// Taken from offline/AnalysisTrain/Run11VtxAna/Run11VtxAna.C, needed by
// pass_converstionVeto() above
//
////////////////////////////////////////////////////////////////////////////////////////////
static int findCloseSvxHits(SvxClusterList *d_svxcluslist,
                            int layer0, float phi0, float zhit0,
                            float vphi[], float vz[], int vclsid[])
{
  //const static float philimit[4]={0.08,0.08,0.08,0.08};
  const static float philimit[4] = {0.12, 0.12, 0.12, 0.12};
  if (d_svxcluslist == NULL) return 0;
  int nfound = 0;

  int nhit = d_svxcluslist->get_nClusters();
  for (int ihit = 0; ihit < nhit; ihit++)
  {
    SvxCluster *svxhit = d_svxcluslist->get_Cluster(ihit);
    int   cid  = svxhit->get_hitID();
    float xhit = svxhit->get_xyz_global(0);
    float yhit = svxhit->get_xyz_global(1);
    float zhit = svxhit->get_xyz_global(2);
    float rhit = sqrt(xhit * xhit + yhit * yhit);
    float phi  = atan2(yhit, xhit);
    if (phi < -1.5708) phi += 6.283184;
    int layer  = layerNumber(rhit);
    if (layer == layer0)
    {
      if (fabs(phi - phi0) < philimit[layer] && fabs(zhit - zhit0) < 0.4)
      {
        vphi[nfound]   = phi;
        vz[nfound]     = zhit;
        vclsid[nfound] = cid;
        nfound++;
      }
    }
    if (nfound >= MAXHIT) return nfound;
  }
  //  cout <<"nfound="<<nfound<<endl;
  return nfound;
}

////////////////////////////////////////////////////////////////////////////////////////////
//
// Taken from offline/AnalysisTrain/Run11VtxAna/Run11VtxAna.C, needed by
// pass_converstionVeto() above
//
////////////////////////////////////////////////////////////////////////////////////////////
static int layerNumber(float rhit)
{
  if (rhit < 3.0) return 0;
  else if (rhit < 6.0) return 1;
  else if (rhit < 14.0) return 2;
  else return 3;
}





