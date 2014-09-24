///////////////////////////////////////////////////////////////////
//
// This Analysis Module is designed to provide information for
// looking at the efficiency of matching CNT tracks to
// VTX information using simulated DST's
//
///////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 3-1-2013
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
#include <stdio.h>
#include <stdlib.h>

#include "AnaVTXCluster.h"

//phenix libraries
#include "getClass.h"

#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"

#include "RunHeader.h"
#include "PHCentralTrack.h"
#include "SvxCentralTrackList.h"
#include "SvxCentralTrack.h"
#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxSegmentList.h"
#include "SvxSegment.h"
#include "VtxOut.h"
#include "PHPoint.h"
#include "PHGlobal.h"
#include "PreviousEvent.h"

//root libraries
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNtuple.h"


float get_centrality(int bbcc);


//===========================================================================
AnaVTXCluster::AnaVTXCluster() :
    m_print_pevent(true),
    m_runnumber(0),
    nEvent(0),
    OutFileName("AnaVTXCluster_output.root"),
    OutputFile(NULL)
{
    ThisName = "AnaVTXCluster";

    reset_variables();
    
    event_offset = 0;
}

//===========================================================================
AnaVTXCluster::AnaVTXCluster(std::string filename) :
    m_print_pevent(true),
    m_runnumber(0),
    nEvent(0),
    OutFileName(filename),
    OutputFile(NULL)
{
    ThisName = "AnaVTXCluster";

    reset_variables();
    
    event_offset = 0;

    //extract segNumber
    char*a = new char[filename.length()+1];
    strcpy(a,filename.c_str());

    char*b = strtok(a,"-");
    b = strtok(0,"-");
    b = strtok(0,"-");
    b = strtok(0,"-");
    b = strtok(0,"-");
    b = strtok(b,".");
    event_offset = atoi(b) + 0;
    std::cout<<"event offset is : "<<event_offset<<std::endl;
}


//===========================================================================
int AnaVTXCluster::Init(PHCompositeNode *topNode)
{
    std::cout << "AnaVTXCluster::Init()" << std::endl;

    OutputFile = new TFile(OutFileName.c_str(), "RECREATE");
    std::cout << "--> Opened file - " << OutFileName << std::endl;

    //Information for all clusters
    //clusntuple = new TNtuple("clusntuple", "All clusters", "run:event:clusterID:layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size");
    seg_clusntuple = new TNtuple("seg_clusntuple", "All clusters Associated with StandAloneTracks", "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid:event");
    cnt_clusntuple = new TNtuple("cnt_clusntuple", "All clusters Associated with SvxCentralTracks", "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid:trkphi:trktheta:event");

    standalone_tracks_per_event = new TH1F("standalone_tracks_per_event", "standalone tracks per event", 51, -0.5, 50.5);
    cnt_tracks_per_event = new TH1F("cnt_tracks_per_event", "cnt tracks per event", 51, -0.5, 50.5);


    //make the tree for event information
    ntp_event = new TTree("ntp_event", "Tree containing event information");
    ntp_event->Branch("run", &run, "run/I");
    ntp_event->Branch("event", &event, "event/I");
    ntp_event->Branch("vtx", &vtx, "vtx[3]/F");
    ntp_event->Branch("vtxE", &vtxE, "vtxE[3]/F");
    ntp_event->Branch("vtxW", &vtxW, "vtxW[3]/F");
    ntp_event->Branch("which_vtx", &which_vtx);
    ntp_event->Branch("centrality", &centrality, "centrality/F");
    ntp_event->Branch("bbcq", &bbcq, "bbcq[2]/F");
    ntp_event->Branch("bbcz", &bbcz, "bbcz/F");
    ntp_event->Branch("tickCut", &tickCut, "tickCut/O");



    //reset event counter
    nEvent = 0;

    return 0;
}

//===========================================================================
int AnaVTXCluster::End(PHCompositeNode *topNode)
{

    std::cout << "AnaVTXCluster::End()" << std::endl;

    std::cout << "--> Writing results to " << OutFileName << std::endl;
    OutputFile->Write();

    OutputFile->Close();
    delete OutputFile;



    return 0;
}

//===========================================================================
int AnaVTXCluster::InitRun(PHCompositeNode *topNode)
{

    RunHeader *runhdr = findNode::getClass<RunHeader>(topNode, "RunHeader");
    if (!runhdr)
    {
        std::cout << "AnaSvxCentralTracksTree::InitRun() : No RunHeader, do nothing and return!" << std::endl;
        return 1;
    }


    m_runnumber = runhdr->get_RunNumber();

    //reset event counter
    nEvent = 0;
    seg_trkid = 0;
    cnt_trkid = 0;
    return 0;

}

//===========================================================================
int AnaVTXCluster::process_event(PHCompositeNode *topNode)
{


    //get the PHCentralTrack node
    PHCentralTrack *phcnttrk = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");

    //get the SvxCentralTrackList node
    SvxCentralTrackList *svxcnttrklist = findNode::getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");


    PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
    if (!global)
    {
        std::cout << PHWHERE << " ERROR::PHGlobal not found" << std::endl;
        return -1;
    }

    //get the SvxSegmentList node
    SvxSegmentList *svxseglist = findNode::getClass<SvxSegmentList>(topNode, "SvxSegmentList");
    if (!svxseglist)
    {
        std::cout << "ERROR!! Can't find SvxSegmentList!" << std::endl;
        return -1;
    }

    //get the SvxClusterList node
    SvxClusterList *svxcluslist = findNode::getClass<SvxClusterList>(topNode, "SvxClusterList");
    if (!svxcluslist)
    {
        std::cout << "ERROR!! Can't find SvxClusterList." << std::endl;
        return -1;
    }

    //----------------------- VtxOut -----------------------------------//
    VtxOut *vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut");
    if (!vtxout)
    {
        std::cout << "ERROR!! Can't find VtxOut!" << std::endl;
        return -1;
    }
    //------------------------------------------------------------------//

    //------------------ PreviousEvent ---------------------------------//
    PreviousEvent *pevent    = findNode::getClass<PreviousEvent>(topNode, "PreviousEvent");
    if (!pevent && m_print_pevent)
    {
        std::cout << "WARNING!! Can't find PreviousEvent! (suppressing further warnings)" << std::endl;
        m_print_pevent = false;
        //return -1;
    }
    //------------------------------------------------------------------//


    //---------------------------------------------//
    // GET CENTRALITY
    //---------------------------------------------//
    centrality  = global->getCentrality();
    float bbc_qn      = global->getBbcChargeN();
    float bbc_qs      = global->getBbcChargeS();

    if (centrality == -9999)
    {
        centrality = get_centrality(bbc_qn + bbc_qs);

        //std::cout<<"error: no centrality; using bbc_q centrality "<<std::endl;
    }

    //---------------------------------------------//
    // CHECK TICK CUT
    //---------------------------------------------//
    bool pass_tick = false;

    if (pevent)
    {
        int pticks[3] = {0};
        for ( int i = 0; i < 3; i++ )
            pticks[i] = pevent->get_clockticks(i);

        pass_tick =  !( ( 50 < pticks[0] && pticks[0] < 120) ||
                        (700 < pticks[1] && pticks[1] < 780) );
    }

    //---------------------------------------------//
    // GET THE EVENT VERTEX
    //---------------------------------------------//
    PHPoint vtxpos;
    vtxpos = vtxout->get_Vertex();

    //get the East vertex (if available, else fill with 0's)
    PHPoint vtxposE;
    vtxposE = vtxout->get_Vertex("SVX_PRECISEE");

    //get the West vertex (if available, else fill with 0's)
    PHPoint vtxposW;
    vtxposW = vtxout->get_Vertex("SVX_PRECISEW");

    //make sure a precise vertex was found
    std::string s_vtx = vtxout->which_Vtx();


    //---------------------------------------------//
    // MAKE RECONSTRUCTED Z VERTEX CUT
    // NOTE: this isn't in ProdEventCutter
    // because the vertex is not yet reconstructed
    //---------------------------------------------//

    if (vtxpos.getZ() > 10 || vtxpos.getZ() < -10)
        return -1;

    //---------------------------------------------//
    // COUNT THE EVENT
    //---------------------------------------------//
    nEvent++;
    if (nEvent > 500000) return -1;
    if (nEvent % 1000 == 0 )
        std::cout << "--> Event #" << nEvent << std::endl;



    //---------------------------------------------//
    // FILL THE EVENT TREE
    //---------------------------------------------//
    reset_variables();
    run = m_runnumber;
    event = nEvent;
    vtx[0] = vtxpos.getX();
    vtx[1] = vtxpos.getY();
    vtx[2] = vtxpos.getZ();
    vtxE[0] = vtxposE.getX();
    vtxE[1] = vtxposE.getY();
    vtxE[2] = vtxposE.getZ();
    vtxW[0] = vtxposW.getX();
    vtxW[1] = vtxposW.getY();
    vtxW[2] = vtxposW.getZ();
    which_vtx = s_vtx;
    //centrality = gcent;
    bbcq[0] = bbc_qn;
    bbcq[1] = bbc_qs;
    bbcz = (vtxout->get_Vertex("BBC")).getZ();
    //tick cut
    //false - passes the tick cut
    //true  - fails the tick cut
    tickCut = !pass_tick;
    ntp_event->Fill();


    //---------------------------------------------//
    // Fill Ntuple for clusters associated
    // with SvxCentralTracks
    //---------------------------------------------//

    //    std::cout<<"n cnt tracks: "<<svxcnttrklist->get_nCentralTracks()<<std::endl;
    if (svxcnttrklist && phcnttrk)
    {
      int ngood_cent_tracks = 0;
      
      for (int itrk = 0; itrk < svxcnttrklist->get_nCentralTracks(); itrk++)
        {
	  SvxCentralTrack *svxtrk = svxcnttrklist->getCentralTrack(itrk);
	  
	  unsigned int cntindex = svxtrk->getDchIndex();
	  if (cntindex >= phcnttrk->get_npart())
	    {
	      std::cout << "WARNING No cnt track for index = " << itrk
			<< " (cntindex = " << cntindex << ")" << std::endl;
	      continue;
	    }
	  
	  int dchqual = phcnttrk->get_quality(cntindex);
	  float the0  = phcnttrk->get_the0(cntindex);
	  float phi0  = phcnttrk->get_phi0(cntindex);
	  // float mom   = phcnttrk->get_mom(cntindex);
	  // float px    = mom * sin(the0) * cos(phi0);
	  // float py    = mom * sin(the0) * sin(phi0);
	  // float pt    = sqrt(px * px + py * py);
	  // int nhits   = svxtrk->getNhits();
	  int score   = svxtrk->getLinkScore();
	  int ndf     = svxtrk->getNDF();
	  float chi2  = ndf < 1 ? 1e9 : svxtrk->getChiSquare() / ndf;
	  int Nclus = svxtrk->getNhits();
	  
	  if ((dchqual == 31 || dchqual == 63) && /*pt > 0.5 && */score > 70 && chi2 < 4)
            {
	      //      std::cout<<"found a good central track"<<std::endl;
	      bool clusters_bad = false;
	      for (int ihit = 0; ihit < Nclus; ihit++)
		{
		  SvxClusterInfo *svxclusinfo = svxtrk->getClusterInfo(ihit);
		  int clusxsize = svxclusinfo->getXZSize(0);
		  int cluszsize = svxclusinfo->getXZSize(1);
		  if(clusxsize > 5 || cluszsize > 5)
		    {
		      clusters_bad = true;
		      break;
		    }
		}
	      
	      if(clusters_bad) continue;
	      
	      ngood_cent_tracks++;
	      
	      for (int ihit = 0; ihit < Nclus; ihit++)
		{
		  SvxClusterInfo *svxclusinfo = svxtrk->getClusterInfo(ihit);
		  
		  if (!svxclusinfo)
                    {
		      std::cout << "ERROR: No ClusInfo Central Track" << std::endl;
		      continue;
                    }
		  float cx   = svxclusinfo->getPosition(0);
		  float cy   = svxclusinfo->getPosition(1);
		  float cz   = svxclusinfo->getPosition(2);
		  float rcyl = sqrt(cx * cx + cy * cy);
		  float zp   = svxclusinfo->getdphi();
		  float zs   = rcyl * zp;
		  float zz   = svxclusinfo->getdz();
		  int layer  = svxclusinfo->getLayer();
		  int ladder = svxclusinfo->getLadder();
		  int sensor = svxclusinfo->getSensor();
		  int clusxsize = svxclusinfo->getXZSize(0);
		  int cluszsize = svxclusinfo->getXZSize(1);
		  
		  SvxCluster *svxclus = svxcluslist->get_Cluster(svxclusinfo->getClusterId());
		  if (!svxclus)
                    {
		      std::cout << "ERROR: No Clus Central Track" << std::endl;
		      continue;
                    }
		  float clx  = svxclus->get_xyz_local(0);
		  float cly  = svxclus->get_xyz_local(1);
		  float clz  = svxclus->get_xyz_local(2);
		  
		  float varr[17] =
                    {
		      layer,
		      ladder,
		      sensor,
		      clx,
		      cly,
		      clz,
		      cx,
		      cy,
		      cz,
		      clusxsize,
		      cluszsize,
		      zz,
		      zs,
		      cnt_trkid,
		      phi0,
		      the0,
		      nEvent+event_offset*500000
                    };
		  //std::cout<<"filling cnt_clusntuple"<<std::endl;
		  cnt_clusntuple->Fill(varr);
                }
            }
            cnt_trkid++;
        }
        cnt_tracks_per_event->Fill(ngood_cent_tracks);
    }
    

    //---------------------------------------------//
    // Fill Ntuple for clusters associated
    // with SvxSegments
    //---------------------------------------------//

    int ngood_seg_tracks = 0;
    //std::cout<<"n segments: "<<svxseglist->get_nSegments()<<std::endl;
    for (int itrk = 0; itrk < svxseglist->get_nSegments(); itrk++)
    {
        SvxSegment *svxseg = svxseglist->get_segment(itrk);
        if (!svxseg)
        {
            continue;
        }

        float score = svxseg->getSegmentScore();
        if ( score < 70 )
            continue;

	bool clusters_bad = false;
	
        for (int ilayer = 0; ilayer < 4; ilayer++)
	  {
	    if (svxseg->getNhits(ilayer) == 1)
	      {

		SvxCluster *svxclus = svxcluslist->get_Cluster(svxseg->getClusterID(ilayer, 0));
		if (!svxclus)
		  {
		    std::cout << "ERROR: No Clus Segment" << std::endl;
		    continue;
		  }
		
		float clusxsize = svxclus->get_xz_size(0);
		float cluszsize = svxclus->get_xz_size(1);		
	
		if(clusxsize > 5 || cluszsize > 5)
		  {
		    clusters_bad = true;
		    break;
		  }
	      }
	  }
	
	if(clusters_bad) continue;

	ngood_seg_tracks++;

        for (int ilayer = 0; ilayer < 4; ilayer++)
        {
            if (svxseg->getNhits(ilayer) == 1)
            {
                    
	      SvxCluster *svxclus = svxcluslist->get_Cluster(svxseg->getClusterID(ilayer, 0));
	      if (!svxclus)
                {
		  std::cout << "ERROR: No Clus Segment" << std::endl;
		  continue;
                }
	      float ladder = svxclus->get_ladder();
	      float sensor = svxclus->get_sensor();
	      float clx = svxclus->get_xyz_local(0);
	      float cly = svxclus->get_xyz_local(1);
	      float clz = svxclus->get_xyz_local(2);
	      float cx = svxclus->get_xyz_global(0);
	      float cy = svxclus->get_xyz_global(1);
	      float cz = svxclus->get_xyz_global(2);
	      float clusxsize = svxclus->get_xz_size(0);
	      float cluszsize = svxclus->get_xz_size(1);
	      
	      float varr[15] =
                {
		  ilayer,
		  ladder,
		  sensor,
		  clx,
		  cly,
		  clz,
		  cx,
		  cy,
		  cz,
		  clusxsize,
		  cluszsize,
		  0,
		  0,
		  seg_trkid,
		  nEvent+event_offset*500000
                };
	      seg_clusntuple->Fill(varr);
            }
	    
        }
        seg_trkid++;
    }
    standalone_tracks_per_event->Fill(ngood_seg_tracks);
    
    
    return 0;
}


float get_centrality(int bbcc)
{
    float centrality;
    if      (bbcc <= 2 )
    {
        centrality = 91.0;
    }
    else if (bbcc <= 5 )
    {
        centrality = 87.5;
    }
    else if (bbcc <= 13 )
    {
        centrality = 82.5;
    }
    else if (bbcc <= 24 )
    {
        centrality = 77.5;
    }
    else if (bbcc <= 38 )
    {
        centrality = 72.5;
    }
    else if (bbcc <= 59 )
    {
        centrality = 67.5;
    }
    else if (bbcc <= 89 )
    {
        centrality = 62.5;
    }
    else if (bbcc <= 130 )
    {
        centrality = 57.5;
    }
    else if (bbcc <= 181 )
    {
        centrality = 52.5;
    }
    else if (bbcc <= 245 )
    {
        centrality = 47.5;
    }
    else if (bbcc <= 323 )
    {
        centrality = 42.5;
    }
    else if (bbcc <= 416 )
    {
        centrality = 37.5;
    }
    else if (bbcc <= 532 )
    {
        centrality = 32.5;
    }
    else if (bbcc <= 660 )
    {
        centrality = 27.5;
    }
    else if (bbcc <= 808 )
    {
        centrality = 22.5;
    }
    else if (bbcc <= 987 )
    {
        centrality = 17.5;
    }
    else if (bbcc <= 1199 )
    {
        centrality = 12.5;
    }
    else if (bbcc <= 1448 )
    {
        centrality =  7.5;
    }
    else
    {
        centrality =  2.5;
    }

    return centrality;

}


//===========================================================================
void AnaVTXCluster::reset_variables()
{
    run = -9999;
    event = -9999;
    centrality = -9999.;
    tickCut = true;
    bbcz = -9999.;

    for (int i = 0; i < 3; i++)
    {
        vtx[i] = -9999.;
        vtxE[i] = -9999.;
        vtxW[i] = -9999.;
    }

    which_vtx = "";

    for (int i = 0; i < 2; i++)
    {
        bbcq[i] = -9999.;
    }
}


