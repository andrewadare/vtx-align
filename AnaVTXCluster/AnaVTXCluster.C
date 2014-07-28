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

//root libraries
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNtuple.h"


float get_centrality(int bbcc);


//===========================================================================
AnaVTXCluster::AnaVTXCluster() :
    m_runnumber(0),
    nEvent(0)
{
    ThisName = "AnaVTXCluster";

    OutFileName = "AnaVTXCluster_output.root";
    std::cout<<"new build"<<std::endl;

}

//===========================================================================
AnaVTXCluster::AnaVTXCluster(std::string filename) :
    m_runnumber(0),
    nEvent(0)
{
    ThisName = "AnaVTXCluster";

    OutFileName = filename;


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

    standalone_tracks_per_event = new TH1F("standalone_tracks_per_event","standalone tracks per event",51,-0.5,50.5);
    cnt_tracks_per_event = new TH1F("cnt_tracks_per_event","cnt tracks per event",51,-0.5,50.5);
    //-- TTree containing track & associated cluster information for SvxCentralTracks --//
    /*svxcntclus = new TTree("svxcntclus", "Track and associated cluster information for SvxCentralTracks");
    svxcntclus->Branch("run", &run, "run/I");
    svxcntclus->Branch("event", &event, "event/I");
    svxcntclus->Branch("svxindex", &svxindex, "svxindex/I");
    svxcntclus->Branch("dchindex", &dchindex, "dchindex/I");
    svxcntclus->Branch("dchquality", &dchquality, "dchquality/I");
    svxcntclus->Branch("pT", &pT, "pT/F");
    svxcntclus->Branch("charge", &charge, "charge/F");
    svxcntclus->Branch("chisq", &chisq, "chisq/F");
    svxcntclus->Branch("ndf", &ndf, "ndf/F");
    svxcntclus->Branch("score", &score, "score/F");
    svxcntclus->Branch("vtx", &vtx, "vtx[3]/F");
    svxcntclus->Branch("Nclus", &Nclus, "Nclus/I");
    svxcntclus->Branch("hitPattern", &hitPattern, "hitPattern/B");
    svxcntclus->Branch("Nclus_layer", &Nclus_layer, "Nclus_layer[4]/I");
    svxcntclus->Branch("clus_layer", &clus_layer, "clus_layer[Nclus]/I");
    svxcntclus->Branch("clus_ladder", &clus_ladder, "clus_ladder[Nclus]/I");
    svxcntclus->Branch("clus_sensor", &clus_sensor, "clus_sensor[Nclus]/I");
    svxcntclus->Branch("clus_lx", &clus_lx, "clus_lx[Nclus]/F");
    svxcntclus->Branch("clus_ly", &clus_ly, "clus_ly[Nclus]/F");
    svxcntclus->Branch("clus_lz", &clus_lz, "clus_lz[Nclus]/F");
    svxcntclus->Branch("clus_gx", &clus_gx, "clus_gx[Nclus]/F");
    svxcntclus->Branch("clus_gy", &clus_gy, "clus_gy[Nclus]/F");
    svxcntclus->Branch("clus_gz", &clus_gz, "clus_gz[Nclus]/F");
    svxcntclus->Branch("clus_xsize",&clus_xsize,"clus_xsize[Nclus]/F");
    svxcntclus->Branch("clus_zsize",&clus_zsize,"clus_zsize[Nclus]/F");
    svxcntclus->Branch("centrality",&centrality,"centrality/F");

    cntntuple = new TTree("cnt_ntuple", "CNT Track information for alignment");
    cntntuple->Branch("gl_itrack",&gl_itrack,"gl_itrack/F");
    cntntuple->Branch("gl_phi0",&gl_phi0,"gl_phi0/F");
    cntntuple->Branch("gl_the0",&gl_the0,"gl_the0/F");
    cntntuple->Branch("gl_tr_mom",&gl_tr_mom,"gl_tr_mom/F");
    cntntuple->Branch("gl_pt",&gl_pt,"gl_pt/F");
    cntntuple->Branch("gl_charge",&gl_charge,"gl_charge/F");
    cntntuple->Branch("gl_emcdz",&gl_emcdz,"gl_emcdz/F");
    cntntuple->Branch("gl_emcdphi",&gl_emcdphi,"gl_emcdphi/F");
    cntntuple->Branch("gl_pc3dz",&gl_pc3dz,"gl_pc3dz/F");
    cntntuple->Branch("gl_pc3dphi",&gl_pc3dphi,"gl_pc3dphi/F");
    cntntuple->Branch("gl_pc2dz",&gl_pc2dz,"gl_pc2dz/F");
    cntntuple->Branch("gl_pc2dphi",&gl_pc2dphi,"gl_pc2dphi/F");
    cntntuple->Branch("gl_trk_quality",&gl_trk_quality,"gl_trk_quality/F");
    cntntuple->Branch("gl_n0",&gl_n0,"gl_n0/F");
    cntntuple->Branch("gl_ecore",&gl_ecore,"gl_ecore/F");
    cntntuple->Branch("gl_alpha",&gl_alpha,"gl_alpha/F");
    cntntuple->Branch("gl_zed",&gl_zed,"gl_zed/F");

    */
    //-- Track & associated cluster information for SvxSegment --//
    /*  svxsegclus = new TTree("svxsegclus", "Track & associated cluser information for SvxSegment");
    svxsegclus->Branch("run", &run, "run/I");
    svxsegclus->Branch("event", &event, "event/I");
    svxsegclus->Branch("svxindex", &svxindex, "svxindex/I");
    svxsegclus->Branch("pT", &pT, "pT/F");
    svxsegclus->Branch("phi",&phi,"phi/F");
    svxsegclus->Branch("theta",&theta,"theta/F");
    svxsegclus->Branch("charge", &charge, "charge/F");
    svxsegclus->Branch("chisq", &chisq, "chisq/F");
    svxsegclus->Branch("ndf", &ndf, "ndf/F");
    svxsegclus->Branch("score", &score, "score/F");
    svxsegclus->Branch("vtx", &vtx, "vtx[3]/F");
    svxsegclus->Branch("Nclus", &Nclus, "Nclus/I");
    svxsegclus->Branch("Nclus_layer", &Nclus_layer, "Nclus_layer[4]/I");
    svxsegclus->Branch("clus_layer", &clus_layer, "clus_layer[Nclus]/I");
    svxsegclus->Branch("clus_ladder", &clus_ladder, "clus_ladder[Nclus]/I");
    svxsegclus->Branch("clus_sensor", &clus_sensor, "clus_sensor[Nclus]/I");
    svxsegclus->Branch("clus_lx", &clus_lx, "clus_lx[Nclus]/F");
    svxsegclus->Branch("clus_ly", &clus_ly, "clus_ly[Nclus]/F");
    svxsegclus->Branch("clus_lz", &clus_lz, "clus_lz[Nclus]/F");
    svxsegclus->Branch("clus_gx", &clus_gx, "clus_gx[Nclus]/F");
    svxsegclus->Branch("clus_gy", &clus_gy, "clus_gy[Nclus]/F");
    svxsegclus->Branch("clus_gz", &clus_gz, "clus_gz[Nclus]/F");
    svxsegclus->Branch("clus_xsize",&clus_xsize,"clus_xsize[Nclus]/F");
    svxsegclus->Branch("clus_zsize",&clus_zsize,"clus_zsize[Nclus]/F");
    svxsegclus->Branch("centrality",&centrality,"centrality/F");
*/
      //        svxcntclus = new TNtuple("svxcntclus", "svxcntclus","event:cntindex:pT:Nhits:chisq:NDF:Nhits_B0:Nhits_B1:Nhits_B2:Nhits_B3:clusterID:layer:ladder:sensor:dphi:dz:lx:ly:lz");

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
    if (!phcnttrk)
    {
        std::cout << "ERROR!! Can't find PHCentralTrack." << std::endl;
        return -1;
    }

    //get the SvxCentralTrackList node
    SvxCentralTrackList *svxcnttrklist = findNode::getClass<SvxCentralTrackList>(topNode, "SvxCentralTrackList");
    if (!svxcnttrklist)
    {
        std::cout << "ERROR!! Can't find SvxCentralTrackList." << std::endl;
        return -1;
    }


    PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
    if(!global)
      {
	std::cout<<PHWHERE<<" ERROR::PHGlobal not found"<<std::endl;
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


    //output some event info


    centrality  = global->getCentrality();
    float bbc_qn      = global->getBbcChargeN();
    float bbc_qs      = global->getBbcChargeS();


    if(centrality==-9999)
      {
	centrality = get_centrality(bbc_qn+bbc_qs);

	//std::cout<<"error: no centrality; using bbc_q centrality "<<std::endl;
      }

    //    std::cout<<"centrality: "<<centrality<<std::endl;

    if(centrality < 50 )
      return -1;

    //---------------------------------------------//
    // GET THE EVENT VERTEX
    //---------------------------------------------//
    PHPoint vtxpos;
    vtxpos = vtxout->get_Vertex();
    float vtx_z = vtxpos.getZ();

    //std::cout<<"vtx_z: "<<vtx_z<<std::endl;

    if(vtx_z > 10 || vtx_z < -10)
      return -1;


    nEvent++;
    if(nEvent>500000) return -1;
    if (nEvent % 1000 == 0 )
        std::cout << "--> Event #" << nEvent << std::endl;
    /*
    //---------------------------------------------//
    // Fill the Cluster ntuple
    //---------------------------------------------//
    Float_t nclus[14];
    for (unsigned int iclus = 0; iclus < svxcluslist->get_nClusters(); iclus++)
    {
        SvxCluster *svxclus = svxcluslist->get_Cluster(iclus);

        if (svxclus)
        {
            nclus[0]  = m_runnumber;
            nclus[1]  = nEvent;
            nclus[2]  = iclus;//svxclus->getClusterId();
            nclus[3]  = svxclus->get_layer();
            nclus[4]  = svxclus->get_ladder();
            nclus[5]  = svxclus->get_sensor();
            nclus[6]  = svxclus->get_xyz_local(0);
            nclus[7]  = svxclus->get_xyz_local(1);
            nclus[8]  = svxclus->get_xyz_local(2);
            nclus[9]  = svxclus->get_xyz_global(0);
            nclus[10] = svxclus->get_xyz_global(1);
            nclus[11] = svxclus->get_xyz_global(2);
	    nclus[12] = svxclus->get_xz_size(0);
	    nclus[13] = svxclus->get_xz_size(1);
            clusntuple->Fill(nclus);
        }
    }*/
    /*
    unsigned int ncnttrack = phcnttrk->get_npart();

    for (unsigned int itrk = 0; itrk < ncnttrack; itrk++)
      {
	gl_itrack=itrk;
	//cout<<"gl_itrack="<<gl_itrack<<"   itrk="<<itrk<<endl;
	gl_tr_mom = phcnttrk->get_mom(itrk);

	gl_phi0 = phcnttrk->get_phi0(itrk);
	if(gl_phi0 > M_PI) gl_phi0 = gl_phi0 - 2.*M_PI;

	gl_the0 = phcnttrk->get_the0(itrk); //theta at the vertex

	gl_pt = fabs(gl_tr_mom)*sin(gl_the0);//pt ok

	gl_charge = phcnttrk->get_charge(itrk);//charge

	gl_emcdz = phcnttrk->get_emcdz(itrk);//emcz

	gl_emcdphi = phcnttrk->get_emcdphi(itrk);//emcdphi

	gl_pc3dz = phcnttrk->get_pc3dz(itrk);//pc3dz

	gl_pc3dphi = phcnttrk->get_pc3dphi(itrk);//pc3dphi

	gl_pc2dz = phcnttrk->get_pc2dz(itrk);//pc2 dz

	gl_pc2dphi = phcnttrk->get_pc2dphi(itrk);//pc2 dphi

	gl_zed= phcnttrk->get_zed(itrk);// Z coordinate at which the track crosses PC1

	gl_trk_quality = phcnttrk->get_quality(itrk);

	gl_n0 = phcnttrk->get_n0(itrk);

	gl_ecore=phcnttrk->get_ecore(itrk); //EMC "shower core" energy.

	gl_alpha = phcnttrk->get_alpha(itrk);

	cntntuple->Fill();

      }
    */
    //---------------------------------------------//
    // Fill Ntuple for clusters associated
    // with CNT tracks
    //---------------------------------------------//

    //    std::cout<<"n cnt tracks: "<<svxcnttrklist->get_nCentralTracks()<<std::endl;

    cnt_tracks_per_event->Fill(svxcnttrklist->get_nCentralTracks());
    for (unsigned int itrk = 0; itrk < svxcnttrklist->get_nCentralTracks(); itrk++)
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
	float mom   = phcnttrk->get_mom(cntindex);
	float the0  = phcnttrk->get_the0(cntindex);
	float phi0  = phcnttrk->get_phi0(cntindex);
	float px    = mom * sin(the0) * cos(phi0);
	float py    = mom * sin(the0) * sin(phi0);
	float pt    = sqrt(px * px + py * py);
	// int nhits   = svxtrk->getNhits();
	int score   = svxtrk->getLinkScore();
	int ndf     = svxtrk->getNDF();
	float chi2  = ndf < 1 ? 1e9 : svxtrk->getChiSquare() / ndf;
        int Nclus = svxtrk->getNhits();

	if ((dchqual == 31 || dchqual == 63) && /*pt > 0.5 && */score > 70 && chi2 < 4)
	  {
	    //	    std::cout<<"found a good central track"<<std::endl;
	    for (int ihit = 0; ihit < Nclus; ihit++)
	      {
		SvxClusterInfo *svxclusinfo = svxtrk->getClusterInfo(ihit);

		if(!svxclusinfo)
		  {
		    std::cout<<"ERROR: No ClusInfo Central Track"<<std::endl;
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
		if(!svxclus)
		  {
		    std::cout<<"ERROR: No Clus Central Track"<<std::endl;
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
		    nEvent
		  };
		//std::cout<<"filling cnt_clusntuple"<<std::endl;
		cnt_clusntuple->Fill(varr);
	      }
	  }
	cnt_trkid++;
      }

    standalone_tracks_per_event->Fill(svxseglist->get_nSegments());
    //std::cout<<"n segments: "<<svxseglist->get_nSegments()<<std::endl;
    for (unsigned int itrk = 0; itrk < svxseglist->get_nSegments(); itrk++)
      {
        SvxSegment *svxseg = svxseglist->get_segment(itrk);
        if (!svxseg)
	  {
            continue;
	  }

        run = m_runnumber;
        event = nEvent;
        svxindex = itrk;
        float px = svxseg->get3Momentum(0);
        float py = svxseg->get3Momentum(1);
        float pz = svxseg->get3Momentum(2);
        pT = sqrt(px * px + py * py);
        phi = TMath::ATan2(py,px);
        theta = TMath::ACos(pz/pT);

        charge = 1;
        if (!svxseg->IsPositive())
	  charge = -1;
        chisq = svxseg->getChiSq();
        ndf = svxseg->getNDF();
        score = svxseg->getSegmentScore();
        vtx[0] = vtxpos.getX();
        vtx[1] = vtxpos.getY();
	if(/*pT < 0.4 || */score < 70 )
	  continue;

	for(int ilayer = 0; ilayer < 4; ilayer++)
	  {
	    if(svxseg->getNhits(ilayer) == 1)
	      {
		SvxCluster *svxclus = svxcluslist->get_Cluster(svxseg->getClusterID(ilayer, 0));
		if(!svxclus)
		  {
		    std::cout<<"ERROR: No Clus Segment"<<std::endl;
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
		    nEvent
		  };
		seg_clusntuple->Fill(varr);
	      }

	  }
	seg_trkid++;
      }



    /*
    for (unsigned int itrk = 0; itrk < svxcnttrklist->get_nCentralTracks(); itrk++)
    {
        SvxCentralTrack *svxtrk = svxcnttrklist->getCentralTrack(itrk);

        //this index links the svx data to cnt traks
        int cntindex = svxtrk->getDchIndex();

        if (cntindex >= phcnttrk->get_npart())
        {
            std:: cout << "WARNING!! Can not find cnt track for index=" << itrk << " (cntindex=" << cntindex << ")" << std::endl;
            continue;
        }

        run = m_runnumber;
        event = nEvent;
        svxindex = itrk;
        dchindex = cntindex;


        float px = phcnttrk->get_mom(cntindex) * sin(phcnttrk->get_the0(cntindex)) * cos(phcnttrk->get_phi0(cntindex));
        float py = phcnttrk->get_mom(cntindex) * cos(phcnttrk->get_the0(cntindex)) * sin(phcnttrk->get_phi0(cntindex));
        pT = sqrt(px * px + py * py);
        charge = phcnttrk->get_charge(cntindex);

        chisq = svxtrk->getChiSquare();
        ndf = svxtrk->getNDF();
        score = svxtrk->getLinkScore();
        vtx[0] = vtxpos.getX();
        vtx[1] = vtxpos.getY();
        vtx[2] = vtxpos.getZ();


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

            if (svxclusinfo)
            {
                SvxCluster *svxclus = svxcluslist->get_Cluster(svxclusinfo->getClusterId());

                if (svxclus)
                {

                    //clus_layer[ihit] = svxclusinfo->getLayer();
                    //clus_ladder[ihit] = svxclusinfo->getLadder();
                    //clus_sensor[ihit] = svxclusinfo->getSensor();
                    //clus_gx[ihit] = svxclusinfo->getPosition(0);
                    //clus_gy[ihit] = svxclusinfo->getPosition(1);
                    //clus_gz[ihit] = svxclusinfo->getPosition(2);

                    clus_layer[ihit]  = svxclus->get_layer();
                    clus_ladder[ihit]  = svxclus->get_ladder();
                    clus_sensor[ihit]  = svxclus->get_sensor();
                    clus_lx[ihit]  = svxclus->get_xyz_local(0);
                    clus_ly[ihit]  = svxclus->get_xyz_local(1);
                    clus_lz[ihit]  = svxclus->get_xyz_local(2);
                    clus_gx[ihit]  = svxclus->get_xyz_global(0);
                    clus_gy[ihit]  = svxclus->get_xyz_global(1);
                    clus_gz[ihit]  = svxclus->get_xyz_global(2);
		    clus_xsize[ihit]   = svxclus->get_xz_size(0);
		    clus_zsize[ihit]   = svxclus->get_xz_size(1);

                }
            }
        }

        svxcntclus->Fill();


    }
    */
    //---------------------------------------------//
    // FILL INFO FOR CLUSTERS ASSOCIATED WITH
    // SVXSEGMENTS
    //---------------------------------------------//
    /*    for (unsigned int itrk = 0; itrk < svxseglist->get_nSegments(); itrk++)
    {
        SvxSegment *svxseg = svxseglist->get_segment(itrk);
        if (!svxseg)
        {
            continue;
        }

        run = m_runnumber;
        event = nEvent;
        svxindex = itrk;

        float px = svxseg->get3Momentum(0);
        float py = svxseg->get3Momentum(1);
	float pz = svxseg->get3Momentum(2);
        pT = sqrt(px * px + py * py);

	phi = TMath::ATan2(py,px);
	theta = TMath::ACos(pz/pT);

        charge = 1;
        if (!svxseg->IsPositive())
            charge = -1;

        chisq = svxseg->getChiSq();
        ndf = svxseg->getNDF();
        score = svxseg->getSegmentScore();
        vtx[0] = vtxpos.getX();
        vtx[1] = vtxpos.getY();
        vtx[2] = vtxpos.getZ();

        Nclus = 0;
        for (int ilayer = 0; ilayer < 4; ilayer++)
            Nclus += svxseg->getNhits(ilayer);

        Nclus_layer[0] = svxseg->getNhits(0);
        Nclus_layer[1] = svxseg->getNhits(1);
        Nclus_layer[2] = svxseg->getNhits(2);
        Nclus_layer[3] = svxseg->getNhits(3);



        for (int ilayer = 0; ilayer < 4; ilayer++)
        {
            for (int ihit = 0; ihit < 1; ihit++)
            {
	      if(Nclus_layer[ilayer]==1)
		{
                SvxCluster *svxclus = svxcluslist->get_Cluster(svxseg->getClusterID(ilayer, ihit));

                if (svxclus)
		  {

		      //clus_layer[ihit] = svxclusinfo->getLayer();
		      //clus_ladder[ihit] = svxclusinfo->getLadder();
		      //clus_sensor[ihit] = svxclusinfo->getSensor();
		      //clus_gx[ihit] = svxclusinfo->getPosition(0);
		     // clus_gy[ihit] = svxclusinfo->getPosition(1);
		    //  clus_gz[ihit] = svxclusinfo->getPosition(2);

		                      clus_layer[ilayer]  = svxclus->get_layer();
                    clus_ladder[ilayer]  = svxclus->get_ladder();
                    clus_sensor[ilayer]  = svxclus->get_sensor();
                    clus_lx[ilayer]  = svxclus->get_xyz_local(0);
                    clus_ly[ilayer]  = svxclus->get_xyz_local(1);
                    clus_lz[ilayer]  = svxclus->get_xyz_local(2);
                    clus_gx[ilayer]  = svxclus->get_xyz_global(0);
                    clus_gy[ilayer]  = svxclus->get_xyz_global(1);
                    clus_gz[ilayer]  = svxclus->get_xyz_global(2);
		    clus_xsize[ilayer] = svxclus->get_xz_size(0);
                    clus_zsize[ilayer] = svxclus->get_xz_size(1);
		    //cout<<"cluster: "<<ilayer<<" radius"<<sqrt(clus_gx

		  }
		}
	      else
		{
		  clus_layer[ilayer]  = -9999;
		  clus_ladder[ilayer]  = -9999;
		  clus_sensor[ilayer]  =-9999;
		  clus_lx[ilayer]  = -9999;
		  clus_ly[ilayer]  = -9999;
		  clus_lz[ilayer]  = -9999;
		  clus_gx[ilayer]  = -9999;
		  clus_gy[ilayer]  = -9999;
		  clus_gz[ilayer]  = -9999;
		  clus_xsize[ilayer] = -9999;
		  clus_zsize[ilayer] = -9999;
		}

            }

	}
	svxsegclus->Fill();
    }
*/


    return 0;
}


float get_centrality(int bbcc)
{
  float centrality;
  if      (bbcc <= 2 ) {centrality = 91.0; }
  else if (bbcc <= 5 ) {centrality = 87.5; }
  else if (bbcc <= 13 ) {centrality = 82.5; }
  else if (bbcc <= 24 ) {centrality = 77.5; }
  else if (bbcc <= 38 ) {centrality = 72.5; }
  else if (bbcc <= 59 ) {centrality = 67.5; }
  else if (bbcc <= 89 ) {centrality = 62.5; }
  else if (bbcc <= 130 ) {centrality = 57.5; }
  else if (bbcc <= 181 ) {centrality = 52.5; }
  else if (bbcc <= 245 ) {centrality = 47.5; }
  else if (bbcc <= 323 ) {centrality = 42.5; }
  else if (bbcc <= 416 ) {centrality = 37.5; }
  else if (bbcc <= 532 ) {centrality = 32.5; }
  else if (bbcc <= 660 ) {centrality = 27.5; }
  else if (bbcc <= 808 ) {centrality = 22.5; }
  else if (bbcc <= 987 ) {centrality = 17.5; }
  else if (bbcc <= 1199 ) {centrality = 12.5; }
  else if (bbcc <= 1448 ) {centrality =  7.5; }
  else {centrality =  2.5; }

  return centrality;

}



