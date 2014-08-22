
#include <iostream>

#include "ProdEventCutter.h"

//phenix libraries
#include "getClass.h"
#include "PHTypedNodeIterator.h"
#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <TrigLvl1.h>
#include <TriggerHelper.h>
#include <VtxOut.h>
#include <PHPoint.h>
#include "RunHeader.h"
#include "PHGlobal.h"
#include "PreviousEvent.h"
#include "TriggerHelper.h"


//===========================================================================
ProdEventCutter::ProdEventCutter() :
    m_pass(0),
    m_central(false),
    m_bbcq(200)
{
    ThisName = "ProdEventCutter";

}

//===========================================================================
ProdEventCutter::ProdEventCutter(std::string filename) :
    m_pass(0),
    m_central(false),
    m_bbcq(200)
{
    ThisName = "ProdEventCutter";

}


//===========================================================================
int ProdEventCutter::Init(PHCompositeNode *topNode)
{
    std::cout << PHWHERE << std::endl;
    if (!m_central)
        std::cout << PHWHERE << " - cutting on events bbcq < " << m_bbcq << std::endl;
    else if (m_central)
        std::cout << PHWHERE << " - cutting on events bbcq > " << m_bbcq << std::endl;

    return 0;
}

//===========================================================================
int ProdEventCutter::End(PHCompositeNode *topNode)
{

    std::cout << PHWHERE << " - Passed " << m_pass << " events." << std::endl;

    return 0;
}

//===========================================================================
int ProdEventCutter::InitRun(PHCompositeNode *topNode)
{

    return 0;

}

//===========================================================================
int ProdEventCutter::process_event(PHCompositeNode *topNode)
{

    //---------------------------------------------//
    // GET REQUIRED NODES
    //---------------------------------------------//


    //---------------------- PHGlobal ----------------------------------//
    PHGlobal *d_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
    if (!d_global)
    {
        std::cout << "ERROR!! Can't find PHGlobal" << std::endl;
        return -1;
    }
    //------------------------------------------------------------------//

    //------------------ PreviousEvent ---------------------------------//
    PreviousEvent *pevent    = findNode::getClass<PreviousEvent>(topNode, "PreviousEvent");
    if (!pevent)
    {
        std::cout << "ERROR!! Can't find PreviousEvent! (suppressing further warnings)" << std::endl;
        return -1;
    }
    //------------------------------------------------------------------//





    //---------------------------------------------//
    // MAKE CENTRALITY CUT
    //---------------------------------------------//

    float bbc_qn      = d_global->getBbcChargeN();
    float bbc_qs      = d_global->getBbcChargeS();
    float bbcq = bbc_qn + bbc_qs;

    if (verbosity > 0)
        std::cout << PHWHERE << " - bbcq=" << bbcq << std::endl;

    if (!m_central && (bbc_qn + bbc_qs) > m_bbcq)
        return ABORTEVENT;
    else if (m_central && (bbc_qn + bbc_qs) < m_bbcq)
        return ABORTEVENT;


    //---------------------------------------------//
    // MAKE TICK CUT
    //---------------------------------------------//
    int pticks[3] = {0};
    for ( int i = 0; i < 3; i++ )
        pticks[i] = pevent->get_clockticks(i);

    bool pass_tick =  !( ( 50 < pticks[0] && pticks[0] < 120) ||
                         (700 < pticks[1] && pticks[1] < 780) );

    if (!pass_tick)
        return ABORTEVENT;


    //---------------------------------------------//
    // CHECK TRIGGERS
    // ABORT IF NOT FOUND
    //---------------------------------------------//
    TriggerHelper d_trghelp(topNode);

    //set up for Run 14 AuAu 200
    bool isnarrow = d_trghelp.didLevel1TriggerGetScaled("BBCLL1(>1 tubes) narrowvtx") ||
                    d_trghelp.didLevel1TriggerGetScaled("BBCLL1(>1 tubes) narrowvtx CopyA") ||
                    d_trghelp.didLevel1TriggerGetScaled("BBCLL1(>1 tubes) narrowvtx Copy B");

    if (!isnarrow) return 0;


    //---------------------------------------------//
    // EVENT IS OK
    //---------------------------------------------//
    if (verbosity > 0) std::cout << PHWHERE << " Passed Event selection" << std::endl;
    m_pass++;

    return EVENT_OK;
}

