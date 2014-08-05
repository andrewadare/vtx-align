
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



//===========================================================================
ProdEventCutter::ProdEventCutter() :
    m_pass(0)
{
    ThisName = "ProdEventCutter";

}

//===========================================================================
ProdEventCutter::ProdEventCutter(std::string filename) :
    m_pass(0)
{
    ThisName = "ProdEventCutter";

}


//===========================================================================
int ProdEventCutter::Init(PHCompositeNode *topNode)
{
    std::cout << PHWHERE << std::endl;

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






    //---------------------------------------------//
    // MAKE EVENT CUT
    //---------------------------------------------//

    float bbc_qn      = d_global->getBbcChargeN();
    float bbc_qs      = d_global->getBbcChargeS();
    float bbcq = bbc_qn + bbc_qs;

    if (verbosity > 0)
        std::cout << PHWHERE << " - bbcq=" << bbcq << std::endl;

    if ( (bbc_qn + bbc_qs) >= 245)
        return ABORTEVENT;



    //---------------------------------------------//
    // EVENT IS OK
    //---------------------------------------------//
    if (verbosity > 0) std::cout << PHWHERE << " Passed Event selection" << std::endl;
    m_pass++;

    return EVENT_OK;
}

