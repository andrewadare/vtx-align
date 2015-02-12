#ifndef __ProdEventCutter_H__
#define __ProdEventCutter_H__


#include "SubsysReco.h"

/**
 * This module makes an event selection
 * Intended for use during prduction
 */
class ProdEventCutter: public SubsysReco
{

public:

    ProdEventCutter();
    ProdEventCutter(std::string filename);
    virtual ~ProdEventCutter() {}

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);
    void SetCentral(bool cent){m_central = cent;};
    void SetBBCqCut(float qcut){m_bbcq = qcut;};

protected:

    int m_pass;
    bool m_central;
    float m_bbcq;
    int m_runnumber;

};



#endif //__ProdEventCutter_H__
