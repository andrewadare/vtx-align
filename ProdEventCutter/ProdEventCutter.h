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

protected:

    int m_pass;

};



#endif //__ProdEventCutter_H__
