#include <stdio.h>
#include <time.h>
#include <string>
#include <iostream>

char output[200];
vector<string> outfiles;

void MakeOutput(int runnumber, int segment, const char *file = "", const char *dir = ".")
{
    char mkdircmd[100];
    sprintf(mkdircmd, "mkdir -p %s", dir);
    sprintf(output, "%s/%s_%s-%010d-%04d.root", dir, file, gSystem->Getenv("PRODTAG"), runnumber, segment);
    gSystem->Exec(mkdircmd);
    outfiles.push_back(output);
    return output;
    //  sprintf(output,"%s_%s-%010d-%04d.root",file,gSystem->Getenv("PRODTAG"),runnumber,segment);
}

void MakePRDFOutput(int runnumber, int segment, const char *file = "")
{
    sprintf(output, "%s_%s-%010d-%04d.PRDFF", file, gSystem->Getenv("PRODTAG"), runnumber, segment);
    char mkdircmd[100];
    sprintf(mkdircmd, "mkdir -p %s", file);
    gSystem->Exec(mkdircmd);
    outfiles.push_back(output);
    //  sprintf(output,"%s_%s-%010d-%04d.PRDF",file,gSystem->Getenv("PRODTAG"),runnumber,segment);
}

void DST_SVX(const int runnumber, const int segment, const char *dir, const char *file, const char *trgsel = 0)
{

    MakeOutput(runnumber, segment, file, dir);
    Fun4AllServer *se = Fun4AllServer::instance();
    Fun4AllDstOutputManager *OutManager  = new Fun4AllDstOutputManager(file, output);

    if (trgsel)
    {
        OutManager->AddEventSelector(trgsel);
    }
    OutManager->AddNode("EventHeader");
    OutManager->AddNode("Sync");
    OutManager->AddNode("TrigLvl1");
    OutManager->AddNode("PreviousEvent");
    OutManager->AddNode("ErtOut");
    OutManager->AddNode("VtxOut");
    OutManager->AddNode("BbcOut");
    OutManager->AddNode("BbcRaw");
    OutManager->AddNode("PHGlobal");
    OutManager->AddNode("PHGlobal_CENTRAL");
    OutManager->AddNode("DchHit_VarArray");
    OutManager->AddNode("EmcHit_VarArray");
    OutManager->AddNode("Pc1Hit_VarArray");
    OutManager->AddNode("Pc2Hit_VarArray");
    OutManager->AddNode("Pc3Hit_VarArray");
    OutManager->AddNode("TofeHit_VarArray");
    OutManager->AddNode("TofwHit_VarArray");
    OutManager->AddNode("CrkHit_VarArray");
    OutManager->AddNode("AccHit_VarArray");
    OutManager->AddNode("CglTrackHits_VarArray");
    OutManager->AddNode("CglTrackBackHits_VarArray");
    OutManager->AddNode("TrackProjection_VarArray");
    OutManager->AddNode("TrackLineProjection_VarArray");
    OutManager->AddNode("TrackPathLength_VarArray");
    OutManager->AddNode("emcHitContainer");
    OutManager->AddNode("SvxEventInfo");
    OutManager->AddNode("SvxClusterList");
    OutManager->AddNode("SvxSegmentList");
    OutManager->AddNode("SvxCentralTrackList");

    //OutManager->AddNode("SvxHit_VarArray");
    //OutManager->AddNode("SvxTrack_VarArray");
    //OutManager->AddNode("SvxCentralTrack_VarArray");
    //OutManager->AddNode("SvxCentralTrackBG_VarArray");

    //  OutManager->AddNode("SvxCentralTrackBackList");


    se->registerOutputManager(OutManager);
}

void SVXQA_IOManager(const int runnumber, const int segment, const char *dir)
{
    MakeOutput(runnumber, segment, "SVXQA", dir);
    Fun4AllServer *se = Fun4AllServer::instance();
    Fun4AllHistoManager *hm = se->getHistoManager("SVXQA");
    if (hm) hm->dumpHistos(output);
}
