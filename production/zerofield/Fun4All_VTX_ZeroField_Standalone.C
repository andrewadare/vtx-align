#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

void Fun4All_VTX_ZeroField_Standalone(int nEvents = 0,
                                      const char *inputfile = "/phenix/u/theok/hhj2/PRDF/ZEROFDATA_P00-0000411768-0000.PRDFF",
                                      const char *outputDir = ".",
                                      const char *parfile = "/direct/phenix+u/dcm07e/work/vtx-align/geom/svxPISA-ideal.par",
                                      float beamcenter_x = 0,
                                      float beamcenter_y = 0,
                                      float eastToWest_x = 0,
                                      float eastToWest_y = 0,
                                      float eastToWest_z = 0,
                                      float vtxToCNT_x = 0,
                                      float vtxToCNT_y = 0,
                                      float vtxToCNT_z = 0,
                                      const char *pixel_refmap = "",
                                      const char *pixel_diffmap = "",
                                      const char *pixel_chipmap = "",
                                      const char *strip_deadchannel = "",
                                      const char *strip_deadRCC = ""
                                     )
{
    //Tell root to really crash when something goes wrong not start it's
    //signal handling.
    for (int i = 0; i < kMAXSIGNALS; i++)
    {
        gSystem->IgnoreSignal((ESignals)i);
    }

    gSystem->Exec("/bin/env");

    char ifile[5000];
    strcpy(ifile, inputfile);
    strtok(ifile, "-");
    int runnumber = atoi(strtok(0, "-"));
    int segnumber = atoi(strtok(strtok(0, "-"), "."));
    cout << "run num: " << runnumber << " seg num: " << segnumber << endl;
    //gSystem->Load("libjprof.so");
    //prof *Pr = new prof;

    ///////////////////////////////////////////
    // Load Libraries
    //////////////////////////////////////////
    gSystem->Load("libfun4all.so");
    gSystem->Load("libfun4allfuncs.so");
    gSystem->Load("libcompactCNT.so");
    //  gSystem->Load("librecal");
    gSystem->Load("libSvxDstQA.so");
    //gSystem->Load("libSvxAlignment.so");
    gSystem->Load("libProdEventCutter.so");


    gROOT->ProcessLine(".L OutputManager.C");
    gROOT->ProcessLine(".L rawdatacheck.C");
    gROOT->ProcessLine(".L TrigSelect.C");

    gSystem->ListLibraries();

    SetCvsTag();

    ///////////////////////////////////////////
    // recoConsts setup
    //////////////////////////////////////////
    recoConsts *rc = recoConsts::instance();
    rc->set_FloatFlag("EASTMAXSAG", -0.017);
    rc->set_FloatFlag("WESTMAXSAG", -0.017);

    // set the cvstag from the build to store in the dst
    //SetCvsTag();

    ///////////////////////////////////////////
    // Make the Server
    //////////////////////////////////////////
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);

    ///////////////////////////////////////////
    // limit insert time for calibrations
    //////////////////////////////////////////
    PdbBankManager *bankManager = PdbBankManager::instance();
    //bankManager->SetMaxInsertTime(1390820000) ; // Mon Jan 27 05:53:20 2014

    ///////////////////////////////////////////
    // Make and register the Raw Data Checker
    //////////////////////////////////////////
    RawDataCheck *raw = rawdatacheck();

    ///////////////////////////////////////////
    // Make the Synchronization Object
    ///////////////////////////////////////////
    SubsysReco *sync = new SyncReco();

    //////////////////////////////////////////
    // Central arms
    //////////////////////////////////////////
    HeadReco *head    = new HeadReco();
    head->SetRawDataCheck(raw); // add the rawdatacheck pointer so a list of
    // bad packets is added to the EventHeader
    SubsysReco *trig    = new TrigReco();
    SubsysReco *peve    = new PreviousEventReco();

    //////////////////////////////////
    // Define the triggers
    //////////////////////////////////

    //////////////////////////////////
    // Accounting
    //////////////////////////////////
    SubsysReco *trigacc = new TriggerAccounting();
    SubsysReco *outacc = new OutputAccounting();

    //////////////////////////////////////////
    // Event
    //////////////////////////////////////////
    BbcReco *bbc     = new BbcReco();
    bbc->setBbcVtxError( 0.5 );

    SubsysReco *vtx     = new VtxReco();
    SubsysReco *global  = new GlobalReco();
    SubsysReco *global_central  = new GlobalReco_central();

    SubsysReco *cutter = new ProdEventCutter();
    //cutter->Verbosity(1);

    //////////////////////////////////////////
    // VTX
    //////////////////////////////////////////



    SvxParManager *svxpar = new SvxParManager();
    svxpar->set_BeamCenter(beamcenter_x, beamcenter_y);
    svxpar->set_OffsetVtxToCnt(vtxToCNT_x, vtxToCNT_y, vtxToCNT_z);
    svxpar->set_OffsetEastToWest(eastToWest_x, eastToWest_y, eastToWest_z);
    svxpar->set_ReadGeoParFromFile(1);
    svxpar->set_GeometryFileName(parfile);

    svxpar->set_UseRefDiffPixelMap(1);
    svxpar->set_ReadRefDiffPixelMapFromFile(1);
    svxpar->set_RefDiffPixelMapFiles(pixel_refmap, pixel_diffmap, pixel_chipmap);

    svxpar->set_ReadStripHotDeadFromFile(1);
    svxpar->set_StripHotDeadFileName(strip_deadchannel);
    svxpar->set_StripHotDeadReadoutsFileName(strip_deadRCC);


    SvxDecode *svxdecode = new SvxDecode();
    svxdecode->includePixel(true);
    svxdecode->includeStripixel(true);
    svxdecode->setAdcOffset(24);
    svxdecode->setAdcCutoff(-24);

    SubsysReco *svxhotdead                   = new SvxApplyHotDead();
    SubsysReco *svxrec                       = new SvxReco();
    SvxPriVertexSeedFinder *svxvtxseedfinder = new SvxPriVertexSeedFinder();

    SvxStandAloneReco *svxstandalone         = new SvxStandAloneReco();
    //svxstandalone->setProjectionFlag(false);
    svxstandalone->setZerofieldFlag(true);
    //svxstandalone->setVertexRecoFlag(2);
    svxstandalone->setPPFlag(true);
    //svxstandalone->setWindowScale(20.0);

    SvxPrimVertexFinder *svxprimvtxfinder    = new SvxPrimVertexFinder();
    SubsysReco *svxprimvtxfinder_west = new SvxPrimVertexFinder("SVXPRIMVTXFINDERW", 1);
    SubsysReco *svxprimvtxfinder_east = new SvxPrimVertexFinder("SVXPRIMVTXFINDERE", 2);

    //////////////////////////////////
    // Register SubSystems
    //////////////////////////////////
    se->registerSubsystem(head);
    se->registerSubsystem(sync);
    se->registerSubsystem(trig);
    se->registerSubsystem(trigacc);

    ///////////////////////////////////////////
    /// Trigger dicing now in loaded TrigSelect.C macro
    ///////////////////////////////////////////
    TrigSelect();

    se->registerSubsystem(outacc);

    se->registerSubsystem(peve);
    se->registerSubsystem(bbc);
    se->registerSubsystem(vtx);
    se->registerSubsystem(global);
    se->registerSubsystem(global_central);
    se->registerSubsystem(cutter);

    se->registerSubsystem(svxpar);
    se->registerSubsystem(svxdecode);
    se->registerSubsystem(svxhotdead);
    se->registerSubsystem(svxrec);
    se->registerSubsystem(svxvtxseedfinder);
    se->registerSubsystem(svxstandalone);
    se->registerSubsystem(svxprimvtxfinder);
    se->registerSubsystem(svxprimvtxfinder_east);
    se->registerSubsystem(svxprimvtxfinder_west);

    //////////////////////////////////
    // Output
    //////////////////////////////////
    // trgsel is a global vector from TrigSelect.C which
    // contains the names of the trigger selectors
    // which are used to determine the filenames
    char dstname[100];
    for (int i = 0; i < trgsel.size(); i++)
    {
        //        sprintf(dstname, "v3_CNT_%s", trgsel[i].c_str());
        //        CNT_Compact(runnumber, segnumber, dstname, trgsel[i].c_str());
        //        sprintf(dstname, "v3_DST_EVE_%s", trgsel[i].c_str());
        //        DST_EVE(runnumber, segnumber, dstname, trgsel[i].c_str());
        sprintf(dstname, "DST_SVX_%s", trgsel[i].c_str());
        DST_SVX(runnumber, segnumber, outputDir, dstname, trgsel[i].c_str());
        //      sprintf(dstname, "DST_MPC_%s", trgsel[i].c_str());
        //      DST_MPC(runnumber, segnumber, dstname, trgsel[i].c_str());
    }

    ///////////////////////////////////////////
    // Analyze the Data.
    //////////////////////////////////////////
    gSystem->Exec("ps -o sid,ppid,pid,user,comm,vsize,rssize,time");
    /*
    DeathToMemoryHogs * dt = new DeathToMemoryHogs();
    dt->event_frequency(100); // check every hundred events
    dt->SaveHisto();  // save data in histogram
    se->registerSubsystem(dt);
    */
    Fun4AllInputManager *in = new Fun4AllPrdfInputManager("PRDFin");
    in->fileopen(inputfile);
    se->registerInputManager(in);
    se->run(nEvents);

    se->End();
    PrintTrigSelect();

    int evts = se->PrdfEvents();

    std::cout << "Total Events:  " << evts << std::endl;


    cout << "Successfully Completed Analysis." << endl;
    PHTimeServer::get()->print_stat();

    delete se;

    cout << "Fun4All successfully completed " << endl;
    gSystem->Exit(0);

}
