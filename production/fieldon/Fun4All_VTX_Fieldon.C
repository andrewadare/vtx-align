#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

void Fun4All_VTX_Fieldon(int nEvents = 0,
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
                         const char *strip_deadRCC = "",
                         float rotphcntphi_east = 0,
                         float rotphcntphi_west = 0,
                         float rotphcntthe_east = 0,
                         float rotphcntthe_west = 0
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
    //rc->set_IntFlag("PADVTX_PC3_MULT_CUT",10);  // padvtx cut
    // No MPC in CuAu
    //  rc->set_IntFlag("MPC_RECO_MODE",0x3);
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
    // Central arms
    //////////////////////////////////////////
    BbcReco *bbc     = new BbcReco();
    bbc->setBbcVtxError( 0.5 );

    SubsysReco *ert     = new ErtReco();
    SubsysReco *zdc     = new ZdcReco();
    SubsysReco *t0      = new T0Reco();
    SubsysReco *pad     = new PadReco();
    SubsysReco *vtx     = new VtxReco();
    SubsysReco *emc     = new EmcReco3();
    SubsysReco *emcres  = new EmcClusterContainerResurrector();
    //SubsysReco *padvtx  = new PadVtxReco();
    SubsysReco *tof     = new TofReco();
    SubsysReco *aero    = new AccReco();
    SubsysReco *dch     = new DchReco();
    SubsysReco *crk     = new CrkReco();
    CglReco *cgl     = new CglReco();
    SubsysReco *aerocl  = new AccclusterReco();
    SubsysReco *ring    = new RingReco();
    SubsysReco *tofw    = new TofwReco();
    SubsysReco *global  = new GlobalReco();
    SubsysReco *global_central  = new GlobalReco_central();
    //  SubsysReco *mpc     = new MpcReco();

    SubsysReco *cutter = new ProdEventCutter();



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
    svxstandalone->setVertexRecoFlag(2);

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
    se->registerSubsystem(zdc);
    se->registerSubsystem(t0);
    se->registerSubsystem(pad);
    se->registerSubsystem(vtx);
    se->registerSubsystem(emc);
    se->registerSubsystem(emcres);
    //se->registerSubsystem(padvtx);
    se->registerSubsystem(ert);
    se->registerSubsystem(tof);
    se->registerSubsystem(aero);
    se->registerSubsystem(dch);
    se->registerSubsystem(crk);
    se->registerSubsystem(tofw);

    se->registerSubsystem(cgl);

    se->registerSubsystem(aerocl);
    se->registerSubsystem(ring);
    //  se->registerSubsystem(mpc);

    se->registerSubsystem(global);
    se->registerSubsystem(global_central);
    //se->registerSubsystem(cutter);


    
    se->registerSubsystem(svxpar);
    se->registerSubsystem(svxdecode);
    se->registerSubsystem(svxhotdead);
    se->registerSubsystem(svxrec);
    se->registerSubsystem(svxvtxseedfinder);
    se->registerSubsystem(svxstandalone);
    se->registerSubsystem(svxprimvtxfinder);
    se->registerSubsystem(svxprimvtxfinder_east);
    se->registerSubsystem(svxprimvtxfinder_west);

    //=========================================
    // These fill the compactCNT storage nodes
    //=========================================
    SubsysReco *fillprojections = new FillTrackProjections();
    SubsysReco *filllineprojections = new FillTrackLineProjections();
    SubsysReco *fillpl = new FillTrackPathLengths();
    SubsysReco *filltrkhits = new FillTrackHits();
    SubsysReco *fillpadhits = new FillPadHits();
    SubsysReco *filldchits = new FillDchHits();
    SubsysReco *filltofehits = new FillTofeHits();
    SubsysReco *filltofwhits = new FillTofwHits();
    SubsysReco *fillcrkhits = new FillCrkHits();
    SubsysReco *fillacchits = new FillAccHits();
    SubsysReco *fillemchits = new FillEmcHits();

    se->registerSubsystem(fillprojections);
    se->registerSubsystem(filllineprojections);
    se->registerSubsystem(fillpl);
    se->registerSubsystem(filltrkhits);
    se->registerSubsystem(filldchits);
    se->registerSubsystem(fillpadhits);
    se->registerSubsystem(filltofehits);
    se->registerSubsystem(filltofwhits);
    se->registerSubsystem(fillcrkhits);
    se->registerSubsystem(fillacchits);

    // This one requires that EmcClusterContainer is already on the node tree
    se->registerSubsystem(fillemchits);

    //==============================================
    // These modules read the compactCNT nodes
    // and create hits objects for each subsystem
    //================================================
    se->registerSubsystem(new RecoverTrackProjections());
    se->registerSubsystem(new RecoverTrackLineProjections());
    se->registerSubsystem(new RecoverTrackPathLengths());
    se->registerSubsystem(new RecoverTrackHits());
    se->registerSubsystem(new RecoverDchHits());
    se->registerSubsystem(new RecoverPadHits());
    se->registerSubsystem(new RecoverTofeHits());
    se->registerSubsystem(new RecoverTofwHits());
    se->registerSubsystem(new RecoverCrkHits());
    se->registerSubsystem(new RecoverAccHits());
    se->registerSubsystem(new RecoverEmcHits());

    //========================
    // Creates PHCentralTrack
    //========================

    se->registerSubsystem(new CreateCNT());

    //=================================================
    // These modules re-associate hits with tracks and
    // fill the PHCentralTrack fields
    //==================================================
    se->registerSubsystem(new FillCNT_TrackProjections());
    se->registerSubsystem(new FillCNT_TrackPathLengths());
    se->registerSubsystem(new FillCNT_TrackHits());
    se->registerSubsystem(new FillCNT_DchHits());
    se->registerSubsystem(new FillCNT_TofeHits());
    se->registerSubsystem(new FillCNT_TofwHits());
    se->registerSubsystem(new FillCNT_PadHits());
    se->registerSubsystem(new FillCNT_CrkHits());
    se->registerSubsystem(new FillCNT_AccHits());
    // This one needs EmcClusterContainer also
    se->registerSubsystem(new FillCNT_EmcHits());

    /// SvxCentralTrackReco should be called after PHCentralTrack is reconstructed.
    SvxCentralTrackReco *svxcentraltrack = new SvxCentralTrackReco();
    //svxcentraltrack->setSearchWindowFlag(2);
    svxcentraltrack->setShiftPHCentralTrack(true);
    svxcentraltrack->setPHCentralTrackPhiThetaShift(0, rotphcntphi_east, rotphcntthe_east);
    svxcentraltrack->setPHCentralTrackPhiThetaShift(1, rotphcntphi_west, rotphcntthe_west);
    se->registerSubsystem(svxcentraltrack);


    //In this macro we want ALL clusters
    //SvxSelectClusters* svxselect = new SvxSelectClusters();
    //se->registerSubsystem(svxselect);

    // svx compactCNTs
    //FillSvxHits *fillsvxhits = new FillSvxHits();
    //fillsvxhits->Verbosity(0);
    //fillsvxhits->setSaveOnlySelectedHits(true);
    //se->registerSubsystem(fillsvxhits);

    // Reaction Plane
    //SubsysReco *svxrp   = new SvxRpSumXYReco();
    //se->registerSubsystem(svxrp);

    //////////////////////////////////
    // Central arm output
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

    SubsysReco *svxproductionqa  = new SvxQAForProduction();
    se->registerSubsystem(svxproductionqa);

    SubsysReco *svxalignmentqa  = new SvxAlignment_QA();
    se->registerSubsystem(svxalignmentqa);

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
    SVXQA_IOManager(runnumber, segnumber, outputDir);
    //  se->dumpHistos(histofile);
    PrintTrigSelect();

    int evts = se->PrdfEvents();

    std::cout << "Total Events:  " << evts << std::endl;


    cout << "Successfully Completed Analysis." << endl;
    PHTimeServer::get()->print_stat();

    delete se;
    /*
        char sql[1000];
        sprintf(sql, "echo \"update %s set nevents=%d,time='now' where runnumber=%d and sequence=%d and prodtag = \'%s\';\" | isql phnxreco_odbc", gSystem->Getenv("PRODDBTABLE"), evts, runnumber, segnumber, gSystem->Getenv("PRODTAG"));
        std::cout << "Update DB statement: " << sql << std::endl;
        gSystem->Exec(sql);

        gSystem->Exec("ps -o sid,ppid,pid,user,comm,vsize,rssize,time");
    */
    cout << "Fun4All successfully completed " << endl;
    gSystem->Exit(0);

}
