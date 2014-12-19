analyze_run14_test_prod_output(int nevents = 0,
                                  const char *dstFile = "",
                                  const char *outfile = "test_testvtxproduction.root"
                                 )
{
    gROOT->SetBatch(true);

    gSystem->Load("libfun4all.so");
    gSystem->Load("libfun4allfuncs.so");
    gSystem->Load("libCNT.so");
    gSystem->Load("libcompactCNT.so");
    gSystem->Load("libsvx.so");
    gSystem->Load("libsimreco.so");
    gSystem->Load("librecal.so");


    gSystem->Load("libAnaAlignmentProd.so");
    gSystem->Load("libRun11VTXana.so");

    gSystem->ListLibraries();

    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);
    recoConsts *rc = recoConsts::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);


    //Run14EventTrigger *trig = new Run14EventTrigger();
    //trig->Verbosity(1);
    //trig->set_bbczcut(10);
    //se->registerSubsystem(trig);

    AnaAlignmentProd *svxtest = new AnaAlignmentProd(outfile);
    svxtest->SetZeroField(true);
    se->registerSubsystem(svxtest);


    Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
    in1->Verbosity(0);
    se->registerInputManager(in1);

    cout << "Analysis started " << endl;
    gSystem->Exec("date");
    se->fileopen("DSTin1", dstFile);
    se->run(nevents);

    se->End();

    cout << "Analysis finished " << endl;



}
