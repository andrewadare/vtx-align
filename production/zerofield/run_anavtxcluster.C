
void run_anavtxcluster(int nevents = 0,
		       int segNumber = 0,
                       const char *infile = "",
                       const char *outfile = "")
{
    gSystem->Load("libfun4all.so");
    gSystem->Load("libfun4allfuncs.so");
    gSystem->Load("libCNT.so");
    gSystem->Load("libcompactCNT.so");
    gSystem->Load("libsvx.so");
    gSystem->Load("libsvxeval.so");
    gSystem->Load("libsimreco.so");
    gSystem->Load("librecal.so");
    gSystem->Load("libcteval");

    gSystem->Load("libAnaVTXCluster.so");

    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(0);
    recoConsts *rc = recoConsts::instance();

    MasterRecalibratorManager *mr = new MasterRecalibratorManager();
    se->registerSubsystem(mr);

    SubsysReco *ana = new AnaVTXCluster(outfile);
    ana->SetEventOffset(segNumber);
    se->registerSubsystem(ana);

    // Input manager(s)
    Fun4AllInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTin1", "DST");
    in1->Verbosity(0);
    se->registerInputManager(in1);

    cout << "Analysis started " << endl;
    gSystem->Exec("date");

    se->fileopen("DSTin1",infile);
    se->run(nevents);
    se->fileclose("DSTin1");

    se->End();

    cout << "Analysis finished " << endl;

}
