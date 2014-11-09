analyze_run14_test_prod_output(int nevents = 0,
                               const char *dstFile = "",
                               const char *outfile = "test_testvtxproduction.root"
                                   //const char *treefile = "test_svxcnttree.root",
                                   //const char *histfile = "test_svxeventhist.root",
                                   //const char *eventfile = "test_svxeventtree.root"
                              )
{
  gROOT->SetBatch(true);

  //const char *infile1 = "/direct/phenix+prod05/phnxreco/run14_online/run14_online_ca/run_0000402000_0000403000/CNT/CNT_run14_online_ca-0000402782-0000.root";
  //const char *infile2  = "/direct/phenix+prod05/phnxreco/run14_online/run14_online_ca/run_0000402000_0000403000/DST_SVX/DST_SVX_run14_online_ca-0000402782-0000.root";
  //const char *infile3  = "/direct/phenix+prod05/phnxreco/run14_online/run14_online_ca/run_0000402000_0000403000/DST_EVE/DST_EVE_run14_online_ca-0000402782-0000.root";


  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libCNT.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libsvx.so");
  gSystem->Load("libsimreco.so");
  gSystem->Load("librecal.so");

  gSystem->Load("libRun14EventTrigger.so");
  //gSystem->Load("libRun11AuAuSvxCentralTrackTree.so");
  //gSystem->Load("libRun11AuAuSvxEventHists.so");
  //gSystem->Load("libRun11AuAuSvxEventTree.so");
  gSystem->Load("libTestVTXProduction.so");
  gSystem->Load("libRun11VTXana.so");


  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  recoConsts *rc = recoConsts::instance();

  MasterRecalibratorManager *mr = new MasterRecalibratorManager();
  se->registerSubsystem(mr);


  //Run14EventTrigger *trig = new Run14EventTrigger();
  //trig->Verbosity(1);
  //trig->set_bbczcut(10);
  //se->registerSubsystem(trig);

  SubsysReco *svxtest = new TestVTXProduction(outfile);
  se->registerSubsystem(svxtest);


  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("DSTin1", "DST");
  in1->Verbosity(0);
  se->registerInputManager(in1);

  cout << "Analysis started " << endl;
  gSystem->Exec("date");
  se->fileopen("DSTin1", dstFile);
  se->run(nevents);

  se->End();

  //Fun4AllHistoManager *hm;
  //hm = se->getHistoManager("Run11AuAuSvxEventHists");
  //hm->dumpHistos(histfile);

  cout << "Analysis finished " << endl;



}
