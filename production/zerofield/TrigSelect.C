#include <vector>
vector<string> trgsel;
vector<string>::const_iterator striter;

void TrigSelect()
{
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  Fun4AllServer *se = Fun4AllServer::instance();

  trgselname = "MB";
  TrigSelect *minbias = new TrigSelect(trgselname);
  trgsel.push_back(trgselname);
  minbias->AddTrigger( "BBCLL1(>1 tubes)" );
  minbias->AddTrigger( "BBCLL1(>1 tubes) narrowvtx" );
  minbias->AddTrigger( "BBCLL1(>1 tubes) narrowvtx CopyA" );
  minbias->AddTrigger( "BBCLL1(>1 tubes) narrowvtx CopyB" );
  minbias->AddTrigger( "BBCLL1(>1 tubes) novertex" );


  TrigSelect *ppg      = new TrigSelect("PPG");      // This will veto ppg triggers
  ppg->AddVetoTrigger("PPG(Laser)");
  ppg->AddVetoTrigger("PPG(Pedestal)");
  ppg->AddVetoTrigger("PPG(Test Pulse)");
  ppg->SetReturnCode("ABORT");

// This will take all triggers from the previous trigselect modules
// and will only save events which were triggered by other triggers
// It will save the event even if it was also triggered by a previously
// selected trigger
  /*
  trgselname = "OT";
  TrigSelect *others = new TrigSelect(trgselname);
  trgsel.push_back(trgselname);
  others->AddNoSaveTrigger(erttrig);
  others->AddNoSaveTrigger(minbias);
  others->AddNoSaveTrigger(mutrig);
  */

  se->registerSubsystem(ppg);
  //se->registerSubsystem(erttrig);
  //se->registerSubsystem(mutrig);
  se->registerSubsystem(minbias);
  //se->registerSubsystem(others);
}

void PrintTrigSelect()
{

  Fun4AllServer *se = Fun4AllServer::instance();
  for (int i = 0; i < trgsel.size(); i++)
    {
      TrigSelect *trg = (TrigSelect *) se->getSubsysReco(trgsel[i].c_str());
      trg->Print();
      cout << endl;
    }
}
