// These routines need to stay together.
// The raw data check now puts bad events into a DB and
// we do not want to declare an event bad for muon reconstruction
// only
// This rarely happens anyway - the only way to get this is if
// a granule drops out of sync temporarily and I have not seen
// that in Run4
//
// chp

RawDataCheck *rawdatacheck()
{
// just in case this lib gets disentangled
  gSystem->Load("libRawDataCheck.so"); 
  Fun4AllServer *se = Fun4AllServer::instance();
  RawDataCheck *chk = RawDataCheck::instance();

  //  Create Check Modules
  GranuleCheck *Acc = new CheckAcc();
  GranuleCheck *Bbc = new CheckBbc();
  GranuleCheck *DchEast = new CheckDch("EAST");
  GranuleCheck *DchWest = new CheckDch("WEST");
  GranuleCheck *EmcEastBottom = new CheckEmc("EAST", "BOTTOM");
  GranuleCheck *EmcEastTop = new CheckEmc("EAST", "TOP");
  GranuleCheck *EmcWestBottom = new CheckEmc("WEST", "BOTTOM");
  GranuleCheck *EmcWestTop = new CheckEmc("WEST", "TOP");
  GranuleCheck *ErtEast = new CheckErt("EAST");
  GranuleCheck *ErtWest = new CheckErt("WEST");
  GranuleCheck *Mpc = new CheckFcal("MPC");
  GranuleCheck *MuidNorth = new CheckMuid("NORTH");
  GranuleCheck *MuidSouth = new CheckMuid("SOUTH");
  GranuleCheck *MutrNorth = new CheckMutr("NORTH");
  GranuleCheck *MutrSouth = new CheckMutr("SOUTH");
  GranuleCheck *PadEast = new CheckPad("EAST");
  GranuleCheck *PadWest = new CheckPad("WEST");
  GranuleCheck *RichEast = new CheckRich("EAST");
  GranuleCheck *RichWest = new CheckRich("WEST");
  GranuleCheck *TofEast = new CheckTof("EAST");
  GranuleCheck *TofWest = new CheckTof("WEST");
  GranuleCheck *Vtxp = new CheckVtxp("VTXP");
  GranuleCheck *Vtxs = new CheckVtxs("VTXS");
  GranuleCheck *Zdc = new CheckZdc();

  //  OK, register them all...
  chk->registerGranule(Acc);
  chk->registerGranule(Bbc);
  chk->registerGranule(DchEast);
  chk->registerGranule(DchWest);
  chk->registerGranule(EmcEastBottom);
  chk->registerGranule(EmcEastTop);
  chk->registerGranule(EmcWestBottom);
  chk->registerGranule(EmcWestTop);
  chk->registerGranule(ErtEast);
  chk->registerGranule(ErtWest);
  chk->registerGranule(Mpc);
  chk->registerGranule(MuidNorth);
  chk->registerGranule(MuidSouth);
  chk->registerGranule(MutrNorth);
  chk->registerGranule(MutrSouth);
  chk->registerGranule(PadEast);
  chk->registerGranule(PadWest);
  chk->registerGranule(RichEast);
  chk->registerGranule(RichWest);
  chk->registerGranule(TofEast);
  chk->registerGranule(TofWest);
  chk->registerGranule(Vtxp);
  chk->registerGranule(Vtxs);
  chk->registerGranule(Zdc);

  // register everything with Fun4All
  se->registerSubsystem(chk);
  //  chk->UpdateDB(0); // prevent db updates in test mode, 
  chk->UpdateDB(1); // prevent db updates in test mode, 1 to write bad events into database
  return chk;
}

// this one is looking into the rebuild.info file
// the cvs tag is extracted from the lines
//  CVS tag:
// -r CVSTAG
// So here it first looks for the "tag:" keyword,
// then checks for the -r in the next line (if it is from the
// the nightly rebuild, this keyword is missing)
// and uses the next string as cvs tag and sets the rc flag accordingly
// if no cvs tag is found this way, the default is "none"


// build the filename from OFFLINE_MAIN
void SetCvsTag()
{
  char filename[200];
  char *offline = getenv("OFFLINE_MAIN");
  sprintf(filename,"%s/rebuild.info",offline);
  SetCvsTagFromFile(filename);
  return ;
}

void SetCvsTagFromFile(const char *filename)
{
  gSystem->Load("libphflag.so");
  recoConsts *rc = recoConsts::instance();
  FILE *rebuildinfo = fopen(filename, "r");
  if (! rebuildinfo)
    {
      rc->set_CharFlag("CVSTAG", "newbuild");
    }
  else
    {
      rc->set_CharFlag("CVSTAG", "none");
      char line[101];
      int foundtag = 0;
      int foundr = 0;
      recoConsts *rc = recoConsts::instance();
      while (fscanf(rebuildinfo, "%100s", line) != EOF)
        {
          if (foundr)
            {
              string cvstag = line;
              rc->set_CharFlag("CVSTAG", cvstag.c_str());
              break;
            }
          if (foundtag)
            {
              if (strcmp(line, "-r") == 0)
                {
                  foundr = 1;
                  foundtag = 0;
                }
            }
          if (strcmp(line, "tag:") == 0)
            {
              foundtag = 1;
            }
        }
      fclose(rebuildinfo);
    }
  cout << "Setting CVSTAG to " << rc->get_CharFlag("CVSTAG") << endl;
  return ;
}

