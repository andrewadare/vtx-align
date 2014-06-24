/*
   rootlogon.C
   This script is automatically invoked whenever ROOT is started.
   Add session-level configurations as needed.
*/

void rootlogon()
{
  Printf("Starting ROOT version %s.", gROOT->GetVersion());
  Printf("Running %s/rootlogon.C on %s.",
         gSystem->Getenv("PWD"),
         gSystem->HostName());


  bool atrcf = false;

  if (atrcf)
    gSystem->Load("libsvxgeo"); // src: offline/packages/svxgeo
  else
  {
    // For dev work / remote use: build shared libs with ACLiC.
    // Set this to your source path
    const char *geo = "/Users/adare/phenix/svxgeo";

    // There are 2 include paths: one for ACLiC, and one for CINT or CLING.
    // Type .include (ROOT 5) or .I (ROOT 6) at the ROOT REPL to list.
    gSystem->AddIncludePath(Form("-I%s ", geo));
#ifdef __CINT__
    gROOT->ProcessLine(Form(".include %s", geo));
#endif
#ifdef __CLING__
    gROOT->ProcessLine(Form(".I %s", geo));
#endif
    gROOT->LoadMacro(Form("%s/SvxTGeo.C+", geo));
    gROOT->LoadMacro(Form("%s/SvxGeoTrack.C+", geo));
    gROOT->LoadMacro(Form("%s/SvxProj.C+", geo));
  }

  // Build library from Mille source
  const char *mpd = "mille";
  gSystem->AddIncludePath(Form("-I%s ", mpd));
#ifdef __CINT__
  gROOT->ProcessLine(Form(".include %s", mpd));
#endif
#ifdef __CLING__
  gROOT->ProcessLine(Form(".I %s", mpd));
#endif

  gROOT->LoadMacro(Form("%s/Mille.cc+", mpd));

  return;
}
