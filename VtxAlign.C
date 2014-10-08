#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "ParameterDefs.h"
#include "MilleFunctions.h"
#include "ConstraintBuilder.h"
#include "VtxIO.h"

using namespace std;

// Globals
const double BField = 0.0;
const bool drawResults = true;
const bool useVtxTracks = true;
const bool useCntTracks = false;
const double regFactor = 1.0;

// Global "presigma" for regularization. Smaller = stronger reg (but 0=none).
// Written to steering file.
// Entries in constraints file(s) take precedence over this.
const double defaultPreSigma = 0.01;

const char *pedeSteerFile = "pede-steer.txt";
const char *ladderConstFile = "ladder_constraints.txt";
const char *hluConstFile = "halflayer_constraints.txt";
string pedeBinFileStd = "standalone.bin";
string pedeBinFileCnt = "svxcnttrks.bin";

void EventLoop(string binfile, geoEvents &events, vecs &sgpars, vecs &zgpars,
               TGraphErrors *bc = 0, TString opt = "");
void CorrectFromFile(const char *filename,
                     SvxTGeo *tgeo,
                     geoEvents &vtxevents,
                     geoEvents &cntevents);

void VtxAlign(int run = 411768,    // Run number of PRDF segment(s)
              int prod = 0,        // Production step. Starts at 0.
              int subiter = 1,     // Geometry update step. Starts at 0.
              TString alignMode = "ladder") // "ladder","halflayer" (+"sim")
{
  // No point in continuing if Millepede II is not installed...
  if (TString(gSystem->GetFromPipe("which pede")).IsNull())
  {
    Printf("\"which pede\" returns nothing. Exiting.");
    gSystem->Exit(-1);
  }

  // Inputs:
  TString rootFileIn = Form("rootfiles/%d-%d-%d.root", run, prod, subiter);
  TString pisaFileIn = Form("geom/%d-%d-%d.par", run, prod, subiter);
  // TString bcFileIn   = Form("rootfiles/bc-%d-%d-%d.root", run, prod, subiter);
  TString bcFileIn   = Form("rootfiles/bc-%d.root", run);
  TFile *bcf         = new TFile(bcFileIn.Data(), "read");
  TGraphErrors *gbc  = (TGraphErrors *) bcf->Get("gbc");

  // Outputs:
  TString rootFileOut = Form("rootfiles/%d-%d-%d.root", run, prod, subiter + 1);
  TString pisaFileOut = Form("geom/%d-%d-%d.par", run, prod, subiter + 1);

  // Assign free global parameter coords for ds=r*dphi residuals here
  // Set includes "s", "x", "y", "r".
  // Also assign presigma list for these coordinates (trumps defaultPreSigma).

  // vecs sgpars {"s", "x", "y"};
  // vecd sgpresigma {0,0,0};

  vecs sgpars {"x", "y"};
  vecd sgpresigma {0,0};

  // vecs sgpars {"s"};
  // vecd sgpresigma {0};

  // Assign free global parameter coords for dz residuals here
  // Set includes "x", "y", "z", "r".
  vecs zgpars {"z"};
  vecd zgpresigma {0};
  // vecs zgpars {"z", "r"};
  // vecd zgpresigma {0, 1e-4};

  TFile *inFile = new TFile(rootFileIn.Data(), "read");
  assert(inFile);

  TTree *svxseg = (TTree *)inFile->Get("vtxtrks");
  TTree *svxcnt;
  assert(svxseg);

  if (useCntTracks > 0)
  {
    svxcnt = (TTree *)inFile->Get("cnttrks");
    assert(svxcnt);
  }

  SvxTGeo *tgeo = VTXModel(pisaFileIn.Data());

  // Write constraints to text file
  vecs constfiles;
  vecs binfiles;
  if (alignMode.Contains("halflayer") && useCntTracks == false)
  {
    constfiles.push_back(hluConstFile);
    WriteHLConstraints(hluConstFile, sgpars, zgpars,
                       sgpresigma, zgpresigma, tgeo, alignMode);
  }
  else if (alignMode.Contains("ladder"))
  {
    constfiles.push_back(ladderConstFile);
    WriteLadderConstraints(ladderConstFile, sgpars, zgpars,
                           sgpresigma, zgpresigma, tgeo, alignMode);
  }

  // Write Millepede steering file
  if (useVtxTracks)
    binfiles.push_back(pedeBinFileStd);
  if (useCntTracks)
    binfiles.push_back(pedeBinFileCnt);
  WriteSteerFile(pedeSteerFile, binfiles, constfiles, regFactor, defaultPreSigma);

  geoEvents vtxevents; // Events containing SVX standalone tracks
  geoEvents cntevents; // Events containing SvxCentralTracks

  if (useVtxTracks)
  {
    GetEventsFromTree(svxseg, tgeo, vtxevents, -1);

    if (gbc)
    {
      Info("", "Using beamcenter from %s", bcf->GetName());
      Info("", " E: (%.3f, %.3f)", gbc->GetX()[0], gbc->GetY()[0]);
      Info("", " W: (%.3f, %.3f)", gbc->GetX()[1], gbc->GetY()[1]);
    }
    // FitTracks(vtxevents, gbc, "fit_to_bc,find_vertex,calc_dca");
    EventLoop(pedeBinFileStd, vtxevents, sgpars, zgpars, gbc, alignMode);
  }
  if (useCntTracks)
  {
    GetEventsFromTree(svxcnt, tgeo, cntevents);
    EventLoop(pedeBinFileCnt, cntevents, sgpars, zgpars, 0,
              alignMode + ", ext");
  }

  // Shell out to pede executable
  gSystem->Exec(Form("pede %s", pedeSteerFile));

  // Copy geometry adjustments (with their errors) to logs/ directory.
  gSystem->Exec(Form("cp millepede.res logs/dp-%d-%d-%d.txt",
                     run, prod, subiter));

  // 1. Change geometry: apply position corrections from pede output file.
  // 2. Update stored global hit (x,y,z) positions in geoEvents vectors.
  // 3. Update residuals on all hits:
  //  - CNT track residuals are updated in CorrectFromFile().
  //  - VTX track residuals are updated by refitting.
  CorrectFromFile("millepede.res", tgeo, vtxevents, cntevents);
  Printf("Refitting in post-alignment geometry to update residuals.");
  FitTracks(vtxevents, gbc, "fit_to_bc, find_vertex, calc_dca");

  cout << "Filling output tree(s)..." << flush;
  TFile *outFile = new TFile(rootFileOut.Data(), "recreate");
  TTree *vtxtrks = CreateTree("vtxtrks");
  TTree *cnttrks = CreateTree("cnttrks");
  FillTree(vtxevents, vtxtrks);
  FillTree(cntevents, cnttrks);
  Printf("done.");

  Printf("Writing %s", pisaFileOut.Data());
  tgeo->WriteParFile(pisaFileOut.Data());

  Printf("Writing %s", rootFileOut.Data());
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);

  Printf("Done!");

  if (drawResults)
  {
    Printf("\n\nDrawing results...\n\n");
    if (useVtxTracks)
    {
      const char *macro = Form("DrawResults.C(%d,%d,%d,%d,%d,\"vtxtrks\")",
                               run, prod, subiter, prod, subiter + 1);
      gROOT->Macro(macro);
    }
    if (useCntTracks)
    {
      const char *macro = Form("DrawResults.C(%d,%d,%d,%d,%d,\"cnttrks\")",
                               run, prod, subiter, prod, subiter + 1);
      gROOT->Macro(macro);
    }
  }
  return;
}

void
EventLoop(string binfile, geoEvents &events, vecs &sgpars, vecs &zgpars,
          TGraphErrors *bc, TString opt)
{
  // Call Mille::Mille() in a loop. See MilleFunctions.h
  // Options:
  // - "ext": Assume residuals were computed from external information.
  // - "halflayer": Align half-layer positions instead of ladders.

  Printf("Running EventLoop(). option(s): [%s]", opt.Data());
  if (opt.Contains("ext"))
    Printf(" - Using external tracks.");
  else
    Printf(" - Using VTX standalone tracks.");
  printf(" - Coordinate(s) available for r*phi residual minimization: ( ");
  for (unsigned int ic = 0; ic < sgpars.size(); ++ic)
    cout << sgpars[ic] << " ";
  cout << ")" << endl;
  printf(" - Coordinate(s) available for z residual minimization: ( ");
  for (unsigned int ic = 0; ic < zgpars.size(); ++ic)
    cout << zgpars[ic] << " ";
  cout << ")" << endl;
  if (opt.Contains("ladder"))
    Printf(" - ladder mode");
  if (opt.Contains("halflayer"))
    Printf(" - halflayer mode");
  if (opt.Contains("sim"))
    Printf(" - chi^2 tuned for simulated data. No ladders masked.");
  else
    Printf(" - chi^2 tuned for real data. Bad ladders will be masked.");

  // If asBinary is false, write a text file instead of binary file.
  // For debugging only - text file is not readable by pede.
  bool asBinary = true;
  Mille m(binfile.c_str(), true);

  for (unsigned int ev = 0; ev < events.size(); ev++)
    for (unsigned int t = 0; t < events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (opt.Contains("ext"))
        MilleCnt(m, trk, sgpars, zgpars, bc, opt);
      else
        MilleVtx(m, trk, sgpars, zgpars, bc, opt);
    }

  Printf("Writing to %s...", binfile.c_str());
  // Mille object must go out of scope for output file to close properly.
  return;
}

void
CorrectFromFile(const char *filename,
                SvxTGeo *tgeo,
                geoEvents &vtxevents,
                geoEvents &cntevents)
{
  // Retrieve Millepede's corrections to the global parameters
  // Key is label, value is correction.
  map<int, double> mpc;
  GetCorrections(filename, mpc);

  // Ladder position corrections
  for (int i = 0; i < tgeo->GetNLayers(); i++)
    for (int j = 0; j < tgeo->GetNLadders(i); j++)
    {
      int l = -1;

      // Ladder translation corrections
      l = Label(i, j, "x");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, mpc[l], 0., 0.);
      l = Label(i, j, "y");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, 0., mpc[l], 0.);
      l = Label(i, j, "z");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, 0. , 0., mpc[l]);

      // Ladder Phi correction from ds
      l = Label(i, j, "s");
      if (mpc.find(l) != mpc.end())
        tgeo->RotateLadderRPhi(i, j, mpc[l]);

      // Ladder Radial correction
      l = Label(i, j, "r");
      if (mpc.find(l) != mpc.end())
        tgeo->MoveLadderRadially(i, j, mpc[l]);
    }

  // Half-layer position corrections
  for (int i = 0; i < tgeo->GetNLayers(); i++)
    for (int j = 0; j < 2; j++) // Arm. 0 = E; 1 = W.
    {
      int l = -1;
      // Half-layer translation corrections
      l = HalfLayerLabel(i, j, "x");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateHalfLayer(i, j, mpc[l], 0., 0.);
      l = HalfLayerLabel(i, j, "y");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateHalfLayer(i, j, 0., mpc[l], 0.);
      l = HalfLayerLabel(i, j, "z");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateHalfLayer(i, j, 0. , 0., mpc[l]);

      // Half-ladder phi correction from ds
      // s is converted to phi which is used for the half-layer rotation
      l = HalfLayerLabel(i, j, "s");
      if (mpc.find(l) != mpc.end())
        tgeo->RotateHalfLayerRPhi(i, j, mpc[l]);
    }

  for (unsigned int ev = 0; ev < vtxevents.size(); ev++)
    for (unsigned int t = 0; t < vtxevents[ev].size(); t++)
      vtxevents[ev][t].UpdateHits();
  for (unsigned int ev = 0; ev < cntevents.size(); ev++)
    for (unsigned int t = 0; t < cntevents[ev].size(); t++)
    {
      for (int ihit = 0; ihit < cntevents[ev][t].nhits; ihit++)
        cntevents[ev][t].hits[ihit].UpdateResiduals();

      cntevents[ev][t].UpdateHits();
    }
  return;
}
