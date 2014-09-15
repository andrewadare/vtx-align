
#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "ParameterDefs.h"
#include "MilleFunctions.h"
#include "ConstraintBuilder.h"
#include "VtxIO.h"
#include "VtxVis.h"

using namespace std;

// Globals
const double BField = 0.0;
const int nStdTracks = 999999;
const int nCntTracks = 0;
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

void VtxAlign(int run = 123456,    // Run number of PRDF segment(s)
              int prod = 0,        // Production step. Starts at 0.
              int subiter = 0,     // Geometry update step. Starts at 0.
              TString alignMode = "halflayer") // "ladder" or "halflayer"
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
  TString bcFileIn   = Form("rootfiles/bc-%d-%d-%d.root", run, prod, subiter);

  // Outputs:
  TString rootFileOut = Form("rootfiles/%d-%d-%d.root", run, prod, subiter + 1);
  TString pisaFileOut = Form("geom/%d-%d-%d.par", run, prod, subiter + 1);

  // Select global parameters and derivatives for s residuals
  // Mixing rotation (s) and translation (x,y) is not recommended. They do not
  // commute and should be handled as mutually exclusive.
  vecs sgpars;
  if (0) sgpars.push_back("s");
  if (1) sgpars.push_back("x");
  if (1) sgpars.push_back("y");
  if (0) sgpars.push_back("r");

  vecs zgpars;
  if (1) zgpars.push_back("x");
  if (1) zgpars.push_back("y");
  if (1) zgpars.push_back("z");
  if (0) zgpars.push_back("r");

  TFile *inFile = new TFile(rootFileIn.Data(), "read");
  assert(inFile);

  TNtuple *svxseg = (TNtuple *)inFile->Get("vtxhits");
  TNtuple *svxcnt = (TNtuple *)inFile->Get("cnt_clusntuple");
  assert(svxseg);
  if (nCntTracks>0)
    assert(svxcnt);

  SvxTGeo *tgeo = VTXModel(pisaFileIn.Data());

  // Write constraints to text file
  vecs constfiles;
  vecs binfiles;
  if (alignMode == "halflayer")
  {
    constfiles.push_back(hluConstFile);
    WriteHLConstraints(hluConstFile, sgpars, zgpars, tgeo);
  }
  else if (alignMode == "ladder")
  {
    constfiles.push_back(ladderConstFile);
    WriteLadderConstraints(ladderConstFile, sgpars, zgpars, tgeo);
  }

  // Write Millepede steering file
  if (nStdTracks > 0)
    binfiles.push_back(pedeBinFileStd);
  if (nCntTracks > 0)
    binfiles.push_back(pedeBinFileCnt);
  WriteSteerFile(pedeSteerFile, binfiles, constfiles);

  // Original ladder positions
  vecd x0; vecd y0; vecd z0;
  GetLadderXYZ(tgeo, x0, y0, z0);

  geoEvents vtxevents; // Events containing SVX standalone tracks
  geoEvents cntevents; // Events containing SvxCentralTracks

  if (nStdTracks > 0)
  {
    GetEventsFromTree(svxseg, tgeo, vtxevents, -1);

    TFile *bcf = new TFile(bcFileIn.Data(), "read");
    TGraphErrors *gbc = (TGraphErrors *) bcf->Get("gbc");

    // gbc->SetPointError(0, gbc->GetEX()[0], gbc->GetEY()[0]);
    // gbc->SetPointError(1, gbc->GetEX()[1], gbc->GetEY()[1]);
    // FitTracks(vtxevents, gbc);// Uses beam center. Doesn't work as well.

    FitTracks(vtxevents);
    EventLoop(pedeBinFileStd, vtxevents, sgpars, zgpars, gbc, alignMode);
  }
  if (nCntTracks > 0)
  {
    GetEventsFromTree(svxcnt, tgeo, cntevents);
    EventLoop(pedeBinFileCnt, cntevents, sgpars, zgpars);
  }

  // Shell out to pede executable
  gSystem->Exec(Form("pede %s", pedeSteerFile));

  // Copy geometry adjustments (with their errors) to logs/ directory.
  gSystem->Exec(Form("cp millepede.res logs/dp-%d-%d-%d.txt",
                     run, prod, subiter));

  // Apply position corrections from pede output file.
  // Also update stored global positions in geoEvents vectors.
  CorrectFromFile("millepede.res", tgeo, vtxevents, cntevents);

  // Record positions after alignment
  vecd x1; vecd y1; vecd z1;
  GetLadderXYZ(tgeo, x1, y1, z1);

  Printf("Refitting in post-alignment geometry to get updated residuals.");
  FitTracks(vtxevents);
  cout << "Filling output tree(s)..." << flush;
  TFile *outFile = new TFile(rootFileOut.Data(), "recreate");
  TNtuple *ht = new TNtuple("vtxhits", "VTX hit variables", HITVARS);
  FillNTuple(vtxevents, ht);
  FillNTuple(cntevents, ht); // Eventually: use a different tree
  Printf("done.");

  // Draw changes in ladder positions
  const int NC = 3;
  TCanvas *c[NC];
  c[0] = DrawXY(tgeo, "vtx_xy", "VTX ladders", "L, dead");
  c[1] = DrawXY(tgeo, "millepede_ds", "#Delta(x,y) corrections", "L, dead, faint");
  DrawDiffs(x0,y0,z0,x1,y1,z1,"s");
  c[2] = DrawXY(tgeo, "millepede_dz", "#Deltaz corrections", "L, dead, faint");
  DrawDiffs(x0,y0,z0,x1,y1,z1,"z");

  Printf("Writing %s", pisaFileOut.Data());
  tgeo->WriteParFile(pisaFileOut.Data());

  Printf("Writing %s", rootFileOut.Data());
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);
  for (int i=0; i<NC; i++)
    c[i]->Write();

  Printf("Done!");
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

  Printf("Calling mille() in EventLoop(). Write to %s...", binfile.c_str());
  Mille m(binfile.c_str());
  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (opt.Contains("ext"))
        MilleCnt(m, trk, sgpars, zgpars, bc);
      else
        MilleVtx(m, trk, sgpars, zgpars, bc, opt);
    }

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
  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      int l = -1;

      // Ladder translation corrections
      l = Label(i,j,"x");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, mpc[l], 0., 0.);
      l = Label(i,j,"y");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, 0., mpc[l], 0.);
      l = Label(i,j,"z");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, 0. ,0., mpc[l]);

      // Ladder Phi correction from ds
      l = Label(i,j,"s");
      if (mpc.find(l) != mpc.end())
        tgeo->RotateLadderRPhi(i, j, mpc[l]);

      // Ladder Radial correction
      l = Label(i,j,"r");
      if (mpc.find(l) != mpc.end())
        tgeo->MoveLadderRadially(i, j, mpc[l]);
    }

  // Half-layer position corrections
  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<2; j++) // Arm. 0 = E; 1 = W.
    {
      int l = -1;
      // Half-layer translation corrections
      l = HalfLayerLabel(i,j,"x");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateHalfLayer(i, j, mpc[l], 0., 0.);
      l = HalfLayerLabel(i,j,"y");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateHalfLayer(i, j, 0., mpc[l], 0.);
      l = HalfLayerLabel(i,j,"z");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateHalfLayer(i, j, 0. ,0., mpc[l]);

      // Half-ladder phi correction from ds
      // s is converted to phi which is used for the half-layer rotation
      l = HalfLayerLabel(i,j,"s");
      if (mpc.find(l) != mpc.end())
        tgeo->RotateHalfLayerRPhi(i, j, mpc[l]);
    }

  for (unsigned int ev=0; ev<vtxevents.size(); ev++)
    for (unsigned int t=0; t<vtxevents[ev].size(); t++)
      vtxevents[ev][t].UpdateHits();
  for (unsigned int ev=0; ev<cntevents.size(); ev++)
    for (unsigned int t=0; t<cntevents[ev].size(); t++)
      cntevents[ev][t].UpdateHits();

  return;
}
