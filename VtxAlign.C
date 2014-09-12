
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
const char *pedeConstFile = "pede-const.txt";
string pedeBinFileStd = "standalone.bin";
string pedeBinFileCnt = "svxcnttrks.bin";

void EventLoop(string binfile, geoEvents &events, vecs &sgpars, vecs &zgpars,
               TGraphErrors *bc = 0, TString opt = "");
void WriteConstFile(const char *filename, vecs &sgpars, vecs &zgpars,
                    SvxTGeo *geo, TString opt = "");
void CorrectFromFile(const char *filename,
                     SvxTGeo *tgeo,
                     geoEvents &vtxevents,
                     geoEvents &cntevents);

void VtxAlign(int run = 123456,    // Run number of PRDF segment(s)
              int prod = 0,        // Production step. Starts at 0.
              int subiter = 0)     // Geometry update step. Starts at 0.
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
  vecs sgpars;
  if (0) sgpars.push_back("s");
  if (1) sgpars.push_back("x");
  if (0) sgpars.push_back("y");
  if (0) sgpars.push_back("r");

  // Select global parameters and derivatives for z residuals
  vecs zgpars;
  if (1) zgpars.push_back("x");
  if (0) zgpars.push_back("y");
  if (1) zgpars.push_back("z");
  if (0) zgpars.push_back("r");

  vecs constfiles;
  vecs binfiles;
  constfiles.push_back(pedeConstFile);
  if (nStdTracks > 0)
    binfiles.push_back(pedeBinFileStd);
  if (nCntTracks > 0)
    binfiles.push_back(pedeBinFileCnt);

  TFile *inFile = new TFile(rootFileIn.Data(), "read");
  assert(inFile);

  TNtuple *svxseg = (TNtuple *)inFile->Get("vtxhits");
  TNtuple *svxcnt = (TNtuple *)inFile->Get("cnt_clusntuple");
  assert(svxseg);
  if (nCntTracks>0)
    assert(svxcnt);

  SvxTGeo *tgeo = VTXModel(pisaFileIn.Data());

  WriteConstFile(pedeConstFile, sgpars, zgpars, tgeo, ""); // "empty" writes empty file
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
    EventLoop(pedeBinFileStd, vtxevents, sgpars, zgpars, gbc);
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
  c[1] = DrawXY(tgeo, "millepede_ds", "#Deltas corrections", "L, dead, faint");
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

  Printf("Calling mille() in EventLoop(). Write to %s...", binfile.c_str());
  Mille m(binfile.c_str());
  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (opt.Contains("ext"))
        MilleCnt(m, trk, sgpars, zgpars, bc);
      else
        MilleVtx(m, trk, sgpars, zgpars, bc);
    }

  // Mille object must go out of scope for output file to close properly.
  return;
}

void
WriteConstFile(const char *filename, vecs &sgpars, vecs &zgpars,
               SvxTGeo *geo, TString opt)
{
  cout << "Writing " << filename << "..." << flush;
  ofstream fs(filename);

  if (opt.Contains("empty"))
  {
    fs.close();
    return;
  }

  bool sdof = In(string("s"), sgpars);
  bool zdof = In(string("z"), zgpars);
  bool xdof = In(string("x"), sgpars) || In(string("x"), zgpars);
  bool ydof = In(string("y"), sgpars) || In(string("y"), zgpars);
  bool rdof = In(string("r"), sgpars) || In(string("r"), zgpars);

  fs << "! Write linear constraints such that sum_l (f_l * p_l) = c"
     << endl;
  fs << "! where the f_l factors provide coord. dependence to constrain shear."
     << endl;
  fs << "! The format is as follows:" << endl;
  fs << "! Constraint c" << endl;
  fs << "! 123 1.0   ! Add a p_123 * f_123 term" << endl;
  fs << "! 234 2.345 ! Add a p_234 * f_234 term" << endl;
  fs << "! ..." << endl;

  /*
    // Fix or regulate ladders. Example: 104 0.0 -1 (fixed)
    fs << "Parameter ! Columns: label value presigma (=0: free; >0: regularized; <0: fixed)"
       << endl;
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      for (int ldr=0; ldr<geo->GetNLadders(lyr); ldr++)
        if (Locked(lyr,ldr))
        {
          double presig = 1e-4; // Smaller values --> stronger regularization
          fs << Form("%4d 0.0 %g ! B%dL%d(s)",
                     Label(lyr, ldr, "s"), presig, lyr, ldr) << endl;
          fs << Form("%4d 0.0 %g ! B%dL%d(z)",
                     Label(lyr, ldr, "z"), presig, lyr, ldr) << endl;
        }
    fs << endl;
  */
  // Global parameter labels for z and s coordinates in east and west arms
  veci wx = LadderLabels(geo, "w", "x");
  veci wy = LadderLabels(geo, "w", "y");
  veci wz = LadderLabels(geo, "w", "z");
  veci ws = LadderLabels(geo, "w", "s");
  veci wr = LadderLabels(geo, "w", "r");
  veci ex = LadderLabels(geo, "e", "x");
  veci ey = LadderLabels(geo, "e", "y");
  veci ez = LadderLabels(geo, "e", "z");
  veci es = LadderLabels(geo, "e", "s");
  veci er = LadderLabels(geo, "e", "r");

  veci wtz = EdgeLadderLabels(geo, "w", "z", "top");
  veci wts = EdgeLadderLabels(geo, "w", "s", "top");
  veci wbz = EdgeLadderLabels(geo, "w", "z", "bottom");
  veci wbs = EdgeLadderLabels(geo, "w", "s", "bottom");
  veci etz = EdgeLadderLabels(geo, "e", "z", "top");
  veci ets = EdgeLadderLabels(geo, "e", "s", "top");
  veci ebz = EdgeLadderLabels(geo, "e", "z", "bottom");
  veci ebs = EdgeLadderLabels(geo, "e", "s", "bottom");

  // Prevent net displacements/distortions, but only if the coordinate 
  // is being used as a global parameter (otherwise the system will be 
  // rank-deficient)
  if (sdof)
  {
    AddConstraints(ws, es, geo, Ones,      fs, "s translation");
    AddConstraints(ws, es, geo, PhiAngles, fs, "s phi shear");
    AddConstraints(ws, es, geo, RPhi,      fs, "s r-phi shear");
    AddConstraints(ws, es, geo, Radii,     fs, "s r shear");
    AddConstraints(wts, ets, geo, Radii,   fs, "top s r shear");
    AddConstraints(wbs, ebs, geo, Radii,   fs, "bottom s r shear");
  }
  if (xdof)
  {
    AddConstraints(wx, ex, geo, Ones,      fs, "x translation");
    AddConstraints(wx, ex, geo, PhiAngles, fs, "x phi shear");
    AddConstraints(wx, ex, geo, RPhi,      fs, "x r-phi shear");
    AddConstraints(wx, ex, geo, Radii,     fs, "x r shear");
  }
  if (ydof)
  {
    AddConstraints(wy, ey, geo, Ones,      fs, "y translation");
    AddConstraints(wy, ey, geo, PhiAngles, fs, "y phi shear");
    AddConstraints(wy, ey, geo, RPhi,      fs, "y r-phi shear");
    AddConstraints(wy, ey, geo, Radii,     fs, "y r shear");
  }
  if (zdof)
  {
    AddConstraints(wz, ez, geo, Ones,      fs, "z translation");
    AddConstraints(wz, ez, geo, PhiAngles, fs, "z phi shear");
    AddConstraints(wz, ez, geo, RPhi,      fs, "z r-phi shear");
    AddConstraints(wz, ez, geo, Radii,     fs, "z r shear");
    AddConstraints(wtz, etz, geo, Radii,   fs, "top z r shear");
    AddConstraints(wbz, ebz, geo, Radii,   fs, "bottom z r shear");
  }
  if (rdof)
  {
    AddConstraints(wr, er, geo, Ones,      fs, "r expansion/contraction");
  }

  fs.close();
  Printf("done.");
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

  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      int l = -1;

      // Pure translation corrections
      l = Label(i,j,"x");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, mpc[l], 0., 0.);
      l = Label(i,j,"y");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, 0., mpc[l], 0.);
      l = Label(i,j,"z");
      if (mpc.find(l) != mpc.end())
        tgeo->TranslateLadder(i, j, 0. ,0., mpc[l]);

      // Phi correction from ds
      l = Label(i,j,"s");
      if (mpc.find(l) != mpc.end())
        tgeo->RotateLadderRPhi(i, j, mpc[l]);

      // Radial correction
      l = Label(i,j,"r");
      if (mpc.find(l) != mpc.end())
        tgeo->MoveLadderRadially(i, j, mpc[l]);
    }

  for (unsigned int ev=0; ev<vtxevents.size(); ev++)
    for (unsigned int t=0; t<vtxevents[ev].size(); t++)
      vtxevents[ev][t].UpdateHits();
  for (unsigned int ev=0; ev<cntevents.size(); ev++)
    for (unsigned int t=0; t<cntevents[ev].size(); t++)
      cntevents[ev][t].UpdateHits();

  return;
}
