
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

void EventLoop(string binfile, geoEvents &events, TGraphErrors *bc = 0,
               TString opt = "");
void WriteConstFile(const char *filename, SvxTGeo *geo, TString opt = "");
void CorrectFromFile(const char *filename,
                     SvxTGeo *tgeo,
                     geoEvents &vtxevents,
                     geoEvents &cntevents);

void VtxAlign(int run = 411768,    // Run number of PRDF segment(s)
              int prod = 100,        // Production step. Starts at 0.
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

  vecs constfiles;
  vecs binfiles;
  constfiles.push_back(pedeConstFile);
  if (nStdTracks > 0)
    binfiles.push_back(pedeBinFileStd);
  if (nCntTracks > 0)
    binfiles.push_back(pedeBinFileCnt);

  TFile *inFile  = new TFile(rootFileIn.Data(), "read");
  assert(inFile);

  TNtuple *svxseg = (TNtuple *)inFile->Get("vtxhits");
  TNtuple *svxcnt = (TNtuple *)inFile->Get("cnt_clusntuple");
  assert(svxseg);
  if (nCntTracks>0)
    assert(svxcnt);

  SvxTGeo *tgeo = VTXModel(pisaFileIn.Data());

  WriteConstFile(pedeConstFile, tgeo, ""); // "empty" writes empty file
  WriteSteerFile(pedeSteerFile, binfiles, constfiles);

  // Original ladder positions
  vecd x0; vecd y0; vecd z0;
  GetLadderXYZ(tgeo, x0, y0, z0);

  geoEvents vtxevents; // Events containing SVX standalone tracks
  geoEvents cntevents; // Events containing SvxCentralTracks

  if (nStdTracks > 0)
  {
    GetEventsFromTree(svxseg, tgeo, vtxevents, -10000);

    TFile *bcf = new TFile(bcFileIn.Data(), "read");
    TGraphErrors *gbc = (TGraphErrors *) bcf->Get("gbc");

    // gbc->SetPointError(0, gbc->GetEX()[0], gbc->GetEY()[0]);
    // gbc->SetPointError(1, gbc->GetEX()[1], gbc->GetEY()[1]);
    // FitTracks(vtxevents, gbc);// Uses beam center. Doesn't work as well.

    FitTracks(vtxevents, 0);
    EventLoop(pedeBinFileStd, vtxevents, gbc);
  }
  if (nCntTracks > 0)
  {
    GetEventsFromTree(svxcnt, tgeo, cntevents);
    EventLoop(pedeBinFileCnt, cntevents);
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
EventLoop(string binfile, geoEvents &events, TGraphErrors *bc, TString opt)
{
  // Call Mille::Mille() in a loop. See MilleFunctions.h
  // Options:
  // - "ext": Assume residuals were computed from external information.

  vecs sgpars;
  vecs zgpars;

  // Select global parameters and derivatives for s residuals here
  if (true) sgpars.push_back("s");
  if (true) sgpars.push_back("x");
  if (true) sgpars.push_back("y");
  if (false) sgpars.push_back("r");

  // Select global parameters and derivatives for z residuals here
  if (true) zgpars.push_back("x");
  if (true) zgpars.push_back("y");
  if (true) zgpars.push_back("z");
  if (false) zgpars.push_back("r");

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
WriteConstFile(const char *filename, SvxTGeo *geo, TString opt)
{
  cout << "Writing " << filename << "..." << flush;
  ofstream fs(filename);

  if (opt.Contains("empty"))
  {
    fs.close();
    return;
  }

  fs << "! Write linear constraints such that sum_l (f_l * p_l) = c"
     << endl;
  fs << "! where the f_l factors provide coord. dependence to constrain shear."
     << endl;
  fs << "! The format is as follows:" << endl;
  fs << "! Constraint c" << endl;
  fs << "! 123 1.0   ! Add a p_123 * f_123 term" << endl;
  fs << "! 234 2.345 ! Add a p_234 * f_234 term" << endl;
  fs << "! ..." << endl;

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

  // Global parameter labels for z and s coordinates in east and west arms
  veci wz = LadderLabels(geo, "w", "z");
  veci ws = LadderLabels(geo, "w", "s");
  veci wr = LadderLabels(geo, "w", "r");
  veci ez = LadderLabels(geo, "e", "z");
  veci es = LadderLabels(geo, "e", "s");
  veci er = LadderLabels(geo, "e", "r");

  // // Labels for global parameters to be excluded from sum constraints
  // veci excl;

  bool constrainEdgeLadders = true;

  if (constrainEdgeLadders)
  {
    // Boundary "wedge" units - sequences of 4 edge ladders at top and bottom
    int wb[4] = {0,0,0,0};
    int wt[4] = {4,9,7,11};
    int et[4] = {5,10,8,12};
    int eb[4] = {9,19,15,23};
    veci wtz; // west top z
    veci wts; // west top s
    veci wbz; // west bottom z
    veci wbs; // west bottom s
    veci etz; // east top z
    veci ets; // east top s
    veci ebz; // east bottom z
    veci ebs; // east bottom s
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    {
      wtz.push_back(Label(lyr, wt[lyr], "z")); // west top z
      wts.push_back(Label(lyr, wt[lyr], "s")); // west top s
      wbz.push_back(Label(lyr, wb[lyr], "z")); // west bottom z
      wbs.push_back(Label(lyr, wb[lyr], "s")); // west bottom s
      etz.push_back(Label(lyr, et[lyr], "z")); // east top z
      ets.push_back(Label(lyr, et[lyr], "s")); // east top s
      ebz.push_back(Label(lyr, eb[lyr], "z")); // east bottom z
      ebs.push_back(Label(lyr, eb[lyr], "s")); // east bottom s
    }

    AddConstraint(wtz, geo, Radii, fs, "West top z r shear");
    AddConstraint(etz, geo, Radii, fs, "East top z r shear");
    AddConstraint(wts, geo, Radii, fs, "West top s r shear");
    AddConstraint(ets, geo, Radii, fs, "East top s r shear");
    AddConstraint(wbz, geo, Radii, fs, "West bottom z r shear");
    AddConstraint(ebz, geo, Radii, fs, "East bottom z r shear");
    AddConstraint(wbs, geo, Radii, fs, "West bottom s r shear");
    AddConstraint(ebs, geo, Radii, fs, "East bottom s r shear");

  } // end constrainEdgeLadders block

  // Prevent net shift or rotation of an entire arm
  AddConstraint(wz, geo, Ones, fs, "West z translation");
  AddConstraint(ez, geo, Ones, fs, "East z translation");
  AddConstraint(ws, geo, Ones, fs, "West s translation");
  AddConstraint(es, geo, Ones, fs, "East s translation");

  // Prevent various shear distortions of an entire arm
  AddConstraint(wz, geo, PhiAngles, fs, "West z phi shear");
  AddConstraint(ez, geo, PhiAngles, fs, "East z phi shear");
  AddConstraint(ws, geo, PhiAngles, fs, "West s phi shear");
  AddConstraint(es, geo, PhiAngles, fs, "East s phi shear");
  AddConstraint(wz, geo, RPhi,      fs, "West z r-phi shear");
  AddConstraint(ez, geo, RPhi,      fs, "East z r-phi shear");
  AddConstraint(ws, geo, RPhi,      fs, "West s r-phi shear");
  AddConstraint(es, geo, RPhi,      fs, "East s r-phi shear");
  AddConstraint(wz, geo, Radii,     fs, "West z r shear");
  AddConstraint(ez, geo, Radii,     fs, "East z r shear");
  AddConstraint(ws, geo, Radii,     fs, "West s r shear");
  AddConstraint(es, geo, Radii,     fs, "East s r shear");

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
