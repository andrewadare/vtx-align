
#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "ParameterDefs.h"
#include "ConstraintBuilder.h"
#include "VtxIO.h"
#include "VtxVis.h"

using namespace std;

// Globals
bool allowRShifts = false;   // Switch to add radial d.o.f. to global pars
const double BField = 0.0;
const int nStdTracks = 999999;
const int nCntTracks = 0;
const char *pedeSteerFile = "pede-steer.txt";
const char *pedeConstFile = "pede-const.txt";
string pedeBinFileStd = "standalone.bin";
string pedeBinFileCnt = "svxcnttrks.bin";

void TrackLoop(string binfile, geoTracks &tracks, TNtuple *hitTree = 0,
               TNtuple *trkTree = 0, TString opt = "");
void EventLoop(string binfile, geoEvents &events, TString opt = "");
void WriteConstFile(const char *filename, SvxTGeo *geo, TString opt = "");
void WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfile);
int GetCorrections(const char *resFile, std::map<int, double> &mpc);
void FilterTracks(geoTracks &a, geoTracks &b, double maxdca);
void MilleVtxStandalone(Mille &m, SvxGeoTrack &trk);
void MilleVtxExternal(Mille &m, SvxGeoTrack &trk);

void VtxAlign(int run = 411768,    // Run number of PRDF segment(s)
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
    GetEventsFromTree(svxseg, tgeo, vtxevents, -1);

    TFile *bcf = new TFile(bcFileIn.Data(), "read");
    TGraphErrors *gbc = (TGraphErrors *) bcf->Get("gbc");

    gbc->SetPointError(0, gbc->GetEX()[0], gbc->GetEY()[0]);
    gbc->SetPointError(1, gbc->GetEX()[1], gbc->GetEY()[1]);

    // FitTracks(vtxevents, gbc);// Uses beam center. Doesn't work as well.

    FitTracks(vtxevents, 0);
    EventLoop(pedeBinFileStd, vtxevents);
  }
  if (nCntTracks > 0)
  {
    GetEventsFromTree(svxcnt, tgeo, cntevents);
    EventLoop(pedeBinFileCnt, cntevents);
  }

  // Shell out to pede executable
  gSystem->Exec(Form("pede %s", pedeSteerFile));

  // Retrieve Millepede's corrections to the global parameters
  // Key is label, value is correction.
  map<int, double> mpc;
  GetCorrections("millepede.res", mpc);

  // Copy geometry adjustments (with their errors) to logs/ directory.
  gSystem->Exec(Form("cp millepede.res logs/dp-%d-%d-%d.txt",
                     run, prod, subiter));

  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      // Phi correction from ds
      tgeo->RotateLadderRPhi(i, j, mpc[Label(i,j,"s")]);
      // Longitudinal (z) correction
      tgeo->TranslateLadder(i, j, 0. ,0., mpc[Label(i,j,"z")]);
      // Radial correction
      if (allowRShifts)
        tgeo->MoveLadderRadially(i, j, mpc[Label(i,j,"r")]);
    }

  for (unsigned int ev=0; ev<vtxevents.size(); ev++)
    for (unsigned int t=0; t<vtxevents[ev].size(); t++)
      vtxevents[ev][t].UpdateHits();
  for (unsigned int ev=0; ev<cntevents.size(); ev++)
    for (unsigned int t=0; t<cntevents[ev].size(); t++)
      cntevents[ev][t].UpdateHits();

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
MilleVtxStandalone(Mille &m, SvxGeoTrack &trk)
{
  // Mille::mille() is called once for each residual (so, twice per hit:
  // one for the s coordinate, one for z).
  // These calls fill a binary file (standalone.bin) with info on
  // measured and expected residual sizes, and how the residuals change with
  // respect to local (track) and global (detector position) parameters.
  // The actual fitting is performed afterward by the pede program.

  for (int j=0; j<trk.nhits; j++)
  {
    SvxGeoHit hit = trk.GetHit(j);
    int ls = Label(hit.layer, hit.ladder, "s");
    int lz = Label(hit.layer, hit.ladder, "z");
    int lr = Label(hit.layer, hit.ladder, "r");
    int slabels[2] = {ls, lr};
    int zlabels[2] = {lz, lr};
    float r = hit.x*hit.x + hit.y*hit.y;
    float theta = TMath::ATan2(r, hit.z);

    // Assign sigmas for denominator of chi square function.
    // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
    // If millepede complains that chi^2/ndf is away from 1.0,
    // this is a good place to make adjustments.
    float sigs = 4. * hit.xsigma * ClusterXResolution(hit.layer);
    float sigz = 4. * hit.zsigma * ClusterZResolution(hit.layer);

    if (false)
      Printf("hit.ds %.3g, sigs %.3g, hit.dz %.3g, sigz %.3g",
             hit.ds, sigs, hit.dz, sigz);

    // Local derivatives
    // Split into two independent fits per track:
    // 1. in the r-z plane: z = z0 + slope*r
    // 2. in the x-y plane: y' = y0' + slope*r (where y' \perp r)
    float sderlc[4] = {1.0,   r, 0.0, 0.0}; // dy(r)/dy0, dy(r)/dslope
    float zderlc[4] = {0.0, 0.0, 1.0,   r}; // dz(r)/dz0, dz(r)/dslope

    // Approximate global derivatives
    float dsdr = hit.ds/r;
    float dzdr = TMath::Tan(theta - TMath::PiOver2());

    // Use ngd to control whether derivatives of residuals
    // w.r.t. radial coordinate are included in global fit.
    int ngd = allowRShifts ? 2 : 1;

    float sdergl[2] = {1., dsdr}; // ds/ds, ds/dr
    float zdergl[2] = {1., dzdr}; // dz/dz, dz/dr
    m.mille(4, sderlc, ngd, sdergl, slabels, hit.ds, sigs);
    m.mille(4, zderlc, ngd, zdergl, zlabels, hit.dz, sigz);
  }

  // Write residuals for this track to file & reset for next track.
  m.end();

  return;
}

void
MilleVtxExternal(Mille &m, SvxGeoTrack &trk)
{
  for (int j=0; j<trk.nhits; j++)
  {
    SvxGeoHit hit = trk.GetHit(j);
    int ls = Label(hit.layer, hit.ladder, "s");
    int lz = Label(hit.layer, hit.ladder, "z");
    int lr = Label(hit.layer, hit.ladder, "r");
    int slabels[2] = {ls, lr};
    int zlabels[2] = {lz, lr};

    float sderlc[1] = {1.0}; // Local derivatives
    float zderlc[1] = {1.0}; // Local derivatives
    float sdergl[1] = {1.0}; // Global derivatives
    float zdergl[1] = {1.0}; // Global derivatives

    // Assign sigmas for denominator of chi square function.
    // Note: expecting that hit.{x,z}sigma = {x,z}_size: 1,2,3....
    // If millepede complains that chi^2/ndf is away from 1.0,
    // this is a good place to make adjustments.
    float sigs = 4. * hit.xsigma * ClusterXResolution(hit.layer);
    float sigz = 4. * hit.zsigma * ClusterZResolution(hit.layer);

    // Here, fit clusters individually. Each cluster is treated as 
    // one "local fit object".
    // Give s residual perfect resolution
    m.mille(1, sderlc, 0, sdergl, slabels, 0, .000001);
    m.mille(1, zderlc, 1, zdergl, zlabels, hit.dz, sigz);
    m.end();

    // Give z residual perfect resolution
    m.mille(1, sderlc, 1, sdergl, slabels, hit.ds, sigs);
    m.mille(1, zderlc, 0, zdergl, zlabels, 0, .000001);
    m.end();
  }
  return;
}

void
EventLoop(string binfile, geoEvents &events, TString opt)
{
  // Options:
  // - "ext": Assume residuals were computed from external information.

  Printf("Calling mille() in EventLoop(). Write to %s...", binfile.c_str());
  Mille m(binfile.c_str());
  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      MilleVtxStandalone(m, trk);
      if (opt.Contains("ext"))
        MilleVtxExternal(m, trk);
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

  // Write linear constraints such that \sum_label w_label * p_label = c
  // where w is the weight to be applied for each term.
  //
  // The format is as follows:
  // Constraint c
  // label w_label
  // label w_label
  // ...

  // Fix ladders. Example: 104 0.0 -1
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

  // Apply presigma to regulate radial movement.
  // Only apply radial constraints if radial corrections are allowed.
  if (allowRShifts)
  {
    fs << "Parameter ! Columns: label value presigma (=0: free; >0: regularized; <0: fixed)"
       << endl;
    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
      for (int ldr=0; ldr<geo->GetNLadders(lyr); ldr++)
      {
        double presig = -1.; // Radial DOF fixed by default
        if (allowRShifts)
          presig = 1e-4; // Smaller (positive) values <-> stronger regularization
        fs << Form("%4d 0.0 %g ! B%dL%d(r)",
                   Label(lyr, ldr, "r"), presig, lyr, ldr) << endl;
      }
    fs << endl;
  }

  vecd x;

  // Global parameter labels for z and s coordinates in east and west arms
  veci wz;
  veci ws;
  veci wr;
  veci ez;
  veci es;
  veci er;
  veci vtxr;

  // Labels for global parameters to be excluded from sum constraints
  veci excl;

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

    for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    {
      // West arm
      for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
      {
        ws.push_back(Label(lyr, ldr, "s"));
        wz.push_back(Label(lyr, ldr, "z"));
        wr.push_back(Label(lyr, ldr, "r"));
        vtxr.push_back(Label(lyr, ldr, "r"));
      }
      // East arm
      for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
      {
        es.push_back(Label(lyr, ldr, "s"));
        ez.push_back(Label(lyr, ldr, "z"));
        er.push_back(Label(lyr, ldr, "r"));
        vtxr.push_back(Label(lyr, ldr, "r"));
      }
    }

    Radii(wtz,excl,geo,x);  AddConstraint(wtz, x, fs, "West top z r shear");
    Radii(etz,excl,geo,x);  AddConstraint(etz, x, fs, "East top z r shear");
    Radii(wts,excl,geo,x);  AddConstraint(wts, x, fs, "West top s r shear");
    Radii(ets,excl,geo,x);  AddConstraint(ets, x, fs, "East top s r shear");
    Radii(wbz,excl,geo,x);  AddConstraint(wbz, x, fs, "West bottom z r shear");
    Radii(ebz,excl,geo,x);  AddConstraint(ebz, x, fs, "East bottom z r shear");
    Radii(wbs,excl,geo,x);  AddConstraint(wbs, x, fs, "West bottom s r shear");
    Radii(ebs,excl,geo,x);  AddConstraint(ebs, x, fs, "East bottom s r shear");
  } // end constrainEdgeLadders block

  // Fill vector of excluded ladders
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr); ldr++)
      if (Dead(lyr,ldr) || Locked(lyr,ldr))
      {
        excl.push_back(Label(lyr, ldr, "z"));
        excl.push_back(Label(lyr, ldr, "s"));
      }

  if (true)
  {
    // Prevent net radial displacement over all vtx ladders
    // Ones(vtxr,excl,x);  AddConstraint(vtxr, x, fs, "Full-VTX radial expansion/contraction");

    // Prevent net shift or rotation of an entire arm
    Ones(wz,excl,x);  AddConstraint(wz, x, fs, "West z translation");
    Ones(ez,excl,x);  AddConstraint(ez, x, fs, "East z translation");
    Ones(ws,excl,x);  AddConstraint(ws, x, fs, "West s translation");
    Ones(es,excl,x);  AddConstraint(es, x, fs, "East s translation");

    // Only apply radial constraints if radial corrections are allowed
    if (allowRShifts)
    {
      Ones(wr,excl,x);  AddConstraint(wr, x, fs, "West r expansion/contraction");
      Ones(er,excl,x);  AddConstraint(er, x, fs, "East r expansion/contraction");
    }

    // Prevent various shear distortions of an entire arm
    PhiAngles(wz,excl,geo,x); AddConstraint(wz, x, fs, "West z phi shear");
    PhiAngles(ez,excl,geo,x); AddConstraint(ez, x, fs, "East z phi shear");
    PhiAngles(ws,excl,geo,x); AddConstraint(ws, x, fs, "West s phi shear");
    PhiAngles(es,excl,geo,x); AddConstraint(es, x, fs, "East s phi shear");
    RPhi(wz,excl,geo,x);      AddConstraint(wz, x, fs, "West z r-phi shear");
    RPhi(ez,excl,geo,x);      AddConstraint(ez, x, fs, "East z r-phi shear");
    RPhi(ws,excl,geo,x);      AddConstraint(ws, x, fs, "West s r-phi shear");
    RPhi(es,excl,geo,x);      AddConstraint(es, x, fs, "East s r-phi shear");
    Radii(wz,excl,geo,x);     AddConstraint(wz, x, fs, "West z r shear");
    Radii(ez,excl,geo,x);     AddConstraint(ez, x, fs, "East z r shear");
    Radii(ws,excl,geo,x);     AddConstraint(ws, x, fs, "West s r shear");
    Radii(es,excl,geo,x);     AddConstraint(es, x, fs, "East s r shear");
  }

  fs.close();
  Printf("done.");
  return;
}

void
WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfiles)
{
  cout << "Writing " << filename << "..." << flush;

  ofstream fs(filename);

  fs << Form("* This is %s created by %s", filename, "VtxAlign.C") << endl;
  fs << Form("* Pass this file to pede: pede %s", filename) << endl;
  fs << endl;

  fs << "Fortranfiles  ! Fortran/text inputs listed here:" << endl;
  for (unsigned int i=0; i<constfiles.size(); i++)
    fs << constfiles[i] << " ! constraints file" << endl;
  fs << endl;

  fs << "Cfiles  ! c/c++ binary input files listed here:" << endl;
  for (unsigned int i=0; i<binfiles.size(); i++)
    fs << binfiles[i] << " ! binary data file" << endl;
  fs << endl;

  fs << "method inversion 5 0.0001  ! Gauss. elim., #iterations, tol." << endl;
  fs << "end" << endl;

  fs.close();

  Printf("done.");

  return;
}

int
GetCorrections(const char *resFile, std::map<int, double> &mpc)
{
  std::ifstream filein(resFile);
  int label;
  double p, col3, col4, col5;

  if (!filein)
  {
    Error("GetCorrections() in VtxAlign.C", "Problem opening %s", resFile);
    return -1;
  }
  else
    for (std::string line; std::getline(filein, line);)
    {
      if (filein >> label >> p >> col3 >> col4 >> col5)
        mpc[label] = p;
    }

  return (int)mpc.size();
}
