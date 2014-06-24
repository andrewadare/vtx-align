
#include "VtxAlignmentUtils.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "SvxProj.h"
#include "Mille.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TRandom3.h>
#include <TGeoManager.h>
#include <TGeoTrack.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TPolyLine.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

typedef vector<double> vecd;
typedef vector<int>    veci;
typedef vector<string> vecs;
typedef vector<SvxGeoTrack> geoTracks;

// Globals
const double BField = 0.0;
const int nStdTracks = 999999;
const int nCntTracks = 0;
const char *pedeSteerFile = "pede-steer.txt";
const char *pedeConstFile = "pede-const.txt";
string pedeBinFileStd = "standalone.bin";
string pedeBinFileCnt = "svxcnttrks.bin";

void TrackLoop(string binfile, geoTracks &tracks, TNtuple *hitTree = 0, TNtuple *trkTree = 0, TString opt = "");
bool TrackOk(SvxGeoTrack &t);
void WriteConstFile(const char *filename, SvxTGeo *geo, TString opt = "");
void WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfile);
int Label(int layer, int ladder, string coord);
void ParInfo(int label, int &layer, int &ladder, string &coord);
void GetLadderXYZ(SvxTGeo *tgeo, vecd &x, vecd &y, vecd &z);
int GetCorrections(const char *resFile, std::map<int, double> &mpc);
void UpdateResiduals(SvxGeoTrack &track);
void GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks);
TCanvas *DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt);
void DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2);
void FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple);
bool Dead(int layer, int ladder);

bool In(int val, veci v);
void Ones(veci &labels, veci &xlabels, vecd &x);
void Radii(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x);
void PhiAngles(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x);
void RPhi(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x);
void AddConstraint(veci &labels, vecd &coords, ofstream &fs,
                   string comment = "", double sumto=0.0);

void VtxAlign(int iter = 1)
{
  // No point in continuing if Millepede II is not installed...
  if (TString(gSystem->GetFromPipe("which pede")).IsNull())
  {
    Printf("\"which pede\" returns nothing. Exiting.");
    gSystem->Exit(-1);
  }
  int run = 411768;
  TString pisaFileIn  = (iter==0) ?
                        Form("geom/svxPISA-%d.par", run) :
                        Form("geom/svxPISA-%d.par.%d", run, iter);
  TString pisaFileOut = Form("geom/svxPISA-%d.par.%d", run, iter + 1);
  TString inFileName  = (iter==0) ?
                        Form("rootfiles/%d_cluster.root", run) :
                        Form("rootfiles/%d_cluster.%d.root", run, iter);
  TString outFileName = Form("rootfiles/%d_cluster.%d.root", run, iter + 1);

  vecs constfiles;
  vecs binfiles;

  constfiles.push_back(pedeConstFile);
  if (nStdTracks > 0)
    binfiles.push_back(pedeBinFileStd);
  if (nCntTracks > 0)
    binfiles.push_back(pedeBinFileCnt);

  TFile *inFile  = new TFile(inFileName.Data(), "read");
  TFile *outFile = new TFile(outFileName.Data(), "recreate");
  const char *hitvars =
    "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid";

  TNtuple *ht1 = new TNtuple("ht1", "Pre-alignment SvxGeoHit variables",
                             hitvars);
  TNtuple *ht2 = new TNtuple("ht2", "Post-alignment SvxGeoHit variables",
                             hitvars);
  TNtuple *trktree = new TNtuple("trktree", "Tracks from VtxAlign.C",
                                 "id:y0:z0:phi:theta");
  assert(inFile);
  TNtuple *svxseg = (TNtuple *)inFile->Get("seg_clusntuple");
  TNtuple *svxcnt = (TNtuple *)inFile->Get("cnt_clusntuple");
  if (!svxseg)
  {
    svxseg = (TNtuple *)inFile->Get("ht2");
    svxseg->SetName("seg_clusntuple");
  }
  assert(svxseg);
  if (nCntTracks>0)
    assert(svxcnt);

  SvxTGeo *tgeo = new SvxTGeo;
  tgeo->ReadParFile(pisaFileIn.Data());
  tgeo->MakeTopVolume(100, 100, 100);
  tgeo->AddSensors();

  TGeoManager *mgr = tgeo->GeoManager();
  mgr->CloseGeometry();

  // TGeoVolume* top = mgr->GetTopVolume();
  // TCanvas* ctop = new TCanvas("ctop", "ctop", 1);
  // top->Draw();
  // return;

  WriteConstFile(pedeConstFile, tgeo, ""); // "empty" writes empty file
  WriteSteerFile(pedeSteerFile, binfiles, constfiles);

  // Original ladder positions
  vecd x0; vecd y0; vecd z0;
  GetLadderXYZ(tgeo, x0, y0, z0);
  geoTracks tracks;
  geoTracks cnttracks;

  GetTracksFromTree(svxseg, tgeo, tracks);
  if (nCntTracks > 0)
    GetTracksFromTree(svxcnt, tgeo, cnttracks);

  if (nStdTracks > 0)
    TrackLoop(pedeBinFileStd, tracks,    ht1, trktree, "fit"); // Fit standalone tracks

  if (nCntTracks > 0)
    TrackLoop(pedeBinFileCnt, cnttracks, ht1, trktree, ""); // Don't fit SvxCnt tracks

  // Shell out to pede executable
  gSystem->Exec(Form("pede %s", pedeSteerFile));

  // Retrieve Millepede's corrections to the global parameters
  // Key is label, value is correction.
  map<int, double> mpc;
  GetCorrections("millepede.res", mpc);

  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      // Phi correction from ds
      tgeo->RotateLadderRPhi(i, j, mpc[Label(i,j,"s")]);
      // Longitudinal (z) correction
      tgeo->TranslateLadder(i, j, 0. ,0., mpc[Label(i,j,"z")]);
    }

  for (unsigned int i=0; i<tracks.size(); i++)
    tracks[i].UpdateHits();
  for (unsigned int i=0; i<cnttracks.size(); i++)
    cnttracks[i].UpdateHits();

  // Record positions after alignment
  // Record ladder positions before alignment
  vecd x1; vecd y1; vecd z1;
  GetLadderXYZ(tgeo, x1, y1, z1);

  Printf("Post-alignment refit...");
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    double pars[4] = {0}; /* y0, z0, phi, theta */
    ZeroFieldResiduals(tracks[i], pars);
    FillNTuple(tracks[i], ht2);
  }
  for (unsigned int i=0; i<cnttracks.size(); i++)
  {
    UpdateResiduals(cnttracks[i]);
    FillNTuple(cnttracks[i], ht2);
  }

  // Draw changes in ladder positions
  const int NC = 2;
  TCanvas *c[NC];
  c[0] = DrawXY(tgeo, "vtx_xy", "VTX ladders", "L, dead");
  c[1] = DrawXY(tgeo, "millepede_dp", "Alignment corrections", "L, dead, faint");
  DrawDiffs(x0,y0,z0,x1,y1,z1);

  Printf("Writing %s", pisaFileOut.Data());
  tgeo->WriteParFile(pisaFileOut.Data());

  Printf("Writing %s", outFileName.Data());
  outFile->cd();
  ht1->Write("ht1");
  ht2->Write("ht2");
  trktree->Write();
  for (int i=0; i<NC; i++)
    c[i]->Write();

  Printf("Done.");
  return;
}

bool
TrackOk(SvxGeoTrack &t)
{
  if (t.nhits == 4)
  {
    if (t.hits[0].layer==0 &&
        t.hits[1].layer==1 &&
        t.hits[2].layer==2 &&
        t.hits[3].layer==3)
    {
      return true;
    }
  }
  return false;
}

void
TrackLoop(string binfile, geoTracks &tracks, TNtuple *hitTree, TNtuple *trkTree,
          TString opt)
{
  Printf("Track loop: fit, compute residuals, write to file...");
  Mille m(binfile.c_str());
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    // Perform straight-line fit --> residuals, track parameters
    double pars[4] = {0}; /* y0, z0, phi, theta */

    if (opt.Contains("fit"))
      ZeroFieldResiduals(tracks[i], pars);
    else
    {
      pars[0] = tracks[i].vy;
      pars[1] = tracks[i].vz;
      pars[2] = tracks[i].phi0;
      pars[3] = tracks[i].the0;
      UpdateResiduals(tracks[i]);
    }

    if (hitTree)
      FillNTuple(tracks[i], hitTree);

    if (trkTree)
    {
      // Track id and fit parameters for TTree
      float trkvars[5] = {0};
      trkvars[0] = tracks[i].hits[0].trkid;
      for (int j=0; j<4; j++)
        trkvars[j+1] = pars[j];
      trkTree->Fill(trkvars);
    }

    // Hit loop
    int slabel[1] = {0};
    int zlabel[1] = {0};
    float derlc[1] = {1.0}; // Local derivatives
    float dergl[1] = {1.0}; // Global derivatives
    float sigma_s[4] = {50e-4, 50e-4, 80e-4, 80e-4};
    float sigma_z[4] = {425e-4, 425e-4, 1000e-4, 1000e-4};

    for (int j=0; j<4; j++)
    {
      sigma_s[j] *= 1./TMath::Sqrt(12);
      sigma_z[j] *= 1./TMath::Sqrt(12);
    }

    for (int j=0; j<tracks[i].nhits; j++)
    {
      SvxGeoHit hit = tracks[i].GetHit(j);
      slabel[0] = Label(hit.layer, hit.ladder, "s");
      zlabel[0] = Label(hit.layer, hit.ladder, "z");

      if (opt.Contains("fit"))
      {
        // if (i<1000) Printf("ds %f dz %f", hit.ds, hit.dz);
        float r = hit.x*hit.x + hit.y*hit.y;
        float sderlc[4] = {1.0,   r, 0.0, 0.0};
        float zderlc[4] = {0.0, 0.0, 1.0,   r};
        float sigs = hit.xsigma * sigma_s[hit.layer];
        float sigz = hit.zsigma * sigma_z[hit.layer];
        m.mille(4, sderlc, 1, dergl, slabel, hit.ds, sigs);
        m.mille(4, zderlc, 1, dergl, zlabel, hit.dz, sigz);
      }
      else
      {
        static int cnt = 0;
        if (cnt < 10)
        {
          Printf("CNT: %d ds %f dz %f : %f %f %f ",
                 slabel[0], hit.ds, hit.dz,
                 derlc[0], dergl[0], sigma_s[j]
                );
          cnt ++;
        }
        float sigs = hit.xsigma * sigma_s[hit.layer];
        float sigz = hit.zsigma * sigma_z[hit.layer];
        m.mille(1, derlc, 1, dergl, slabel, hit.ds, sigs);
        m.mille(1, derlc, 0, dergl, slabel, 0, .000001);
        m.end(); // Write residuals for this track & reset for next one
        m.mille(1, derlc, 1, dergl, zlabel, hit.dz, sigz);
        m.mille(1, derlc, 0, dergl, slabel, 0, .000001);
        m.end(); // Write residuals for this track & reset for next one
      }
    }
    if (opt.Contains("fit"))
      m.end(); // Write residuals for this track & reset for next one
  }

  // Mille object must go out of scope for output file to close properly.
  return;
}

void
WriteConstFile(const char *filename, SvxTGeo *geo, TString opt)
{
  Printf("Writing %s", filename);
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

  // Parameter
  // 104 0.0 -1 !

  // Global parameter labels for z and s coordinates in east and west arms
  veci wz;
  veci ws;
  veci ez;
  veci es;

  // Labels for global parameters to be excluded from sum constraints
  veci excl;

  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
  {
    // West arm
    for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
    {
      wz.push_back(Label(lyr, ldr, "z"));
      ws.push_back(Label(lyr, ldr, "s"));
    }
    // East arm
    for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
    {
      ez.push_back(Label(lyr, ldr, "z"));
      es.push_back(Label(lyr, ldr, "s"));
    }
  }

  // Fill vector of excluded ladders
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr); ldr++)
      if (Dead(lyr,ldr))
      {
        excl.push_back(Label(lyr, ldr, "z"));
        excl.push_back(Label(lyr, ldr, "s"));
      }

  vecd x;
  Ones(wz,excl,x);          AddConstraint(wz, x, fs, "West z translation");
  Ones(ez,excl,x);          AddConstraint(ez, x, fs, "East z translation");
  Ones(ws,excl,x);          AddConstraint(ws, x, fs, "West s translation");
  Ones(es,excl,x);          AddConstraint(es, x, fs, "East s translation");
  Radii(wz,excl,geo,x),     AddConstraint(wz, x, fs, "West z r shear");
  Radii(ez,excl,geo,x);     AddConstraint(ez, x, fs, "East z r shear");
  Radii(ws,excl,geo,x);     AddConstraint(ws, x, fs, "West s r shear");
  Radii(es,excl,geo,x);     AddConstraint(es, x, fs, "East s r shear");
  PhiAngles(wz,excl,geo,x); AddConstraint(wz, x, fs, "West z phi shear");
  PhiAngles(ez,excl,geo,x); AddConstraint(ez, x, fs, "East z phi shear");
  PhiAngles(ws,excl,geo,x); AddConstraint(ws, x, fs, "West s phi shear");
  PhiAngles(es,excl,geo,x); AddConstraint(es, x, fs, "East s phi shear");
  RPhi(wz,excl,geo,x);      AddConstraint(wz, x, fs, "West z r-phi shear");
  RPhi(ez,excl,geo,x);      AddConstraint(ez, x, fs, "East z r-phi shear");
  RPhi(ws,excl, geo,x);     AddConstraint(ws, x, fs, "West s r-phi shear");
  RPhi(es,excl, geo,x);     AddConstraint(es, x, fs, "East s r-phi shear");

  fs.close();
  Printf("Done.");
  return;
}

void
WriteSteerFile(const char *filename, vecs &binfiles, vecs &constfiles)
{
  Printf("Writing %s", filename);

  ofstream fs(filename);

  fs << Form("* This is %s created by %s", filename, "SimVtxAlign.C") << endl;
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

  fs << "method inversion 10 0.0001  ! Gauss. elim., #iterations, tol." << endl;
  fs << "end" << endl;

  fs.close();

  Printf("Done.");

  return;
}

int
Label(int layer, int ladder, string coord)
{
  // Return a global parameter label for this layer, ladder, and coordinate.
  // In Millepede II, any unique integer > 0 will do.
  // Labels are not required to be sequential.
  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3
  int ic = 99;

  if (coord.compare("s") == 0)
    ic = 0;
  if (coord.compare("z") == 0)
    ic = 1;
  // Add other coordinates / degrees of freedom here as needed
  // Then ParInfo() would need corresponding modification.

  return 70*ic + start[layer] + ladder + 1;
}

void
ParInfo(int label, int &layer, int &ladder, string &coord)
{
  // Get layer, ladder, and coordinate string from global parameter label.
  // Inverse of Label() function.

  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3
  int l = 0;

  if (label > 70)
  {
    coord = "z";
    l = label - 70;
  }
  else
  {
    coord = "s";
    l = label;
  }

  layer = 3;
  for (int i=3; i>=0; i--)
    if (l < start[i] + 1)
      layer = i - 1;

  ladder = l - start[layer] - 1;

  return;
}

TCanvas *
DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt)
{
  TCanvas *c = new TCanvas(name, title, 900, 900);
  TH1F *hf = c->DrawFrame(-20, -20, 20, 20, title);
  hf->SetXTitle("East                 x [cm]                 West");
  hf->GetXaxis()->CenterTitle();
  hf->SetYTitle("y [cm]");
  hf->GetYaxis()->CenterTitle();

  TLatex label;
  label.SetTextSize(0.02);
  double xyz[3] = {0};
  for (int i=0; i<geo->GetNLayers(); i++)
    for (int j=0; j<geo->GetNLadders(i); j++)
    {

      if (opt.Contains("dead") && Dead(i,j)) // Mark dead ladders
      {
        double xyz[3] = {0};
        geo->GetSensorXYZ(i,j,0,xyz);

        TLatex ltx;
        ltx.SetTextSize(0.04);
        ltx.SetTextAlign(22);
        ltx.SetTextColor(kRed-4);
        ltx.DrawLatex(xyz[0], xyz[1], "#times");
      }

      TPolyLine *s = geo->LadderOutlineXY(i,j);
      s->SetLineColor(opt.Contains("faint") ? kGray : kGray+2);
      s->SetLineWidth(opt.Contains("faint") ? 1 : 2);
      s->Draw("same");

      if (opt.Contains("faint"))
        continue;

      geo->GetSensorXYZ(i, j, 0, xyz);
      int horz = xyz[0] > 0 ? 1 : 3;
      int vert = xyz[1] > 0 ? 1 : 3;
      label.SetTextAlign(10*horz + vert);
      label.DrawLatex(xyz[0], xyz[1], Form(" %d ", j));


    }
  return c;
}

void
GetLadderXYZ(SvxTGeo *tgeo, vecd &x, vecd &y, vecd &z)
{
  int ipt = 0;
  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      double xyz[3] = {0};
      tgeo->GetSensorXYZ(i,j,0,xyz);
      x.push_back(xyz[0]);
      y.push_back(xyz[1]);
      z.push_back(xyz[2]);

      if (false)
        Printf("xyz %f %f %f", x[ipt], y[ipt], z[ipt]);

      ipt++;
    }

  return;
}

void
DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2)
{
  int n = (int)x1.size();

  for (int i=0; i<n; i++)
  {
    double f = 20;
    double x = x1[i], y = y1[i];
    double dx = x2[i] - x;
    double dy = y2[i] - y;
    double dz = z2[i] - z1[i];

    // Draw points showing displacement in z coordinate
    TMarker mkr;
    mkr.SetMarkerStyle(kOpenCircle);
    mkr.SetMarkerColor(kGray+1);
    mkr.SetMarkerSize(0.5);
    //    mkr.DrawMarker(x,y);

    if (TMath::Abs(dz) > 5e-4) // Only label changes > 5 um
    {
      mkr.SetMarkerStyle(kOpenCircle);
      mkr.SetMarkerColor(dz>0 ? kRed-4 : kAzure-4); // Red/blue shift mnemonic
      mkr.SetMarkerSize(100*TMath::Abs(dz)); // 1 = 8px diam, 2 = 16px, ...
      mkr.DrawMarker(x,y);

      TLatex ltx;
      ltx.SetTextSize(0.018);
      ltx.SetTextAlign(22);
      ltx.SetTextColor(dz>0 ? kRed+2 : kAzure+3);
      ltx.DrawLatex(x, y, Form("%.0f", 1e4*dz));
    }

    // Draw arrows showing (significant) displacements in xy plane
    if (dx*dx + dy*dy > 5e-4) // Only label changes > 5 um
    {
      TArrow a;
      a.SetLineWidth(2);
      a.DrawArrow(x, y, x + f*dx, y + f*dy, 0.005);

      TLatex ltx;
      ltx.SetTextSize(0.01);
      ltx.DrawLatex(x + f*dx, y + f*dy, Form("(%.0f, %.0f)", 1e4*dx, 1e4*dy));
    }
  }

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
    Printf("Error opening %s", resFile);
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

void
UpdateResiduals(SvxGeoTrack &track)
{
  static TRandom3 ran;

  for (int i=0; i<track.nhits; i++)
  {
    SvxGeoHit hit = track.GetHit(i);
    double r = TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);
    double xproj = track.vx + r*TMath::Cos(track.phi0);
    double yproj = track.vy + r*TMath::Sin(track.phi0);
    double zproj = track.vz + r/TMath::Tan(track.the0);

    double phiproj = TMath::ATan2(yproj, xproj);
    double phihit = TMath::ATan2(hit.y, hit.x);

    track.hits[i].dz = zproj - hit.z ;
    track.hits[i].ds = r*fmod(phiproj - phihit, TMath::TwoPi()) ;
  }

  return;
}

void
GetTracksFromTree(TNtuple *t, SvxTGeo *geo, geoTracks &tracks)
{
  // This function reads hit ntuple variables of the form
  // "layer:ladder:sensor:xs:ys:zs:x:y:z:xsigma:zsigma:dz:ds:trkid"
  // into the tracks vector.

  assert(t);
  assert(geo);

  long nentries = t->GetEntries();
  int previd = -1;

  Printf("Reading tracks from NTuple (%d entries)...", (int)nentries);
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int id = t->GetLeaf("trkid")->GetValue();
    if (id != previd) // New track
    {
      SvxGeoTrack t;
      // Track info is not stored in tree. If it were, it would go here.
      tracks.push_back(t);
    }

    SvxGeoHit hit;
    hit.layer  = t->GetLeaf("layer")->GetValue();
    hit.ladder = t->GetLeaf("ladder")->GetValue();
    hit.sensor = t->GetLeaf("sensor")->GetValue();
    hit.xs     = t->GetLeaf("lx")->GetValue();
    hit.ys     = t->GetLeaf("ly")->GetValue();
    hit.zs     = t->GetLeaf("lz")->GetValue();
    hit.x      = t->GetLeaf("gx")->GetValue();
    hit.y      = t->GetLeaf("gy")->GetValue();
    hit.z      = t->GetLeaf("gz")->GetValue();
    hit.xsigma = t->GetLeaf("x_size")->GetValue();
    hit.zsigma = t->GetLeaf("z_size")->GetValue();
    hit.dz     = t->GetLeaf("res_z")->GetValue();
    hit.ds     = t->GetLeaf("res_s")->GetValue();
    hit.trkid  = t->GetLeaf("trkid")->GetValue();
    hit.node   = geo->SensorNode(hit.layer, hit.ladder, hit.sensor);
    tracks.back().nhits++;
    tracks.back().hits.push_back(hit);
    previd = id;
  }

  Printf("%d tracks imported.", (int)tracks.size());
  return;
}

void
FillNTuple(SvxGeoTrack &gt, TNtuple *ntuple)
{
  // TNtuple *ntuple = new TNtuple("t", "SvxGeoHit variables",
  // "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:x_size:z_size:res_z:res_s:trkid"


  assert(ntuple);

  if (false)
    Printf("ds (%8.4f,%8.4f,%8.4f,%8.4f),  "
           "dz (%8.4f,%8.4f,%8.4f,%8.4f)",
           gt.hits[0].ds, gt.hits[1].ds, gt.hits[2].ds, gt.hits[3].ds,
           gt.hits[0].dz, gt.hits[1].dz, gt.hits[2].dz, gt.hits[3].dz);

  for (int ihit=0; ihit<gt.nhits; ihit++)
  {
    int nj = 14;
    std::vector<float> vars(nj, 0.);
    int j = 0;
    SvxGeoHit hit = gt.GetHit(ihit);

    vars[j++] = hit.layer  ;
    vars[j++] = hit.ladder ;
    vars[j++] = hit.sensor ;
    vars[j++] = hit.xs     ;
    vars[j++] = hit.ys     ;
    vars[j++] = hit.zs     ;
    vars[j++] = hit.x      ;
    vars[j++] = hit.y      ;
    vars[j++] = hit.z      ;
    vars[j++] = hit.xsigma ;
    vars[j++] = hit.zsigma ;
    vars[j++] = hit.dz     ;
    vars[j++] = hit.ds     ;
    vars[j++] = hit.trkid  ;
    ntuple->Fill(&vars[0]);
  }

  return;
}

bool
Dead(int layer, int ladder)
{
  if (layer==1 && ladder==11) return true;
  if (layer==3 && ladder==10) return true;
  if (layer==3 && ladder==16) return true;
  if (layer==3 && ladder==23) return true;
  return false;
}

bool In(int val, veci v)
{
  for (unsigned int j=0; j<v.size(); j++)
    if (val==v[j])
      return true;
  return false;
}

void
Ones(veci &labels, veci &xlabels, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  for (unsigned int i=0; i<labels.size(); i++)
    x[i] = In(labels[i], xlabels) ? 0.0 : 1.0;
}

void
Radii(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    x[i] = In(labels[i], xlabels) ? 0.0 : geo->SensorRadius(lyr, ldr, 0);
  }
}

void
PhiAngles(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    float phi = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0),
                                        TMath::TwoPi());
    x[i] = In(labels[i], xlabels) ? 0.0 : phi;
  }
}

void
RPhi(veci &labels, veci &xlabels, SvxTGeo *geo, vecd &x)
{
  x.clear();
  x.resize(labels.size(), 0.0);
  int lyr, ldr;
  string s;
  for (unsigned int i=0; i<labels.size(); i++)
  {
    ParInfo(labels[i], lyr, ldr, s);
    float phi = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0),
                                        TMath::TwoPi());
    float rphi = phi*geo->SensorRadius(lyr, ldr, 0);
    x[i] = In(labels[i], xlabels) ? 0.0 : rphi;
  }
}

void
AddConstraint(veci &labels, vecd &coords, ofstream &fs, string comment, double sumto)
{
  fs << "Constraint " << sumto << "  ! " << comment << endl;
  for (unsigned int i=0; i<labels.size(); i++)
    fs << labels[i] << " " << coords[i] << endl;

  return;
}
