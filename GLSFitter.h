#ifndef __GLSFITTER_H__
#define __GLSFITTER_H__

#include "VtxAlignBase.h"
#include <TDecompSVD.h>

using namespace std;
typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L);
TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L, TMatrixD &cov);
void TrackFitL(SvxGeoTrack &gt, TString opt = ""); // Longitudinal component --> z0, theta
void TrackFitT(SvxGeoTrack &gt, TGraphErrors *bc = 0); // Transverse component --> y0', phi
void FitTracks(geoTracks &tracks, TGraphErrors *bc = 0);
void FitTracks(geoEvents &events, TGraphErrors *bc = 0, TString opt = "");
void CalculateDCA(geoTracks &event, TGraphErrors *bc, TVectorD &xpz, TString opt = "");

// It is questionable whether the following functions belong in this file.
bool East(double phi);
double ClusterXResolution(int layer);
double ClusterZResolution(int layer);
double ZVertexResolution();
TVectorD XYCenter(geoTracks &event, TString arm, int ntrk = -1, TString opt="");
TVectorD ZVertexGLS(geoTracks &event, TString arm, int ntrk = -1, TString opt="");
TVectorD IPVec(TVectorD &a, TVectorD &n, TVectorD &p);
TVectorD IPVec(SvxGeoTrack &t, TVectorD &p);
void FindVertex(geoTracks &event, TString opt = "xyz");
void FindVertexEastWest(geoTracks &event, TString opt = "xyz");
TVectorD RetrieveVertex(geoTracks &event, TString opt);

TVectorD
SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L, TMatrixD &cov)
{
  // Simple generalized linear least-squares fitter.
  // Solve y(X) = beta'*X + error(L) for beta (vector of regression coefs.)
  // L is the inverse covariance matrix for y.
  // Least-squares solution involves solving X' L X beta = X' L y

  double tol = 1.e-15;
  TMatrixD XT(X); XT.T();
  TMatrixD A = XT * L * X;
  TVectorD b = XT * L * y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();
  TMatrixD V  = svd.GetV();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  int n = s.GetNrows();
  TMatrixD Sd(n, n);
  for (int i = 0; i < n; i++)
    Sd(i, i) = s(i) > tol ? 1. / s(i) : 0.;

  // Result
  TVectorD beta = V * Sd * UT * b;

  // and covariance matrix
  V *= Sd;
  TMatrixD C(V, TMatrixD::kMultTranspose, V);
  cov.ResizeTo(C);
  cov = C;

  return beta;
}

TVectorD
SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L)
{
  // Simple generalized linear least-squares fitter.
  // Solve y(X) = beta'*X + error(L) for beta (vector of regression coefs.)
  // L is the inverse covariance matrix for y.
  // Least-squares solution involves solving X' L X beta = X' L y

  TMatrixD XT(X); XT.T();
  TMatrixD A = XT * L * X;
  TVectorD b = XT * L * y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  TMatrixD Sd(s.GetNrows(), s.GetNrows());
  for (int i = 0; i < s.GetNrows(); i++)
    Sd(i, i) = s(i) > 0 ? 1. / s(i) : 0.;

  TVectorD beta = svd.GetV() * Sd * UT * b;

  return beta;
}

void
TrackFitL(SvxGeoTrack &gt, TString opt)
{
  // Longitudinal/polar component
  // Fit track using a straight line z(x') = z0 + c*x'.
  // Assigns residuals, z-intercept, and polar angle.

  // m = # measurements; n = # parameters.
  int m = gt.nhits, n = 2;
  if (opt.Contains("fit_to_z_vertex"))
    m += 1;

  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);
  TVectorD xp(m); // This is x', the abcissa in the rotated system.

  for (int ihit = 0; ihit < gt.nhits; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    xp(ihit) = cos(gt.phirot)*hit.x + sin(gt.phirot)*hit.y;
    X(ihit, 0) = 1;
    X(ihit, 1) = xp(ihit);

    // Expecting that hit.zsigma is z_size in integer units 1,2,3....
    double zsigma = hit.zsigma * ClusterZResolution(hit.layer);
    Cinv(ihit, ihit) = (zsigma > 0) ? 1. / zsigma : 1.;
    y(ihit) = hit.z;
  }

  if (opt.Contains("fit_to_z_vertex"))
  {
    // Append vertex x' and z to the system
    X(m-1, 0) = 1;
    X(m-1, 1) = cos(gt.phirot)*gt.vx + sin(gt.phirot)*gt.vy;
    y(m-1) = gt.vz;
    Cinv(m-1, m-1) = 1./ZVertexResolution();
  }

  TVectorD beta = SolveGLS(X, y, Cinv);
  double z0 = beta(0), c = beta(1);

  for (int ihit = 0; ihit < gt.nhits; ihit++)
  {
    // Assign z residual: projected - measured
    SvxGeoHit hit = gt.GetHit(ihit);
    gt.hits[ihit].dz = z0 + c*xp(ihit) - hit.z;
  }

  gt.z0 = z0;
  gt.the0 = TMath::ATan2(1.0, c);

  return;
}

void
TrackFitT(SvxGeoTrack &gt, TGraphErrors *bc)
{
  // Transverse/azimuthal component
  // Perform straight-line fit y' = m'*x' + b' after rotating points
  // to approximately lie along the x axis. Then by construction,
  // x' ~ r, m' ~ 0, and the error is (mostly) in the y' = phi direction.
  // Assign s = r*phi residual and azimuthal angle and intercept info.

  int m = gt.nhits;
  int n = 2;

  if (bc)
    m += 1; // Add beamcenter as a new datapoint

  TMatrixD points(n, m); // Datapoints. Columns are x',y' pairs
  TMatrixD X(m, n);      // Column 0 is 1's, column 1 is x' (~r) hit coords.
  TMatrixD Cinv(m, m);   // Inverse covariance matrix. Currently diagonal.
  TVectorD y(m);         // Dependent variables y'.
  double phirot = 0;     // <phi_cluster> used to rotate clusters

  int outermostLayer = 0;
  for (int ihit = 0; ihit < gt.nhits; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    points(0, ihit) = hit.x;
    points(1, ihit) = hit.y;

    X(ihit, 0) = 1;
    // Note that hit.xsigma is expected to be in integer units 1,2,3....
    double xsigma = hit.xsigma * ClusterXResolution(hit.layer);
    Cinv(ihit, ihit) = (xsigma > 0) ? 1. / xsigma : 1.;

    // Compute rotation angle as phi of outermost hit
    if (hit.layer > outermostLayer)
    {
      outermostLayer = hit.layer;
      phirot = TMath::ATan2(hit.y, hit.x);
    }
  }
  if (phirot < 0)
    phirot += TMath::TwoPi();

  if (bc)
  {
    // Append beamcenter to the system
    int arm = (gt.hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.
    points(0, m-1) = bc->GetX()[arm];
    points(1, m-1) = bc->GetY()[arm];
    X(m-1, 0) = 1;
    float sigbc = bc->GetEX()[arm];
    assert(sigbc > 0);
    Cinv(m-1, m-1) = 1./sigbc;
  }

  // Rotate x, y by -phirot so error is approximately in y' direction only.
  // In rotated frame (xp, yp) = (cx + sy, cy - sx)
  TMatrixD R(2, 2);
  double c = cos(phirot), s = sin(phirot);
  R(0, 0) =  c; R(0, 1) = s;
  R(1, 0) = -s; R(1, 1) = c;

  points = R * points;

  TMatrixDColumn(X, 1) = TMatrixDRow(points, 0);
  y = TMatrixDRow(points, 1);
  // Fit track to get beta prime = [b', m'] in rotated system.
  TVectorD betap = SolveGLS(X, y, Cinv);
  double bp = betap(0), mp = betap(1); // Intercept and slope, both small.

  // // Rotate back to get y-intercept b and slope m of track.
  // double y0 = bp / (c - mp*s);
  // double slope = (s + mp*c) / (c - mp*s);

  double phi = phirot + TMath::ATan(mp);
  if (phi < 0)
    phi += TMath::TwoPi();

  // Assign residuals to hits
  for (int ihit = 0; ihit < gt.nhits; ihit++)
  {
    double xp  = points(0, ihit); // x' in rotated system--always positive
    assert(xp > 0);

    double ds = mp * xp + bp - points(1, ihit); // projected - measured y'
    gt.hits[ihit].ds = ds;
  }

  // Assign track variables
  gt.phi0 = phi;
  gt.phirot = phirot;
  gt.yp0 = bp;

  return;
}

void
FitTracks(geoTracks &tracks, TGraphErrors *bc)
{
  cout << Form("Fitting %lu tracks...", tracks.size()) << flush;

  for (unsigned int i = 0; i < tracks.size(); i++)
  {
    TrackFitT(tracks[i], bc); // Call first
    TrackFitL(tracks[i]); // Call second
  }

  Printf("done.");
  return;
}

void
FindVertex(geoTracks &event, TString opt)
{
  // Compute vertex position from both VTX arms.
  // Requires that tracks have already been fit.

  TVectorD xy = XYCenter(event, "");
  TVectorD rz = ZVertexGLS(event, "");

  for (unsigned int t = 0; t < event.size(); t++)
  {
    SvxGeoTrack &trk = event[t];

    if (opt.Contains("xy"))
    {
      trk.vx = xy(0);
      trk.vy = xy(1);
    }

    if (opt.Contains("z"))
    {
      trk.vz = rz(1);
    }
  }
  return;
}

void
FindVertexEastWest(geoTracks &event, TString opt)
{
  // Compute vertex position from each VTX arm.
  // Assign east (west) vertex to east (west) arm tracks.
  // Requires that tracks have already been fit.

  TVectorD xye = XYCenter(event, "east");
  TVectorD xyw = XYCenter(event, "west");
  TVectorD rze = ZVertexGLS(event, "east");
  TVectorD rzw = ZVertexGLS(event, "west");

  for (unsigned int t = 0; t < event.size(); t++)
  {
    SvxGeoTrack &trk = event[t];
    int arm = (trk.hits[0].x < 0.) ? 0 : 1;   // 0 = East, 1 = West.

    if (opt.Contains("xy"))
    {
      TVectorD vxy = arm ? xyw : xye;
      trk.vx = vxy(0);
      trk.vy = vxy(1);
    }

    if (opt.Contains("z"))
    {
      double vz = arm ? rzw(1) : rze(1);
      trk.vz = vz;
    }
  }
  return;
}

void
CalculateDCA(geoTracks &event, TGraphErrors *bc, TVectorD &xpz, TString opt)
{
  // Compute the 2-D distance of closest approach from reference point to track.
  // The reference point (x,y) is one of the following:
  // 1. A beam center if provided in the TGraph.
  // 2. The primary vertex saved in the track if no bc is provided.
  // The reference point z is provided as xpz(1) unless the option
  // "use_stored_vertex" is passed in.

  TVectorD xye(2);
  TVectorD xyw(2);
  double ze = 0, zw = 0;
  if (bc)
  {
    xye(0) = bc->GetX()[0];
    xye(1) = bc->GetY()[0];
    xyw(0) = bc->GetX()[1];
    xyw(1) = bc->GetY()[1];
  }

  // Compute radial DCA and save to track.
  for (unsigned int t = 0; t < event.size(); t++)
  {
    SvxGeoTrack &trk = event[t];
    int arm = (trk.hits[0].x < 0.) ? 0 : 1;   // 0 = East, 1 = West.
    TVectorD vxy = arm ? xyw : xye;

    if (!bc)
    {
      vxy(0) = trk.vx;
      vxy(1) = trk.vy;
    }

    // rotate vertex into primed frame. xp is distance from origin.
    // compute ds then rotate [0; ds] out of primed frame --> ipvec.
    TMatrixD R(2, 2);
    double c = cos(trk.phirot), s = sin(trk.phirot);
    R(0, 0) = c; R(0, 1) = -s;
    R(1, 0) = s; R(1, 1) = c;
    TMatrixD Rinv(R);
    Rinv(0, 1) = +s;
    Rinv(1, 0) = -s;

    TVectorD vxyp = Rinv * vxy;

    // IP vector in primed frame
    double mp = tan(trk.phi0 - trk.phirot);
    TVectorD ipp(2);
    ipp(1) = mp * vxyp(0) + trk.yp0 - vxyp(1);

    // Rotate back out of the primed frame
    TVectorD ipvec = R * ipp;

    TVectorD trknormal(2);
    double p = trk.phi0 + TMath::PiOver2();
    trknormal(0) = cos(p);
    trknormal(1) = sin(p);

    double vz = xpz(1);
    if (opt.Contains("use_stored_vertex"))
    {
      vz = trk.vz;
    }

    float xydca = trknormal * ipvec;
    float zdca  = trk.z0 + vxyp(0) / tan(trk.the0) - vz;

    trk.xydca = xydca;
    trk.zdca  = zdca;
  }

  return;
}

void
FitTracks(geoEvents &events, TGraphErrors *bc, TString opt)
{
  if (opt.Contains("find_vertex"))
  {
    Info("FitTracks()", "Vertex will be computed and (re)assigned to tracks.");

    // Prevent incompatible option. Vertex shouldn't exist yet.
    assert(!opt.Contains("fit_to_z_vertex"));
  }

  if (opt.Contains("fit_to_bc"))
  {
    assert(bc);
    Info("FitTracks()",
         "Track fits will include fixed beamcenter position (%.4f, %.4f).",
         bc->GetX()[1], bc->GetY()[1]);
  }
  if (opt.Contains("fit_to_z_vertex"))
  {
    // Check the first event for occurrences of trk.vz == zero (the initialized
    // value). Issue a warning if all zeros.
    int nzero = 0;
    for (unsigned int t = 0; t < events[0].size(); t++)
      if (fabs(events[0][t].vz) < 1e-15)
        nzero++;
    Info("FitTracks()", "Track fits will include stored z-vertex position.");

    if (nzero == (int)events[0].size())
      Warning("FitTracks()",
              "\n\nfit_to_z_vertex requested, but z vertex of all tracks in "
              "first event is ZERO!\n\n");

    if (opt.Contains("find_vertex"))
      Warning("FitTracks()",
              "Both fit_to_z_vertex and find_vertex options are requested.\n"
              "Z-vertex used for track fits will be overwritten!\n"
              "Combining these options is not recommended.");
  }

  if (opt.Contains("calc_dca"))
  {
    if (bc)
    {
      // Info("FitTracks()",
      //      "x-y DCA will be computed with respect to beam center.");
      // Info("FitTracks()", " E: (%.3f, %.3f)", bc->GetX()[0], bc->GetY()[0]);
      // Info("FitTracks()", " W: (%.3f, %.3f)", bc->GetX()[1], bc->GetY()[1]);
      Info("FitTracks()",
           "x-y DCA will use fixed beamcenter position (%.4f, %.4f) (east), "
           "(%.4f, %.4f) (west).",
           bc->GetX()[0], bc->GetY()[0], bc->GetX()[1], bc->GetY()[1]);
    }
    else
      Info("FitTracks()",
           "DCA will be computed with respect to stored primary vertex.");
    Info("FitTracks()",
         "z DCA will be computed with respect to global (E+W) z vertex.");
  }

  cout << Form("Fitting tracks in %lu events...", events.size()) << flush;

  for (unsigned int ev = 0; ev < events.size(); ev++)
  {
    for (unsigned int t = 0; t < events[ev].size(); t++)
    {
      // Do transverse fit first
      if (opt.Contains("fit_to_bc"))
        TrackFitT(events[ev][t], bc);
      else
        TrackFitT(events[ev][t]);

      // And longitudinal fit second (uses track data assigned in TrackFitT)
      if (opt.Contains("fit_to_z_vertex"))
        TrackFitL(events[ev][t], "fit_to_z_vertex");
      else
        TrackFitL(events[ev][t]);
    }

    if (opt.Contains("find_vertex"))
    {
      FindVertex(events[ev], "xyz");
    }

    // TODO: provide another option like find_ew_vertex ?
    // FindVertexEastWest(events[ev], "xyz");

    if (opt.Contains("calc_dca"))
    {
      // Calculate z vertex using both arms of VTX.
      // xpz(0) is the x' vertex position (essentially r).
      // xpz(1) is the vertex z position.
      TVectorD xpz = ZVertexGLS(events[ev], "");
      CalculateDCA(events[ev], bc, xpz, "");
    }
  }

  Printf("done.");
  return;
}

TVectorD
XYCenter(geoTracks &tracks, TString arm, int ntrk, TString opt)
{
  // Compute least-squares (x,y) center from ntrk tracks.
  // Used for primary vertex or beam center estimation.
  // Solves the linear system y0[i] = -m[i]*bc0 + bc1 for i..ntrk-1 tracks.
  // m[i] is the slope (tan(phi)), y0[i] is the y-intercept, and (bc0,bc1)
  // is the least-squares vertex or beam (x,y) position.
  // - arm can be "east", "west" or "" for both.
  // - ntrk is an optional upper limit on the number of tracks to use
  //   (to limit cpu time).
  // - opt: 1 or more of the following:
  //   * "" or empty (default).
  //   * "print" Print result.
  //   * "hitw"  Weight 4-hit tracks > 3-hit tracks in the fit.

  assert(ntrk <= (int)tracks.size());
  int n = ntrk > 0 ? ntrk : (int)tracks.size();
  TVectorD xy(2);

  if ((arm == "east" || arm == "west") && ntrk <= 0)
  {
    int nsub = 0;
    for (int i = 0; i < n; i++)
    {
      bool east = East(tracks[i].phi0);
      if ((arm == "east" && east) || (arm == "west" && !east))
        nsub++;
    }
    n = nsub;
  }

  if (n < 2)
  {
    xy(0) = -9999;
    xy(1) = -9999;
    return xy;
  }

  TMatrixD M(n, 2);  // "Design matrix" containing track slopes
  TVectorD y0(n);    // Vector of track y-intercept values
  TMatrixD L(n, n);  // Error (co)variance for track fit parameters y0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dy0 from track fits
  TMatrixD cov(2, 2); // Assigned in SolveGLS()

  int row = 0;
  for (unsigned int i = 0; i < tracks.size(); i++)
  {
    if (row == n)
      break;

    SvxGeoTrack trk = tracks[i];
    double mp = tan(trk.phi0 - trk.phirot);
    double yint = trk.yp0/(cos(trk.phirot) - mp*sin(trk.phirot));
    double phi = trk.phi0;
    bool east = East(phi);

    if (opt.Contains("hitw") && trk.nhits == 4)
      L(row, row) *= 0.5;

    if (arm == "")
    {
      M(row, 0) = -TMath::Tan(phi);
      M(row, 1) = 1;
      y0(row) = yint;
      row++;
    }
    else if ((arm == "east" && east) || (arm == "west" && !east))
    {
      M(row, 0) = -TMath::Tan(phi);
      M(row, 1) = 1;
      y0(row) = yint;
      row++;
    }
  }
  xy = SolveGLS(M, y0, L, cov);

  if (opt.Contains("print"))
  {
    Printf("%s x,y (%.3f +- %.3f, %.3f +- %.3f)",
           arm.Data(),
           xy(0), TMath::Sqrt(cov(0, 0)),
           xy(1), TMath::Sqrt(cov(1, 1)));
    cov.Print();
  }

  return xy;
}

TVectorD
ZVertexGLS(geoTracks &tracks, TString arm, int ntrk, TString opt)
{
  // Compute z vertex using tracks from east arm only ("east"),
  // west arm only ("west"), or both arms ("").
  // Tracks are not modified by calling this function.

  assert(ntrk <= (int)tracks.size());
  int n = ntrk > 0 ? ntrk : (int)tracks.size();
  TVectorD rz(2);

  if ((arm=="east" || arm=="west") && ntrk <= 0)
  {
    int nsub = 0;
    for (int i=0; i<n; i++)
    {
      bool east = East(tracks[i].phi0);
      if ((arm=="east" && east) || (arm=="west" && !east))
        nsub++;
    }
    n = nsub;
  }

  if (n < 2)
  {
    rz(0) = -9999;
    rz(1) = -9999;
    return rz;
  }

  TMatrixD M(n,2);   // "Design matrix" containing track slopes
  TVectorD z0(n);    // Vector of track y-intercept values
  TMatrixD L(n,n);   // Error (co)variance for track fit parameters z0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dz0 from track fits
  TMatrixD cov(2,2); // Assigned in SolveGLS()

  int row=0;
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    if (row==n)
      break;

    SvxGeoTrack trk = tracks[i];
    double slope = 1./tan(trk.the0);
    bool east = East(trk.phi0);

    if (opt.Contains("hitw") && trk.nhits == 4)
      L(row,row) *= 0.5;

    if (arm=="")
    {
      M(row, 0) = -slope;
      M(row, 1) = 1;
      z0(row) = trk.z0;
      row++;
    }
    else if ((arm=="east" && east) || (arm=="west" && !east))
    {
      M(row, 0) = -slope;
      M(row, 1) = 1;
      z0(row) = trk.z0;
      row++;
    }
  }
  rz = SolveGLS(M, z0, L, cov);

  if (opt.Contains("print"))
  {
    Printf("%s x,y (%.3f +- %.3f, %.3f +- %.3f)",
           arm.Data(),
           rz(0), sqrt(cov(0,0)),
           rz(1), sqrt(cov(1,1)));
    cov.Print();
  }

  return rz;
}

bool
East(double phi)
{
  return (phi > TMath::PiOver2() && phi < 3 * TMath::PiOver2());
}

double
ClusterXResolution(int layer)
{
  // Cluster position resolution estimates [cm] in layers 0-3.
  // Based on hardware specs + empirical studies.
  double res = 50e-4 / TMath::Sqrt(12.);
  double xres[4] = {res, 1.2 * res, 2 * res, 2.4 * res};

  return xres[layer];
}

double
ClusterZResolution(int layer)
{
  // Cluster position resolution estimates [cm] in layers 0-3.
  // Based on hardware specs + empirical studies.
  double res = 425e-4 / TMath::Sqrt(12.);
  double zres[4] = {res, res, 2 * res, 2 * res};

  return zres[layer];
}

double
ZVertexResolution()
{
  // Roughly based on std dev of z DCA distribution in data
  // TODO: Make this a settable value
  return 200e-4;
}

TVectorD
IPVec(TVectorD &a, TVectorD &n, TVectorD &p)
{
  // Compute impact parameter vector from point p to line x = a + tn
  // where
  // - x is a straight-line track trajectory
  // - a is a point on vector x (e.g. y-intercept point (0, y0))
  // - n is a unit vector directing x (e.g. (cos(phi), sin(phi)).
  return a - p - ((a - p) * n) * n;
}

TVectorD
IPVec(SvxGeoTrack &trk, TVectorD &p)
{
  double mp = tan(trk.phi0 - trk.phirot);
  // double y0 = trk.yp0 / (cos(trk.phirot) - mp * sin(trk.phirot));

  TVectorD a(2);
  a(0) = 0;
  a(1) = trk.yp0; //trk.yp0/cos(trk.phirot);
  TVectorD n(2);
  n(0) = cos(trk.phi0);
  n(1) = sin(trk.phi0);
  return IPVec(a, n, p);
}

TVectorD
RetrieveVertex(geoTracks &event, TString opt)
{
  // Get the event vertex for either the
  // East or West arm stored in the geoTrack
  TVectorD v(2);

  if (!opt.Contains("east") && !opt.Contains("west"))
  {
    Error("RetrieveVertex() in GLSFitter.h",
          "Option must contain either \"east\" or \"west\"");
    return v;
  }


  for (unsigned int t = 0; t < event.size(); t++)
  {
    int arm = (event[t].hits[0].x < 0.) ? 0 : 1;   // 0 = East, 1 = West.

    if (arm == 0 && opt.Contains("east"))
    {
      v(0) = event[t].vx;
      v(1) = event[t].vy;
      return v;
    }
    if (arm == 1 && opt.Contains("west"))
    {
      v(0) = event[t].vx;
      v(1) = event[t].vy;
      return v;
    }
  }

  //no vertex found
  return v;

}

#endif
