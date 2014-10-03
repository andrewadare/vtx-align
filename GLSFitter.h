#ifndef __GLSFITTER_H__
#define __GLSFITTER_H__

#include "VtxAlignBase.h"
#include <TDecompSVD.h>

using namespace std;
typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L);
TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L, TMatrixD &cov);
void TrackFitL(SvxGeoTrack &gt); // Longitudinal component --> z0, theta
void TrackFitT(SvxGeoTrack &gt); // Transverse component --> y0', phi
void FitTracks(geoTracks &tracks, TGraphErrors *bc = 0);
void FitTracks(geoEvents &events, TGraphErrors *bc = 0, TString opt = "");
void RadialDCA(geoTracks &event, TGraphErrors *bc, TString opt);

// It is questionable whether the following functions belong in this file.
bool East(double phi);
double ClusterXResolution(int layer);
double ClusterZResolution(int layer);
TVectorD XYCenter(geoTracks &event, TString arm, int ntrk = -1, TString opt="");
TVectorD IPVec(TVectorD &a, TVectorD &n, TVectorD &p);
TVectorD IPVec(SvxGeoTrack &t, TVectorD &p);
double ZVertex(geoTracks &event, TString arm);


TVectorD
SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L, TMatrixD &cov)
{
  // Simple generalized linear least-squares fitter.
  // Solve y(X) = beta'*X + error(L) for beta (vector of regression coefs.)
  // L is the inverse covariance matrix for y.
  // Least-squares solution involves solving X' L X beta = X' L y

  double tol = 1.e-15;
  TMatrixD XT(X); XT.T();
  TMatrixD A = XT*L*X;
  TVectorD b = XT*L*y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();
  TMatrixD V  = svd.GetV();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  int n = s.GetNrows();
  TMatrixD Sd(n,n);
  for (int i=0; i<n; i++)
    Sd(i,i) = s(i)>tol ? 1./s(i) : 0.;

  // Result
  TVectorD beta = V * Sd * UT * b;

  // and covariance matrix
  V *= Sd;
  TMatrixD C(V,TMatrixD::kMultTranspose,V);
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
  TMatrixD A = XT*L*X;
  TVectorD b = XT*L*y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  TMatrixD Sd(s.GetNrows(), s.GetNrows());
  for (int i=0; i<s.GetNrows(); i++)
    Sd(i,i) = s(i)>0 ? 1./s(i) : 0.;

  TVectorD beta = svd.GetV() * Sd * UT * b;

  return beta;
}

void
TrackFitL(SvxGeoTrack &gt)
{
  // Longitudinal/polar component
  // Fit track using a straight line z(r) = z0 + c*r.
  // Assigns residuals, z-intercept, and polar angle.

  // m = # measurements; n = # parameters.
  int m = gt.nhits, n = 2;
  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);

  for (int ihit=0; ihit < m; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    X(ihit,0) = 1;
    X(ihit,1) = TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);

    // Expecting that hit.zsigma is z_size in integer units 1,2,3....
    double zsigma = hit.zsigma * ClusterZResolution(hit.layer);
    Cinv(ihit,ihit) = (zsigma > 0) ? 1./zsigma : 1.;
    y(ihit) = hit.z;
  }

  TVectorD beta = SolveGLS(X, y, Cinv);
  double z0 = beta(0), c = beta(1);

  for (int ihit=0; ihit<m; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    double zproj = z0 + c*TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);

    gt.hits[ihit].dz = zproj - hit.z;
  }

  gt.z0 = z0;
  gt.the0 = TMath::ATan2(1.0, c);

  return;
}

void
TrackFitT(SvxGeoTrack &gt)
{
  // Transverse/azimuthal component
  // Perform straight-line fit y' = m'*x' + b' after rotating points
  // to approximately lie along the x axis. Then by construction,
  // x' ~ r, m' ~ 0, and the error is (mostly) in the y' = phi direction.
  // Assign s = r*phi residual and azimuthal angle and intercept info.

  int m = gt.nhits;
  int n = 2;
  TMatrixD points(n, m); // Datapoints. Columns are x',y' pairs
  TMatrixD X(m, n);      // Column 0 is 1's, column 1 is x' (~r) hit coords.
  TMatrixD Cinv(m, m);   // Inverse covariance matrix. Currently diagonal.
  TVectorD y(m);         // Dependent variables y'.
  double phirot = 0;     // <phi_cluster> used to rotate clusters

  int outermostLayer = 0;
  for (int ihit=0; ihit<m; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    points(0,ihit) = hit.x;
    points(1,ihit) = hit.y;

    X(ihit,0) = 1;
    // Note that hit.xsigma is expected to be in integer units 1,2,3....
    double xsigma = hit.xsigma * ClusterXResolution(hit.layer);
    Cinv(ihit,ihit) = (xsigma > 0) ? 1./xsigma : 1.;

    // Compute rotation angle as phi of outermost hit
    if (hit.layer > outermostLayer)
    {
      outermostLayer = hit.layer;
      phirot = TMath::ATan2(hit.y, hit.x);
    }
  }
  if (phirot < 0)
    phirot += TMath::TwoPi();

  // Rotate x, y by -phirot so error is approximately in y' direction only.
  // In rotated frame (xp, yp) = (cx + sy, cy - sx)
  TMatrixD R(2,2);
  double c = TMath::Cos(phirot), s = TMath::Sin(phirot);
  R(0,0) =  c; R(0,1) = s;
  R(1,0) = -s; R(1,1) = c;

  points = R * points;

  TMatrixDColumn(X, 1) = TMatrixDRow(points, 0);
  y = TMatrixDRow(points, 1);
  // Fit track to get beta prime = [b', m'] in rotated system.
  TVectorD betap = SolveGLS(X,y,Cinv);
  double bp = betap(0), mp = betap(1); // Intercept and slope, both small.

  // // Rotate back to get y-intercept b and slope m of track.
  // double y0 = bp / (c - mp*s);
  // double slope = (s + mp*c) / (c - mp*s);

  double phi = phirot + TMath::ATan(mp);
  if (phi < 0)
    phi += TMath::TwoPi();

  // Assign residuals to hits
  for (int ihit=0; ihit<gt.nhits; ihit++)
  {
    double xp  = points(0, ihit); // x' in rotated system--always positive
    assert(xp > 0);

    double ds = mp*xp + bp - points(1, ihit); // projected - measured y'
    gt.hits[ihit].ds = ds;
  }

  // Assign track variables
  gt.phi0 = phi;
  gt.phirot = phirot;
  gt.yp0 = bp;

  return;
}
// void
// TrackFitT(SvxGeoTrack &gt, double *pars, TGraphErrors * /*bc*/)
// {
//   // Perform straight-line fit y' = m'*x' + b' after rotating points
//   // to approximately lie along the x axis. Then by construction,
//   // x' ~ r, m' ~ 0, and the error is (mostly) in the y' = phi direction.
//   // Assign s = r*phi residual to gt.hits[i].ds.
//   // Optionally put [y0, phi] from fit in pars array.
//   // Use beam center info if provided.

//   // int m = (bc) ? gt.nhits + 1 : gt.nhits;
//   int m = gt.nhits;
//   int n = 2;
//   TMatrixD points(n, m); // Datapoints. Columns are x',y' pairs
//   TMatrixD X(m, n);      // Column 0 is 1's, column 1 is x' (~r) hit coords.
//   TMatrixD Cinv(m, m);   // Inverse covariance matrix. Currently diagonal.
//   TVectorD y(m);         // Dependent variables y'.
//   double phirot = 0;     // <phi_cluster> used to rotate clusters

//   for (int ihit=0; ihit<m; ihit++)
//   {
//     SvxGeoHit hit = gt.GetHit(ihit);
//     points(0,ihit) = hit.x;
//     points(1,ihit) = hit.y;

//     X(ihit,0) = 1;
//     // Note that hit.xsigma is expected to be in integer units 1,2,3....
//     double xsigma = hit.xsigma * ClusterXResolution(hit.layer);
//     Cinv(ihit,ihit) = (xsigma > 0) ? 1./xsigma : 1.;
//     phirot += 1./m * TMath::ATan2(hit.y, hit.x);
//   }

//   // Rotate x, y by -phirot so error is approximately in y' direction only.
//   // In rotated frame (xp, yp) = (cx + sy, cy - sx)
//   TMatrixD R(2,2);
//   double c = TMath::Cos(phirot), s = TMath::Sin(phirot);
//   R(0,0) =  c; R(0,1) = s;
//   R(1,0) = -s; R(1,1) = c;

//   points = R * points;

//   TMatrixDColumn(X, 1) = TMatrixDRow(points, 0);
//   y = TMatrixDRow(points, 1);
//   // Fit track to get beta prime = [b', m'] in rotated system.
//   TVectorD betap = SolveGLS(X,y,Cinv);
//   double bp = betap(0), mp = betap(1); // Intercept and slope, both small.

//   // Rotate back to get y-intercept b and slope m of track.
//   double y0 = bp / (c - mp*s);
//   double slope = (s + mp*c) / (c - mp*s);

//   double phi = phirot + TMath::ATan(mp);
//   if (phi < 0)
//     phi += TMath::TwoPi();
//   if (pars)
//   {
//     pars[0] = y0;
//     pars[1] = phi;
//   }

//   // Assign residuals to hits
//   for (int ihit=0; ihit<gt.nhits; ihit++)
//   {
//     double xp  = points(0, ihit); // x' in rotated system--always positive
//     assert(xp > 0);

//     double ds = mp*xp + bp - points(1, ihit); // projected - measured y'
//     gt.hits[ihit].ds = ds;
//   }

//   // Assign rotation angle and y-intercept of rotated track
//   gt.yp0 = bp;
//   gt.phirot = phirot;

//   return;
// }

// void
// ZeroFieldResiduals(SvxGeoTrack &gt,
//                    double *pars, /* y0, z0, yslope, zslope */
//                    TGraphErrors *bc /* bcx, bcy, sigma */)
// {
//   if (gt.nhits < 1)
//   {
//     Printf("ZeroFieldResiduals(): No hits in track. Skipping.");
//     return;
//   }

//   double zpars[2] = {0}, spars[2] = {0};
//   TrackFitL(gt, zpars);
//   TrackFitT(gt, spars, bc);

//   pars[0] = spars[0]; // y0 (y-intercept)
//   pars[1] = zpars[0]; // z0 (z-intercept)
//   pars[2] = spars[1]; // slope in transverse (y vs x) plane
//   pars[3] = zpars[1]; // slope in longitudinal (z vs r) plane

//   return;
// }

void
FitTracks(geoTracks &tracks, TGraphErrors * /*bc*/)
{
  cout << Form("Fitting %lu tracks...", tracks.size()) << flush;

  for (unsigned int i=0; i<tracks.size(); i++)
  {
    TrackFitL(tracks[i]);
    TrackFitT(tracks[i]);
  }

  Printf("done.");
  return;
}

void
RadialDCA(geoTracks &event, TGraphErrors *bc, TString opt)
{

  // Compute the 2-D distance of closest approach from reference point to track.
  // The reference point (x,y) is one of the following:
  // 1. A primary vertex computed here (compute_vertex). Also assigns vx,vy (!)
  // 2. A primary vertex already saved with the track (use_stored_vertex)
  // 3. A beam center provided in the TGraph (opt = "")
  // 4. Fall-through case: DCA will be computed w.r.t. (0,0).

  if (opt.Contains("no_dca"))
    return;

  TVectorD xye(2);
  TVectorD xyw(2);
  if (opt.Contains("compute_vertex"))
  {
    // Compute x,y vertex position from each VTX arm.
    // Requires that tracks have already been fit.
    xye = XYCenter(event, "east");
    xyw = XYCenter(event, "west");
  }
  else if (bc)
  {
    xye(0) = bc->GetX()[0];
    xye(1) = bc->GetY()[0];
    xyw(0) = bc->GetX()[1];
    xyw(1) = bc->GetY()[1];
  }
  else
  {
    // xye, xyw are already initialized to zero; this is just for clarity.
    xye(0) = 0.;
    xye(1) = 0.;
    xyw(0) = 0.;
    xyw(1) = 0.;
  }

  // Compute radial DCA and save to track.
  for (unsigned int t=0; t<event.size(); t++)
  {
    SvxGeoTrack &trk = event[t];
    int arm = (trk.hits[0].x < 0.) ? 0 : 1;   // 0 = East, 1 = West.
    TVectorD vxy = arm ? xyw : xye;

    // Printf("gen,calc x: (%.0f, %.0f) | gen,calc y: (%.0f, %.0f)",
    //        1e4*trk.vx, 1e4*vxy(0), 1e4*trk.vy, 1e4*vxy(1));

    if (opt.Contains("use_stored_vertex"))
    {
      vxy(0) = trk.vx;
      vxy(1) = trk.vy;
    }
    // TVectorD ipvec = IPVec(trk, vxy);         // Impact parameter vector

    ///////////////
    // rotate vertex into primed frame. xp is distance from origin.
    // compute ds then rotate [0; ds] out of primed frame --> ipvec.
    TMatrixD R(2,2);
    double c = cos(trk.phirot), s = sin(trk.phirot);
    R(0,0) = c; R(0,1) = -s;
    R(1,0) = s; R(1,1) = c;
    TMatrixD Rinv(R);
    Rinv(0,1) = +s;
    Rinv(1,0) = -s;

    TVectorD vxyp = Rinv * vxy;

    // IP vector in primed frame
    double mp = tan(trk.phi0 - trk.phirot);
    TVectorD ipp(2);
    ipp(1) = mp*vxyp(0) + trk.yp0 - vxyp(1);

    // Rotate back out of the primed frame

    TVectorD ipvec = R*ipp;

    TVectorD trknormal(2);
    double p = trk.phi0 + TMath::PiOver2();
    trknormal(0) = cos(p);
    trknormal(1) = sin(p);

    float dca = trknormal * ipvec;

    // Assignment
    trk.xydca = dca;
    if (opt.Contains("compute_vertex"))
    {
      trk.vx = vxy(0);
      trk.vy = vxy(1);
    }
  }

  return;
}

void
FitTracks(geoEvents &events, TGraphErrors *bc, TString opt)
{
  if (opt.Contains("compute_vertex"))
    Info("", "Vertex and DCA will be computed after track fits "
         "and (re)assigned to tracks.");
  else if (opt.Contains("use_stored_vertex"))
    Info("", "DCA will be computed after track fits using stored vertices.");
  else if (bc)
    Info("", "DCA will be computed with respect to beam center.");
  else if (opt.Contains("no_dca"))
    Info("", "Skipping DCA calculation.");

  cout << Form("Fitting tracks in %lu events...", events.size()) << flush;

  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      TrackFitL(events[ev][t]);
      TrackFitT(events[ev][t]);
    }

    RadialDCA(events[ev], bc, opt);
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
    xy(0) = -9999;
    xy(1) = -9999;
    return xy;
  }

  TMatrixD M(n,2);   // "Design matrix" containing track slopes
  TVectorD y0(n);    // Vector of track y-intercept values
  TMatrixD L(n,n);   // Error (co)variance for track fit parameters y0, m
  L.UnitMatrix();
  L *= 0.01;         // TODO: Get mean dm, dy0 from track fits
  TMatrixD cov(2,2); // Assigned in SolveGLS()

  int row=0;
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    if (row==n)
      break;

    SvxGeoTrack trk = tracks[i];
    double mp = tan(trk.phi0 - trk.phirot);
    double yint = trk.yp0/(cos(trk.phirot) - mp*sin(trk.phirot));
    // double yint = trk.yp0; // CHECK!! NEED TO ROTATE Y INTERCEPT?
    double phi = trk.phi0;
    bool east = East(phi);

    if (opt.Contains("hitw") && trk.nhits == 4)
      L(row,row) *= 0.5;

    if (arm=="")
    {
      M(row, 0) = -TMath::Tan(phi);
      M(row, 1) = 1;
      y0(row) = yint;
      row++;
    }
    else if ((arm=="east" && east) || (arm=="west" && !east))
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
           xy(0), TMath::Sqrt(cov(0,0)),
           xy(1), TMath::Sqrt(cov(1,1)));
    cov.Print();
  }

  return xy;
}

bool
East(double phi)
{
  return (phi > TMath::PiOver2() && phi < 3*TMath::PiOver2());
}

double
ClusterXResolution(int layer)
{
  // Cluster position resolution estimates [cm] in layers 0-3.
  // Based on hardware specs + empirical studies.
  double res = 50e-4/TMath::Sqrt(12.);
  double xres[4] = {res, 1.2*res, 2*res, 2.4*res};

  return xres[layer];
}

double
ClusterZResolution(int layer)
{
  // Cluster position resolution estimates [cm] in layers 0-3.
  // Based on hardware specs + empirical studies.
  double res = 425e-4/TMath::Sqrt(12.);
  double zres[4] = {res, res, 2*res, 2*res};

  return zres[layer];
}

TVectorD
IPVec(TVectorD &a, TVectorD &n, TVectorD &p)
{
  // Compute impact parameter vector from point p to line x = a + tn
  // where
  // - x is a straight-line track trajectory
  // - a is a point on vector x (e.g. y-intercept point (0, y0))
  // - n is a unit vector directing x (e.g. (cos(phi), sin(phi)).
  return a - p - ((a - p)*n)*n;
}

TVectorD
IPVec(SvxGeoTrack &t, TVectorD &p)
{
  // TVectorD a(2); a(1) = t.yp0; // CHECK!! NEED TO ROTATE Y INTERCEPT?

  double mp = tan(t.phi0 - t.phirot);
  double y0 = t.yp0/(cos(t.phirot) - mp*sin(t.phirot));

  TVectorD a(2);
  a(0) = 0;
  a(1) = t.yp0; //t.yp0/cos(t.phirot);
  TVectorD n(2);
  n(0) = cos(t.phi0);
  n(1) = sin(t.phi0);
  return IPVec(a,n,p);
}

double
ZVertex(geoTracks &tracks, TString arm)
{
  int nz = 0;
  const int maxnz = tracks.size();
  assert(maxnz>0);
  double zs[maxnz];
  double probs[1] = {0.5}; // For median = 50% quantile
  double quantiles[1] = {0};

  for (unsigned int i=0; i<tracks.size(); i++)
  {
    bool east = East(tracks[i].phi0);
    if (arm=="")
      zs[nz++] = tracks[i].vz;
    else if ((arm=="east" && east) || (arm=="west" && !east))
      zs[nz++] = tracks[i].vz;
  }
  TMath::Quantiles(nz, 1, zs, quantiles, probs, false);

  return quantiles[0];
}



#endif
