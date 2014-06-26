
#include "SvxGeoTrack.h"
#include <TDecompSVD.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TNtuple.h>

using namespace std;

TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L);
TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L, TMatrixD &cov);
void TrackFitZResid(SvxGeoTrack &gt, double *pars = 0);
void TrackFitSResid(SvxGeoTrack &gt, double *pars = 0);
void ZeroFieldResiduals(SvxGeoTrack &gt, double *pars /* y0, z0, phi, theta */);
void Residuals(SvxGeoTrack &tt, SvxGeoTrack &mt, TNtuple *t);

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
ZeroFieldResiduals(SvxGeoTrack &gt, double *pars /* y0, z0, yslope, zslope */)
{
  if (gt.nhits < 1)
  {
    Printf("ZeroFieldResiduals(): No hits in track. Skipping.");
    return;
  }

  double zpars[2] = {0}, spars[2] = {0};
  TrackFitZResid(gt, zpars);
  TrackFitSResid(gt, spars);

  pars[0] = spars[0]; // y0 (y-intercept)
  pars[1] = zpars[0]; // z0 (z-intercept)
  pars[2] = spars[1]; // slope in transverse (y vs x) plane
  pars[3] = zpars[1]; // slope in longitudinal (z vs r) plane

  return;
}

void
TrackFitZResid(SvxGeoTrack &gt, double *pars)
{
  // Perform straight-line fit z(r) = z0 + c*r.
  // Assign residuals to gt.hits[i].dz. Put [z0, theta] from fit in pars.
  int m = gt.nhits, n = 2;
  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);

  for (int ihit=0; ihit < m; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    X(ihit,0) = 1;
    X(ihit,1) = TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);
    Cinv(ihit,ihit) = (hit.zsigma > 0) ? 1./hit.zsigma : 1.;
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

  if (pars)
  {
    pars[0] = z0;
    pars[1] = TMath::ATan2(1.0, c);
  }

  return;
}

void
TrackFitSResid(SvxGeoTrack &gt, double *pars)
{
  // Perform straight-line fit y' = m'*x' + b' after rotating points
  // to approximately lie along the x axis. Then by construction,
  // x' ~ r, m' ~ 0, and the error is (mostly) in the y' = phi direction.
  // Assign s = r*phi residual to gt.hits[i].ds.
  // Optionally put [y0, phi] from fit in pars array.

  int m = gt.nhits, n = 2;
  TMatrixD points(n, m); // Columns are x',y' pairs
  TMatrixD X(m, n);      // Column 0 is 1's, column 1 is x' (~r) hit coords.
  TMatrixD Cinv(m, m);   // Inverse covariance matrix. Currently diagonal.
  TVectorD y(m);         // Dependent variables y'.
  double phirot = 0;

  for (int ihit=0; ihit<m; ihit++)
  {
    SvxGeoHit t = gt.GetHit(ihit);
    points(0,ihit) = t.x;
    points(1,ihit) = t.y;
    X(ihit,0) = 1;
    Cinv(ihit,ihit) = (t.xsigma > 0) ? 1./t.xsigma : 1.;
    phirot += 1./m * TMath::ATan2(t.y, t.x);
  }

  // Rotate x, y by -phirot so error is approximately in y' direction only.
  // In rotated frame (xp, yp) = (cx + sy, cy - sx)
  TMatrixD R(2,2);
  double c = TMath::Cos(phirot), s = TMath::Sin(phirot);
  R(0,0) =  c; R(0,1) = s;
  R(1,0) = -s; R(1,1) = c;

  points = R * points;

  TMatrixDColumn(X, 1) = TMatrixDRow(points, 0);
  y = TMatrixDRow(points, 1);

  // Fit track to get [b', m'] in rotated system.
  TVectorD beta = SolveGLS(X,y,Cinv);
  double bp = beta(0), mp = beta(1);

  // Rotate back to get y-intercept b of track in x=0 plane
  // and phi angle of track.
  double y0  = bp / (c - mp*s);
  double phi = phirot + TMath::ATan(mp);
  if (phi < 0)
    phi += TMath::TwoPi();
  if (pars)
  {
    pars[0] = y0;
    pars[1] = phi;
  }

  for (int ihit=0; ihit<m; ihit++)
  {
    double x  = points(0, ihit);
    double ds = mp*x + bp - points(1, ihit); // projected - measured y'
    gt.hits[ihit].ds = (x < 0) ? -ds : +ds;
  }

  return;
}
