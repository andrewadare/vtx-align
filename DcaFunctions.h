#ifndef __DCAFUNCTIONS_H__
#define __DCAFUNCTIONS_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"


TVectorD IPVec(TVectorD &a, TVectorD &n, TVectorD &p);
TVectorD IPVec(SvxGeoTrack &t, TVectorD &p);
// TGraph *DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm,
//                 TH1D *hr=0, int ntracks=10000);
TGraph *DcaDist(geoTracks &tracks, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=10000);
void DcaDist(geoEvents &events, TString arm, TH2D *hxy, TH1D *hr);
void DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
              TH2D *hist, TProfile *prof);
void DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
              TH2D *hist, TProfile *prof);

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
  TVectorD a(2); a(1) = t.vy;
  TVectorD n(2);
  n(0) = TMath::Cos(t.phi0);
  n(1) = TMath::Sin(t.phi0);
  return IPVec(a,n,p);
}

TGraph *
DcaDist(geoTracks &tracks, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  TGraph *g = new TGraph();
  g->SetMarkerStyle(kFullDotMedium);
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    double phi = tracks[i].phi0;
    if ((arm=="east" && East(phi)) || (arm=="west" && !East(phi)))
    {
      TVectorD a(2); a(1) = tracks[i].vy;
      TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
      TVectorD d = IPVec(a,n,bc);  // d = a - bc - ((a - bc)*n)*n;

      if (i<(unsigned int)ntracks)
        g->SetPoint(i, d(0), d(1));
      if (hr)
        hr->Fill(TMath::Sqrt(d*d));
    }
  }
  return g;
}

void
DcaDist(geoEvents &events, TString arm, TH2D *hxy, TH1D *hr)
{
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk < 10)
      continue;

    TVectorD vertex = Vertex(events.at(ev), arm);
    TVectorD vxy(2);
    vxy(0) = vertex(0);
    vxy(1) = vertex(1);

    for (int t=0; t<ntrk; t++)
    {
      SvxGeoTrack trk = events[ev][t];
      double phi = trk.phi0;
      TVectorD a(2); a(1) = trk.vy;
      TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);

      if ((arm=="east" && East(phi)) || (arm=="west" && !East(phi)))
      {
        TVectorD d = IPVec(a,n,vxy);

        hxy->Fill(d(0), d(1));

        if (hr)
          hr->Fill(TMath::Sqrt(d*d));
      }
    }
  }

  return;
}

void
DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
         TH2D *hist, TProfile *prof)
{
  int len = bce.GetNrows();
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    SvxGeoTrack trk = tracks[i];
    double phi = trk.phi0;
    TVectorD a(len);
    a(1) = trk.vy;
    TVectorD n(len);
    n(0) = TMath::Cos(phi);
    n(1) = TMath::Sin(phi);
    if (len==3)
      a(2) = East(phi) ? bce(2) : bcw(2); // So Delta z = 0 (computing 2-D dca)

    TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
    double dmag = TMath::Sqrt(d*d);
    double phir = fmod(TMath::PiOver2()+phi,TMath::TwoPi());
    hist->Fill(phir, dmag);
    prof->Fill(phir, dmag);
  }

  return;
}

void
DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
         TH2D *hist, TProfile *prof)
{
  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      double phi = trk.phi0;
      TVectorD a(2); a(1) = trk.vy;
      TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
      TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
      double dmag = TMath::Sqrt(d*d);
      double phir = fmod(TMath::PiOver2()+phi,TMath::TwoPi());
      hist->Fill(phir, dmag);
      prof->Fill(phir, dmag);
    }

  return;
}

#endif
