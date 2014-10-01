#ifndef __DCAFUNCTIONS_H__
#define __DCAFUNCTIONS_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"

TGraph *DcaDist(geoTracks &tracks, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=10000);
void DcaDist(geoEvents &events, TString arm, TH2D *hxy, TH1D *hr);
void DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
              TH2D *hist, TProfile *prof);
void DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
              TH2D *hist, TProfile *prof);

TGraph *
DcaDist(geoTracks &tracks, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  TGraph *g = new TGraph();
  g->SetMarkerStyle(kFullDotMedium);
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    Printf("**************************DEPRECATE****************************");
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
