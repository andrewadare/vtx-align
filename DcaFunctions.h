#ifndef __DCAFUNCTIONS_H__
#define __DCAFUNCTIONS_H__

#include "GLSFitter.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGeoManager.h>
#include <TProfile.h>

typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

TVectorD IPVec(TVectorD &a, TVectorD &n, TVectorD &p);
TVectorD IPVec(SvxGeoTrack &t, TVectorD &p);
TGraph *DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=10000);
TGraph *DcaDist(geoTracks &tracks, TVectorD &bc, TString arm,
                TH1D *hr=0, int ntracks=10000);
TGraph *DcaDist(geoEvents &events, TVectorD &bc, TString arm,
                TH1D *hr, int ntracks=10000);
// TH2D *DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
//                const char *name, const char *title);
// TH2D *DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
//                const char *name, const char *title);
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
DcaDist(TFile *f, TNtuple *t, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  // This function is no longer being used.
  TTreeReader r(t->GetName(), f);
  TTreeReaderValue<float> ty0(r, "y0");
  TTreeReaderValue<float> phi(r, "phi");
  TGraph *g = new TGraph();
  g->SetMarkerStyle(kFullDotMedium);

  int i=0;
  while (r.Next())
  {
    double m = TMath::Tan(*phi);
    if ((arm=="east" && East(*phi)) || (arm=="west" && !East(*phi)))
    {
      TVectorD a(2); a(1) = *ty0;
      TVectorD n(2); n(0) = TMath::Cos(*phi); n(1) = TMath::Sin(*phi);
      TVectorD d = IPVec(a,n,bc);  // d = a - bc - ((a - bc)*n)*n;

      if (i<ntracks)
        g->SetPoint(i, d(0), d(1));
      if (hr)
        hr->Fill(TMath::Sqrt(d*d));

      i++;
    }
  }
  return g;
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

TGraph *
DcaDist(geoEvents &events, TVectorD &bc, TString arm, TH1D *hr, int ntracks)
{
  int i = 0;
  TGraph *g = new TGraph(ntracks);
  g->SetMarkerStyle(kFullDotMedium);

  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
    {
      SvxGeoTrack trk = events[ev][t];
      double phi = trk.phi0;
      if ((arm=="east" && East(phi)) || (arm=="west" && !East(phi)))
      {
        TVectorD a(2); a(1) = trk.vy;
        TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
        TVectorD d = IPVec(a,n,bc);  // d = a - bc - ((a - bc)*n)*n;

        if (i<ntracks)
        {
          g->SetPoint(i, d(0), d(1));
          i++;
        }

        if (hr)
          hr->Fill(TMath::Sqrt(d*d));
      }
    }

  return g;
}

// TH2D *
// DcaVsPhi(geoTracks &tracks, TVectorD &bce, TVectorD &bcw,
//          const char *name, const char *title)
// {
//   TH2D *h = new TH2D(name, title, 100, 0, TMath::TwoPi(), 100, -0.0, 0.1);

//   for (unsigned int i=0; i<tracks.size(); i++)
//   {
//     double phi = tracks[i].phi0;
//     TVectorD a(2); a(1) = tracks[i].vy;
//     TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
//     TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
//     h->Fill(fmod(TMath::PiOver2()+phi,TMath::TwoPi()), TMath::Sqrt(d*d));
//   }

//   return h;
// }

// TH2D *
// DcaVsPhi(geoEvents &events, TVectorD &bce, TVectorD &bcw,
//          const char *name, const char *title)
// {
//   TH2D *h = new TH2D(name, title, 100, 0, TMath::TwoPi(), 100, -0.0, 0.1);

//   for (unsigned int ev=0; ev<events.size(); ev++)
//     for (unsigned int t=0; t<events[ev].size(); t++)
//     {
//       SvxGeoTrack trk = events[ev][t];
//       double phi = trk.phi0;
//       TVectorD a(2); a(1) = trk.vy;
//       TVectorD n(2); n(0) = TMath::Cos(phi); n(1) = TMath::Sin(phi);
//       TVectorD d = IPVec(a,n, East(phi) ? bce : bcw);
//       h->Fill(fmod(TMath::PiOver2()+phi,TMath::TwoPi()), TMath::Sqrt(d*d));
//     }

//   return h;
// }

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
