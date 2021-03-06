#ifndef __VERTEXFINDER_H__
#define __VERTEXFINDER_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"

TVectorD Vertex(geoTracks &event, TString arm, TString opt = "");
void FillVertexArrays(geoEvents &events, int minmult,
                      vecd &vxe, vecd &vye, vecd &vze,
                      vecd &vxw, vecd &vyw, vecd &vzw,
                      vecd &dvx, vecd &dvy, vecd &dvz,
                      int &n);
void FillVertexArrays(geoEvents &events, int minmult,
                      double *vxe, double *vye, double *vze,
                      double *vxw, double *vyw, double *vzw,
                      double *dvx, double *dvy, double *dvz,
                      int &nfilled);
void FillVertexHists(geoEvents &events, int minmult,
                     TH2D *he, TH2D *hw, TH1D **hdv);

// Code snippet for median calculation:
// root [0] double x[] = {1,2,3,4,5,6,7,8,9,10};
// root [1] double probs[] = {0.25, 0.5, 0.75};
// root [2] double quantiles[3] = {0};
// root [3] TMath::Quantiles(10, 3, x, quantiles, probs, false)
// root [4] quantiles
// (double [3]) { 3.250000e+00, 5.500000e+00, 7.750000e+00 }

TVectorD
Vertex(geoTracks &event, TString arm, TString opt)
{
  // Compute least-squares vertex from tracks.

  int n = event.size();
  assert(n>0);
  TVectorD xy = XYCenter(event, arm, -1, opt);

  TVectorD xyz(3);
  xyz(0) = xy(0);
  xyz(1) = xy(1);
  xyz(2) = ZVertex(event, arm);
  return xyz;
}

void
FillVertexArrays(geoEvents &events, int minmult,
                 vecd &vxe, vecd &vye, vecd &vze,
                 vecd &vxw, vecd &vyw, vecd &vzw,
                 vecd &dvx, vecd &dvy, vecd &dvz,
                 int &n)
{
  // Fill v{x,y,z}{e,w} vectors with vertex positions.
  // No checking is done that arrays are pre-allocated.

  n = 0;
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");

      vxe.push_back(ve(0));
      vye.push_back(ve(1));
      vze.push_back(ve(2));
      vxw.push_back(vw(0));
      vyw.push_back(vw(1));
      vzw.push_back(vw(2));

      dvx.push_back(vxw[n] - vxe[n]);
      dvy.push_back(vyw[n] - vye[n]);
      dvz.push_back(vzw[n] - vze[n]);

      n++;
    }
  }

  return;
}

void
FillVertexArrays(geoEvents &events, int minmult,
                 double *vxe, double *vye, double *vze,
                 double *vxw, double *vyw, double *vzw,
                 double *dvx, double *dvy, double *dvz,
                 int &n)
{
  // Fill v{x,y,z}{e,w} arrays with vertex positions.
  // No checking is done that arrays are pre-allocated.

  n = 0;
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");

      vxe[n] = ve(0);
      vye[n] = ve(1);
      vze[n] = ve(2);
      vxw[n] = vw(0);
      vyw[n] = vw(1);
      vzw[n] = vw(2);

      dvx[n] = vxw[n] - vxe[n];
      dvy[n] = vyw[n] - vye[n];
      dvz[n] = vzw[n] - vze[n];

      n++;
    }
  }

  return;
}

void FillVertexHists(geoEvents &events, int minmult,
                     TH2D *he, TH2D *hw, TH1D **hdv)
{
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");

      he->Fill(ve(0), ve(1));
      hw->Fill(vw(0), vw(1));

      for (int j=0; j<3; j++)
        hdv[j]->Fill(vw(j) - ve(j));
    }
  }

}

#endif
