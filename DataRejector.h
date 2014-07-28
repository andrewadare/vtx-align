#ifndef __DATAREJECTOR_H__
#define __DATAREJECTOR_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"
#include "DcaFunctions.h"
#include "VtxIO.h"

void
FilterData(geoEvents &events,
           double vertexprobmin = 0.01,
           double vertexprobmax = 0.99,
           double maxdca = 0.2,
           double maxres_s = 0.1,
           double maxres_z = 0.1,
           TNtuple* ntuple = 0);

void
FilterData(geoEvents &events,
           double vertexprobmin,
           double vertexprobmax,
           double maxdca,
           double maxres_s,
           double maxres_z,
           TNtuple* ntuple)
{

  Printf("Event multiplicity scan...");
  int minmult = 20;
  int N = 0; // Pre-allocation size
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ne=0, nw=0;
    int mult = events[ev].size();
    if (mult >= minmult)
      N++;
    for (int t=0; t<mult; t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (East(trk.phi0)) ne++;
      else nw++;
    }
  }

  Printf("Allocating vertex arrays...");
  int nfilled = 0;
  double vxe[N], vye[N], vze[N];  // Vertex from east arm
  double vxw[N], vyw[N], vzw[N];  // Vertex from west arm
  double dvx[N], dvy[N], dvz[N];  // West - east offsets

  Printf("Filling vertex arrays...");
  FillVertexArrays(events, minmult,
                   vxe, vye, vze, vxw, vyw, vzw, dvx, dvy, dvz,
                   nfilled);

  Printf("Computing quantiles for {e,w}{x,y}");
  const int nq = 5;
  double probs[nq] = {vertexprobmin, 0.32, 0.50, 0.68, vertexprobmax};
  double qxe[nq] = {0}, qxw[nq] = {0}, qye[nq] = {0}, qyw[nq] = {0};
  TMath::Quantiles(nfilled, nq, vxe, qxe, probs, false);
  TMath::Quantiles(nfilled, nq, vxw, qxw, probs, false);
  TMath::Quantiles(nfilled, nq, vye, qye, probs, false);
  TMath::Quantiles(nfilled, nq, vyw, qyw, probs, false);

  TVectorD bce(2); bce(0) = qxe[2]; bce(1) = qye[2];
  TVectorD bcw(2); bcw(0) = qxw[2]; bcw(1) = qyw[2];

  Printf("Rejecting outliers in event loop...");
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    bool keep = true;
    int ntrk = events[ev].size();

    // Reject outlying vertices (box cut for now)
    if (ntrk > minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");

      keep = false;
      if (ve(0) > qxe[0] && ve(0) < qxe[nq-1]) // East x vertex
        if (ve(1) > qye[0] && ve(1) < qye[nq-1]) // East y vertex
          if (vw(0) > qxw[0] && vw(0) < qxw[nq-1]) // West x vertex
            if (vw(1) > qyw[0] && vw(1) < qyw[nq-1]) // West y vertex
              keep = true;
    }

    if (keep)
      for (int t=0; t<ntrk; t++)
      {
        SvxGeoTrack trk = events[ev][t];
        double phi = trk.phi0;
        TVectorD d = IPVec(trk, East(phi) ? bce : bcw);

        // DCA outlier cut
        if (TMath::Sqrt(d*d) > maxdca)
          continue;

        // Residual magnitude cut
        bool reject = false;
        for (int j=0; j<trk.nhits; j++)
        {
          SvxGeoHit hit = trk.GetHit(j);
          if (TMath::Abs(hit.ds) > maxres_s || TMath::Abs(hit.dz) > maxres_z)
            reject = true;
        }
        if (!reject)
          FillNTuple(trk, ntuple, ev);
      }
  }

  // double rejectFrac = 1.0 - (float)b.size()/a.size();
  // Printf("%.3f%% of tracks rejected by beam center DCA cut of %.0f um"
  //        " and residual outlier cut of %.0f um.\n\n",
  //        100*rejectFrac, 1e4*maxdca, 1e4*resCut);

  return;
}

#endif
