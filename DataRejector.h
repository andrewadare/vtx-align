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
           TNtuple *ntuple = 0);

void
FilterData(geoEvents &events,
           double vertexprobmin,
           double vertexprobmax,
           double maxdca,
           double maxres_s,
           double maxres_z,
           TNtuple *ntuple)
{
  int minmult = 10;

  // x-y vertex distributions
  Printf("Filling mult/vertex distributions...");
  double x0 = -1.0, y0 = -1.0, x1 = +1.0, y1 = +1.0;
  TH2D *hve = new TH2D(Form("%s_ve",ntuple->GetName()),
                       Form("East vertex;x [cm];y [cm]"),
                       500, x0, x1, 500, y0, y1);
  TH2D *hvw = new TH2D(Form("%s_vw",ntuple->GetName()),
                       Form("West vertex;x [cm];y [cm]"),
                       500, x0, x1, 500, y0, y1);
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int mult = events[ev].size();
    if (mult >= minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");
      hve->Fill(ve(0), ve(1));
      hvw->Fill(vw(0), vw(1));
    }
  }

  Printf("Computing quantiles for {e,w}{x,y}");
  const int nq = 5;
  double probs[nq] = {vertexprobmin, 0.32, 0.50, 0.68, vertexprobmax};
  double qxe[nq] = {0}, qxw[nq] = {0}, qye[nq] = {0}, qyw[nq] = {0};
  hve->ProjectionX()->GetQuantiles(nq, qxe, probs);
  hvw->ProjectionX()->GetQuantiles(nq, qxw, probs);
  hve->ProjectionY()->GetQuantiles(nq, qye, probs);
  hvw->ProjectionY()->GetQuantiles(nq, qyw, probs);

  TVectorD bce(2); bce(0) = qxe[2]; bce(1) = qye[2];
  TVectorD bcw(2); bcw(0) = qxw[2]; bcw(1) = qyw[2];

  Printf("Rejecting outliers in event loop...");
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    bool keep = true;
    int ntrk = events[ev].size();

    // Reject outlying vertices
    if (ntrk > minmult)
    {
      keep = false;

      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");
      double dxe = ve(0) - bce(0);
      double dye = ve(1) - bce(1);
      double dxw = vw(0) - bcw(0);
      double dyw = vw(1) - bcw(1);
      double ex = 0.5*(qxe[nq-1] - qxe[0]);
      double ey = 0.5*(qye[nq-1] - qye[0]);
      double wx = 0.5*(qxw[nq-1] - qxw[0]);
      double wy = 0.5*(qyw[nq-1] - qyw[0]);

      if (dxe*dxe/ex/ex + dye*dye/ey/ey < 1.0 &&
          dxw*dxw/wx/wx + dyw*dyw/wy/wy < 1.0)
      {
        keep = true;
      }
      // keep = false;
      // if (ve(0) > qxe[0] && ve(0) < qxe[nq-1]) // East x vertex
      //   if (ve(1) > qye[0] && ve(1) < qye[nq-1]) // East y vertex
      //     if (vw(0) > qxw[0] && vw(0) < qxw[nq-1]) // West x vertex
      //       if (vw(1) > qyw[0] && vw(1) < qyw[nq-1]) // West y vertex
      //         keep = true;
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
