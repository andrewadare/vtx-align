#ifndef __DATAREJECTOR_H__
#define __DATAREJECTOR_H__

#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"
#include "DcaFunctions.h"
#include "VtxIO.h"

#include "TRandom3.h"

void FilterData(geoEvents &segevents,
                geoEvents &cntevents,
                TTree *segtree = 0,
                TTree *cnttree = 0,
                TGraphErrors *gbc = 0,
                double vertexprobmin = 0.01,
                double vertexprobmax = 0.99,
                double maxdca = 0.2,
                double maxres_s = 0.1,
                double maxres_z = 0.1,
                int maxclus = 10,
                int minmult = 10,
                int nhitsmin = 3,
                float frac4hit = -1);

bool RejectOutlierVtx(geoTracks &trks,
                      TVectorD &bce, TVectorD &ebce,
                      TVectorD &bcw, TVectorD &ebcw,
                      int minmult = 10);

bool RejectOutlierTrk(SvxGeoTrack trk,
                      float maxdca = 0.2,
                      float maxres_s = 0.1,
                      float maxres_z = 0.1,
                      int maxclus = 10,
                      int nhitsmin = 3);

void Fix4HitFrac(geoEvents &segeventsout, geoEvents &cnteventsout,
                 geoEvents &segeventsin, geoEvents &cnteventsin,
                 float frac = 0.5);


void
FilterData(geoEvents &segevents, geoEvents &cntevents,
           TTree *segtree, TTree *cnttree,
           TGraphErrors *gbc,
           double vertexprobmin, double vertexprobmax ,
           double maxdca,
           double maxres_s, double maxres_z,
           int maxclus, int minmult, int nhitsmin,
           float frac4hit)
{

  // This function filters data, keeping
  // svxsegments and svxcnt in sync
  // if cntevents is included
  assert(segtree);

  if (cntevents.size() > 0 && !cnttree)
  {
    Error("FilterData() in DataRejector.h",
          "cntevents is non-empty, but cnttree is is NULL");
    return;
  }

  //double check that geoEvents are the same size
  //else they aren't synced and this won't work
  if (cntevents.size() > 0 && segevents.size() != cntevents.size())
  {
    Error("FilterData() in DataRejector.h",
          "Number of segevents != cntevents! Nsegevents=%u, Ncntevents=%u",
          (unsigned int)segevents.size(),
          (unsigned int)cntevents.size());
    return;
  }

  // get the bc & uncertainties for rejecting outlying events
  TVectorD bce(2), bcw(2); //beamcenter values
  TVectorD ebce(2), ebcw(2);//beamcenter error values
  if (gbc)
  {
    // Use a fixed beam center

    //East
    bce(0) = gbc->GetX()[0]; bce(1) = gbc->GetY()[0];
    ebce(0) = 4 * gbc->GetErrorX(0); 
    ebce(1) = 4 * gbc->GetErrorY(0);

    //West
    bcw(0) = gbc->GetX()[1]; bcw(1) = gbc->GetY()[1];
    ebcw(0) = 4 * gbc->GetErrorX(1); 
    ebcw(1) = 4 * gbc->GetErrorY(1);

  }
  else
  {
    // Calculate the beam center

    // x-y vertex distributions
    Printf("Filling mult/vertex distributions...");
    double x0 = -1.0, y0 = -1.0, x1 = +1.0, y1 = +1.0;
    TH2D *hve = new TH2D(Form("%s_ve", segtree->GetName()),
                         Form("East vertex;x [cm];y [cm]"),
                         500, x0, x1, 500, y0, y1);
    TH2D *hvw = new TH2D(Form("%s_vw", segtree->GetName()),
                         Form("West vertex;x [cm];y [cm]"),
                         500, x0, x1, 500, y0, y1);
    for (unsigned int ev = 0; ev < segevents.size(); ev++)
    {
      int mult = segevents[ev].size();
      if (mult >= minmult)
      {
        TVectorD ve = RetrieveVertex(segevents[ev], "east");
        hve->Fill(ve(0), ve(1));

        TVectorD vw = RetrieveVertex(segevents[ev], "west");
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

    bce(0) = qxe[2]; bce(1) = qye[2];
    ebce(0) = 0.5 * (qxe[nq - 1] - qxe[0]);
    ebce(1) = 0.5 * (qye[nq - 1] - qye[0]);

    bcw(0) = qxw[2]; bcw(1) = qyw[2];
    ebcw(0) = 0.5 * (qxw[nq - 1] - qxw[0]);
    ebcw(1) = 0.5 * (qyw[nq - 1] - qyw[0]);

    Printf("  bcE=(%.4f+/-%.4f, %.4f+/-%.4f)", bce(0), ebce(0), bce(1), ebce(1));
    Printf("  bcW=(%.4f+/-%.4f, %.4f+/-%.4f)", bcw(0), ebcw(0), bcw(1), ebcw(1));


  }


  Printf("Rejecting outliers in event loop...");
  int nsegin = 0;
  int ncntin = 0;

  //temporary containers to hold accepted tracks
  geoEvents accsegevents;
  geoEvents acccntevents;

  for (unsigned int ev = 0; ev < segevents.size(); ev++)
  {
    nsegin += segevents.at(ev).size();
    if (cntevents.size() > 0)
      ncntin += cntevents.at(ev).size();

    bool rejectev = RejectOutlierVtx(segevents.at(ev),
                                     bce, ebce,
                                     bcw, ebcw,
                                     minmult);

    if (!rejectev)
    {
      geoTracks tmptrks;
      //loop over segments and reject outliers
      for (unsigned int itrk = 0; itrk < segevents.at(ev).size(); itrk++)
      {
        bool rejecttrk = RejectOutlierTrk(segevents[ev][itrk],
                                          maxdca, maxres_s, maxres_z,
                                          maxclus, nhitsmin);
        if (!rejecttrk)
        {
          tmptrks.push_back(segevents[ev][itrk]);
        }
      }

      if (tmptrks.size() > 0)
      {
        accsegevents.push_back(tmptrks);
        if (cntevents.size() > 0)
          acccntevents.push_back(cntevents[ev]);
      }
    }
  }

  //if desired, filter further to get a given percentage of 4hit segments
  if (frac4hit > 0)
  {

    geoEvents filtsegevents;
    geoEvents filtcntevents;

    Fix4HitFrac(filtsegevents, filtcntevents,
                accsegevents, acccntevents,
                frac4hit);

    FillTree(filtsegevents, segtree);
    if (cntevents.size() > 0)
      FillTree(filtcntevents, cnttree);
  }
  else
  {
    //write out all the desired tracks
    FillTree(accsegevents, segtree);
    if (cntevents.size() > 0)
      FillTree(acccntevents, cnttree);
  }

  double segRejectFrac = 1.0 - (float)segtree->GetEntries() / (float)nsegin;
  Printf("%.3f%% of segments rejected by beam center DCA cut of %.0f um"
         " and residual (s,z) outlier cut of (%.0f,%.0f) um.",
         100 * segRejectFrac, 1e4 * maxdca, 1e4 * maxres_s, 1e4 * maxres_z);

  if (cntevents.size() > 0)
  {
    double cntRejectFrac = 1.0 - (float)cnttree->GetEntries() / (float)ncntin;
    Printf("%.3f%% of svxcnt rejected by beam center DCA cut of %.0f um"
           " and residual (s,z) outlier cut of (%.0f,%.0f) um.",
           100 * cntRejectFrac, 1e4 * maxdca, 1e4 * maxres_s, 1e4 * maxres_z);
  }

  return;
}

bool
RejectOutlierVtx(geoTracks &trks,
                 TVectorD &bce, TVectorD &ebce,
                 TVectorD &bcw, TVectorD &ebcw,
                 int minmult)
{
  // Reject events based on Vtx
  // true = reject event

  if ((int)trks.size() < minmult)
  {
    return true;
  }

  //find the event east/west vertex
  TVectorD ve = RetrieveVertex(trks, "east");
  TVectorD vw = RetrieveVertex(trks, "west");

  double dxe = ve(0) - bce(0);
  double dye = ve(1) - bce(1);
  double dxw = vw(0) - bcw(0);
  double dyw = vw(1) - bcw(1);

  double ex = ebce(0);
  double ey = ebce(1);
  double wx = ebcw(0);
  double wy = ebcw(1);

  /*
  Printf("\n");
  Printf("  bcE=(%.4f+/-%.4f, %.4f+/-%.4f)", bce(0), ebce(0), bce(1), ebce(1));
  Printf("  bcW=(%.4f+/-%.4f, %.4f+/-%.4f)", bcw(0), ebcw(0), bcw(1), ebcw(1));
  Printf("   vE=(%.4f, %.4f)", ve(0), ve(1));
  Printf("   vW=(%.4f, %.4f)", vw(0), vw(1));
  Printf("   fE=%f", dxe * dxe / ex / ex + dye * dye / ey / ey);
  Printf("   fW=%f", dxw * dxw / wx / wx + dyw * dyw / wy / wy);
  */

  if (dxe * dxe / ex / ex + dye * dye / ey / ey < 1.0 &&
      dxw * dxw / wx / wx + dyw * dyw / wy / wy < 1.0)
  {
    return false;
  }
  else
  {
    //Printf("  FAILED RejectOutlierVtx()");
    return true;
  }

  return false;
  // keep = false;
  // if (ve(0) > qxe[0] && ve(0) < qxe[nq-1]) // East x vertex
  //   if (ve(1) > qye[0] && ve(1) < qye[nq-1]) // East y vertex
  //     if (vw(0) > qxw[0] && vw(0) < qxw[nq-1]) // West x vertex
  //       if (vw(1) > qyw[0] && vw(1) < qyw[nq-1]) // West y vertex
  //         keep = true;
}

bool
RejectOutlierTrk(SvxGeoTrack trk,
                 float maxdca, float maxres_s, float maxres_z, int maxclus,
                 int nhitsmin)
{
  // Reject tracks
  // true = reject

  // DCA outlier cut
  if (TMath::Abs(trk.xydca) > maxdca)
    return true;

  // Require at least nhitsmin clusters
  if (trk.nhits < nhitsmin)
    return true;

  // Residual magnitude cut
  for (int j = 0; j < trk.nhits; j++)
  {
    SvxGeoHit hit = trk.GetHit(j);
    if (TMath::Abs(hit.ds) > maxres_s || TMath::Abs(hit.dz) > maxres_z)
      return true;

    // Reject tracks with suspiciously large clusters
    if (hit.xsigma > maxclus || hit.zsigma > maxclus)
      return true;
  }

  return false;
}


void Fix4HitFrac(geoEvents &segeventsout, geoEvents &cnteventsout,
                 geoEvents &segeventsin, geoEvents &cnteventsin,
                 float frac)
{
  // Reject 3 hit tracks randomely until they make up the desired
  // fraction of overall tracks

  // check if we need to keep cnt tracks in sync
  bool synccnt = false;
  if (cnteventsin.size() > 0)
    synccnt = true;


  // count the number of 3/4 hit tracks
  int N3hit = 0;
  int N4hit = 0;
  for (unsigned int i = 0; i < segeventsin.size(); i++)
  {
    for (unsigned int j = 0; j < segeventsin.at(i).size(); j++)
    {
      if (segeventsin[i][j].nhits == 3)
        N3hit++;
      if (segeventsin[i][j].nhits == 4)
        N4hit++;
    }
  }
  Printf("Input segments have %.2f%% 4 hit tracks",
         (float)N4hit / (float)(N4hit + N3hit) * 100.);

  // calculate how many 3 hit tracks we want based on the number of 4 hits
  // n3 = (1 - frac) / frac * n4
  int N3hit_desired = N4hit * ((1 - frac) / frac);

  if (N3hit_desired > N3hit)
  {
    Info("Fix4HitFrac() in DataRejector.h",
         "Fraction of 4 hit track exceeds expectations.");
    segeventsout = segeventsin;
    cnteventsout = cnteventsin;
    return;
  }

  //segeventsout = segeventsin;
  //cnteventsout = cnteventsin;


  float keepfrac = (float)N3hit_desired / (float)N3hit;

  TRandom3 rnd;

  int N3hit_out = 0;
  int N4hit_out = 0;
  for (unsigned int i = 0; i < segeventsin.size(); i++)
  {
    geoTracks tmp;

    for (unsigned int j = 0; j < segeventsin.at(i).size(); j++)
    {

      if (segeventsin[i][j].nhits == 3)
      {
        if (rnd.Uniform(0, 1) < keepfrac)
        {
          tmp.push_back(segeventsin[i][j]);
          N3hit_out++;
        }
      }
      else
      {
        N4hit_out++;
      }

    }

    if (tmp.size() > 0)
    {
      segeventsout.push_back(tmp);
      if (synccnt)
        cnteventsout.push_back(cnteventsin[i]);
    }
  }


  if (!synccnt)
    cnteventsout = cnteventsin;

  Printf("Output segments have %.2f%% 4 hit tracks",
         (float)N4hit_out / (float)(N4hit_out + N3hit_out) * 100.);


  return;


}


#endif
