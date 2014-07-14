#include "UtilFns.h"
#include "VtxIO.h"
#include "GLSFitter.h"
#include "VertexFinder.h"

#include <TGeoManager.h>
#include <TFile.h>
#include <TNtuple.h>

void EastWestOffset()
{
  int run = 411768;
  int iter = 1;

  TString pisaFileIn = Form("geom/svxPISA-%d.par", run);
  if (iter > 0)
    pisaFileIn += Form(".%d", iter);

  TString fnames[] = {"",
                      "july3_parv1_small",
                      "july11_v3_500kevents"
                     };
  TFile *f = new TFile(Form("rootfiles/%d_%s.root", run, fnames[iter].Data()));
  TNtuple *t = (TNtuple *)f->Get("seg_clusntuple");

  SvxTGeo *geo = new SvxTGeo;
  geo->ReadParFile(pisaFileIn.Data());
  geo->MakeTopVolume(100, 100, 100);
  geo->AddSensors();
  geo->GeoManager()->CloseGeometry();

  geoEvents events;
  GetEventsFromTree(t, geo, events, (int)1e6);
  FitTracks(events);

  TH1D *hdv[3];
  for (int k=0; k<3; k++)
    hdv[k] = new TH1D(Form("hdv%d",k), Form("hdv%d",k), 100, -1, 1);

  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > 10)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");
      // Printf("%d tracks, E: (%.3f,%.3f,%.3f) W:(%.3f,%.3f,%.3f)",
      //        ntrk,ve(0),ve(1),ve(2),vw(0),vw(1),vw(2));
      for (int k=0; k<3; k++)
        hdv[k]->Fill(vw(k) - ve(k));
    }
  }
  
  for (int k=0; k<3; k++)
    DrawObject(hdv[k]);

  return;
}


