#include "UtilFns.h"
#include "VtxIO.h"
#include "VertexFinder.h"
#include "DcaFunctions.h"

#include <TGeoManager.h>
#include <TFile.h>
#include <TNtuple.h>

void EastWestOffset()
{
  int run = 411768;
  int iter = 2;
  gStyle->SetPalette(56, 0, 0.5); // Inverted "radiator", 50% transparency
  TObjArray *cList = new TObjArray();

  TString pisaFileIn = Form("geom/svxPISA-%d.par", run);
  if (iter > 0)
    pisaFileIn += Form(".%d", iter+1);

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
  GetEventsFromTree(t, geo, events);
  FitTracks(events);

  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k=0; k<3; k++)
    hdv[k] = new TH1D(Form("hdv%d",k),
                      Form("West - East vertex difference #Delta%s;#Delta%s [cm];events",
                           xyzstr[k], xyzstr[k]),
                      100, -1, 1);

  // DCA vs phi histogram and TProfile
  double dmin = 0.0, dmax = 0.1;
  TH2D *hbp = new TH2D("hbp", "", 100, 0, TMath::TwoPi(), 100, dmin, dmax);
  TProfile *prof = new TProfile("prof", "", 100, 0, TMath::TwoPi(), dmin, dmax);
  prof->SetMarkerStyle(kFullCircle);
  prof->BuildOptions(dmin, dmax, "g");

  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > 20)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");
      // Printf("%d tracks, E: (%.3f,%.3f,%.3f) W:(%.3f,%.3f,%.3f)",
      //        ntrk,ve(0),ve(1),ve(2),vw(0),vw(1),vw(2));

      for (int k=0; k<3; k++)
      {
        double dv = vw(k) - ve(k);
        if (TMath::Abs(dv) > 1e-15 && TMath::Abs(dv) < 1.0)
          hdv[k]->Fill(dv);
      }

      DcaVsPhi(events[ev], ve, vw, hbp, prof);
    }
  }

  for (int k=0; k<3; k++)
  {
    DrawObject(hdv[k], "", Form("cdv%d",k), cList);
  }

  DrawObject(hbp, "colz", "dcavsphi", cList, 900, 500);
  hbp->SetTitle("Zero-field DCA");
  TAxis *ax = hbp->GetXaxis();
  TAxis *ay = hbp->GetYaxis();
  ax->SetTitle("#Leftarrow  west      #phi + #pi/2      east  #Rightarrow");
  ay->SetTitle("DCA [cm]");
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(408);
  ay->SetNdivisions(205);
  gPad->SetGridy();
  prof->Draw("same");

  PrintPDFs(cList, Form("pdfs/iter%d", iter), "");
  PrintPDF(cList, Form("pdfs/dv-iter%d", iter));

  return;
}


