#include "TrackGen.h"
#include "VtxIO.h"

void GenerateEvents()
{
  int run = 123456;
  int prod = 0;
  const int nevents = (int)1e4;
  SvxTGeo *tgeo = VTXModel("geom/svxPISA-ideal.par");
  TString rootFileOut = Form("rootfiles/%d-%d-%d.root", run, prod, 0);
  TString pisaFileOut = Form("geom/%d-%d-%d.par", run, prod, 0);

  Printf("Generating %d events...", nevents);
  geoEvents events(nevents);
  TRandom3 ran;
  TVectorD vertex(3);
  for (int n=0; n<nevents; n++)
  {
    vertex(0) = 0.010*ran.Gaus();
    vertex(1) = 0.010*ran.Gaus();
    vertex(2) = 2*ran.Gaus();

    geoTracks tracks;
    int ntracks = 1 + ran.Integer(50);
    GenerateEvent(tgeo, ntracks, vertex, tracks);
    events[n] = tracks;
  }

  Printf("%lu events generated.", events.size());

  // Misalign half-layer 1W.
  // Then update global hit positions to reflect misalignments.
  // Local hit positions (x,z on sensor) remain unchanged.
  tgeo->TranslateHalfLayer(1, 1, 0.01, 0.01, 0.);
  for (unsigned int ev=0; ev<events.size(); ev++)
    for (unsigned int t=0; t<events[ev].size(); t++)
      events[ev][t].UpdateHits();

  // Here, fit tracks to (re)assign residuals, angles, and intercepts.
  // This puts the tracks on the same footing as the post-aligned refit sample.
  // Without this step, a comparison of residuals before and after alignment
  // would not be apples-to-apples.
  FitTracks(events);

  // Write out (misaligned) geometry to par file
  Printf("Writing %s", pisaFileOut.Data());
  tgeo->WriteParFile(pisaFileOut.Data());

  cout << "Filling output tree(s)..." << flush;
  TFile *outFile = new TFile(rootFileOut.Data(), "recreate");
  TNtuple *ht = new TNtuple("vtxhits", "VTX hit variables", HITVARS);
  FillNTuple(events, ht);
  Printf("done.");
  Printf("Writing %s", rootFileOut.Data());
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);
  Printf("Done!");

  return;
}