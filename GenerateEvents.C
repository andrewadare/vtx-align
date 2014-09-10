#include "TrackGen.h"
#include "VtxIO.h"

void GenerateEvents()
{
  const int nevents = (int)1e4;
  SvxTGeo *tgeo = VTXModel("geom/svxPISA-ideal.par");
  TString rootFileOut = Form("rootfiles/%d-%d-%d.root", 123456, 0, 0);

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

  // Here, fit tracks to (re)assign residuals, angles, and intercepts.
  // This puts the tracks on the same footing as the post-aligned refit sample.
  // Without this step, a comparison of residuals before and after alignment 
  // would not be apples-to-apples.
  FitTracks(events);

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