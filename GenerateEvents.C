#include "TrackGen.h"
#include "VtxIO.h"

void GenerateEvents()
{
  int run = 123456;
  int prod = 0;
  const int nvtxevents = (int)9e1;
  const int ncntevents = (int)1e1;
  SvxTGeo *tgeo = VTXModel("geom/svxPISA-ideal.par");
  TString rootFileOut = Form("rootfiles/%d-%d-%d.root", run, prod, 0);
  TString pisaFileOut = Form("geom/%d-%d-%d.par", run, prod, 0);

  geoEvents vtxevents(nvtxevents);
  geoEvents cntevents(ncntevents);
  TRandom3 ran;
  TVectorD vertex(3);

  Printf("Generating %d VTX standalone events...", nvtxevents);
  for (int n=0; n<nvtxevents; n++)
  {
    vertex(0) = 0.010*ran.Gaus();
    vertex(1) = 0.010*ran.Gaus();
    vertex(2) = 2*ran.Gaus();

    geoTracks tracks;
    int ntracks = 1 + ran.Integer(50);
    GenerateEvent(tgeo, ntracks, vertex, tracks);
    vtxevents[n] = tracks;
  }
  Printf("%lu VTX events generated.", vtxevents.size());

  // Generate fake SvxCentralTracks
  Printf("Generating %d CNT events...", ncntevents);
  for (int n=0; n<ncntevents; n++)
  {
    vertex(0) = 0.010*ran.Gaus();
    vertex(1) = 0.010*ran.Gaus();
    vertex(2) = 2*ran.Gaus();

    geoTracks tracks;
    int ntracks = 1 + ran.Integer(50);
    GenerateEvent(tgeo, ntracks, vertex, tracks);
    cntevents[n] = tracks;
  }
  Printf("%lu CNT events generated.", cntevents.size());

  // ** Apply misalignments here **
  int layer = 1;
  int arm = 1; // East = 0, West = 1
  // tgeo->TranslateHalfLayer(layer, arm, 0.02, 0.02, 0.0);
  tgeo->RotateHalfLayerRPhi(layer, arm, 0.04);

  // Update global hit positions to reflect misalignments.
  // Local hit positions (x,z on sensor) remain unchanged.
  for (unsigned int ev=0; ev<vtxevents.size(); ev++)
    for (unsigned int t=0; t<vtxevents[ev].size(); t++)
      vtxevents[ev][t].UpdateHits();

  for (unsigned int ev=0; ev<cntevents.size(); ev++)
    for (unsigned int t=0; t<cntevents[ev].size(); t++)
    {
      for (int ihit=0; ihit<cntevents[ev][t].nhits; ihit++)
        cntevents[ev][t].hits[ihit].UpdateResiduals();

      cntevents[ev][t].UpdateHits();
    }

  // Fit VTX tracks to (re)assign residuals, angles, and intercepts.
  // This puts the tracks on the same footing as the post-aligned refit sample.
  // Without this step, a comparison of residuals before and after alignment
  // would not be apples-to-apples.
  // CNT tracks are not fit here--their parameters are determined externally.
  FitTracks(vtxevents);

  // Write out (misaligned) geometry to par file
  Printf("Writing %s", pisaFileOut.Data());
  tgeo->WriteParFile(pisaFileOut.Data());

  cout << "Filling output tree(s)..." << flush;

  TFile *outFile = new TFile(rootFileOut.Data(), "recreate");
  TTree *vtxtree = CreateTree("vtxtrks");
  TTree *cnttree = CreateTree("cnttrks");
  FillTree(vtxevents, vtxtree);
  FillTree(cntevents, cnttree);

  Printf("done.");
  Printf("Writing %s", rootFileOut.Data());
  outFile->cd();
  outFile->Write(0, TObject::kOverwrite);
  Printf("Done!");

  return;
}
