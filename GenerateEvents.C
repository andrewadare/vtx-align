#include "TrackGen.h"
#include "VtxIO.h"

void GenerateEvents()
{
  int run = 123456;
  int prod = 0;
  const int nvtxevents = (int)9e3;
  const int ncntevents = (int)1e3;
  SvxTGeo *tgeo = VTXModel("geom/svxPISA-ideal.par");
  TString rootFileOut = Form("rootfiles/%d-%d-%d.root", run, prod, 0);
  TString pisaFileOut = Form("geom/%d-%d-%d.par", run, prod, 0);

  geoEvents vtxevents(nvtxevents);
  geoEvents cntevents(ncntevents);
  TRandom3 ran;
  TVectorD vertex(3);

  // Assign a fixed beam x,y position if available
  double bcx = 0, bcy = 0;  // 0.3532, 0.0528
  TGraphErrors *gbc = 0;
  TFile *bcf = new TFile(Form("rootfiles/bc-%d.root", run));
  if (bcf)
    gbc = (TGraphErrors *)bcf->Get("gbc");
  if (gbc)
  {
    bcx = gbc->GetX()[1];
    bcy = gbc->GetY()[1];
    Info("", "Using a fixed beamcenter from %s", bcf->GetName());
    Info("", " West: (%.4f, %.4f)", bcx, bcy);
  }

  Printf("Generating %d VTX standalone events...", nvtxevents);
  for (int n=0; n<nvtxevents; n++)
  {
    vertex(0) = bcx + 0.010*ran.Gaus();
    vertex(1) = bcy + 0.010*ran.Gaus();
    vertex(2) = 2*ran.Gaus();

    geoTracks tracks;
    int ntracks = 10 + ran.Integer(40);
    GenerateEvent(tgeo, ntracks, vertex, tracks);

    // Save vertex
    for (unsigned int t=0; t<tracks.size(); t++)
    {
      tracks[t].vx = vertex(0);
      tracks[t].vy = vertex(1);
      tracks[t].vz = vertex(2);
      if (n<10 && t==0) Printf("bc %f %f\n %f %f %f", 
                       bcx,bcy,
                       tracks[t].vx,
                       tracks[t].vy,
                       tracks[t].vz);
    }

    vtxevents[n] = tracks;
  }
  Printf("%lu VTX events generated.", vtxevents.size());

  // Generate fake SvxCentralTracks
  Printf("Generating %d CNT events...", ncntevents);
  for (int n=0; n<ncntevents; n++)
  {
    vertex(0) = bcx + 0.010*ran.Gaus();
    vertex(1) = bcy + 0.010*ran.Gaus();
    vertex(2) = 2*ran.Gaus();

    geoTracks tracks;
    int ntracks = 10 + ran.Integer(40);
    GenerateEvent(tgeo, ntracks, vertex, tracks);

    // Save vertex
    for (unsigned int t=0; t<tracks.size(); t++)
    {
      tracks[t].vx = vertex(0);
      tracks[t].vy = vertex(1);
      tracks[t].vz = vertex(2);
    }

    cntevents[n] = tracks;
  }
  Printf("%lu CNT events generated.", cntevents.size());

  // ** Apply misalignments here **
  int layer = 1;
  int arm = 0; // East = 0, West = 1

  // Experiment 1:
  // tgeo->RotateHalfLayerRPhi(layer, arm, 0.04);

  // Experiment 2:
  tgeo->RotateHalfLayerRPhi(0, arm, 0.02);
  tgeo->RotateHalfLayerRPhi(1, arm, -0.02);
  tgeo->RotateHalfLayerRPhi(2, arm, 0.03);
  tgeo->RotateHalfLayerRPhi(3, arm, -0.03);
  tgeo->TranslateArm(arm, -0.05, 0.02, 0.0);
  tgeo->TranslateHalfLayer(0, arm, 0.02, 0.02, -0.01);
  tgeo->TranslateHalfLayer(1, arm, -0.01, 0.01, 0.03);
  tgeo->TranslateHalfLayer(2, arm, 0.01, 0.005, -0.005);
  tgeo->TranslateHalfLayer(3, arm, 0.00, -0.02, 0.0);

  // Apply random Gaussian azimuthal misalignments
  for (int a=0; a<2; a++)
    for (int lyr=0; lyr<4; lyr++)
    {
      int first = -1, last  = -1;
      tgeo->LadderRange(lyr, a, first, last);
      for (int ldr = first; ldr <= last; ldr++)
        tgeo->RotateLadderRPhi(lyr, ldr, 0.002*ran.Gaus()); // 20 microns
    }

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
  FitTracks(vtxevents, gbc, "calc_dca");

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
