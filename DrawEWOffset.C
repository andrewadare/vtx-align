#include "VtxAlignBase.h"
#include "VtxVis.h"
#include "VtxIO.h"
#include "GLSFitter.h"


bool openPDF = true; // Mac only, but simple to adapt to other viewers with CLI

void SetAxisProps(TH1 *h, int xdivs = 208, int ydivs = 208,
                  double labelSize = 0.04, double titleSize = 0.04,
                  double xTitleOffset = 1.1, double yTitleOffset = 1.5,
                  bool centerX = true, bool centerY = true);


void DrawEWOffset(int run = 411768,
                  int prod = 1,
                  int subit = 0,
                  const char *treename = "vtxtrks")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  
  TObjArray *cList = new TObjArray();

  //load the par file
  std::string geoFile = Form("geom/%i-%i-%i.par", run, prod, subit);
  Info("", "--Reading geometry from %s--", geoFile.c_str());
  SvxTGeo *tgeo = VTXModel(geoFile.c_str());
  assert(tgeo);

  //input file
  std::string fName = Form("rootfiles/%d-%d-%d.root", run, prod, subit);
  Info("", "--Reading TTree from %s--", fName.c_str());
  TFile *inFile = new TFile(fName.c_str(), "read");
  assert(inFile);

  //input ttree
  TTree *t = (TTree *) inFile->Get(treename);
  assert(t);

  //read tree into geoEvents
  Info("", "--Reading events from tree %s--", treename);
  geoEvents events;
  GetEventsFromTree(t, tgeo, events);

  //refit w/o vertex and recalculate east/west vertex
  Info("", "--Re-fitting tracks--");
  FitTracks(events, 0, "find_vertex");

  
  //setup histograms
  Info("", "--Intializing histograms--");
  TH1D *hvze = new TH1D(Form("hvze_%d_%d", prod, subit),
                        Form(";east arm z vertex [cm];tracks"), 200, -15, 15);
  TH1D *hvzw = new TH1D(Form("hvzw_%d_%d", prod, subit),
                        Form(";west arm z vertex [cm];tracks"), 200, -15, 15);
  SetAxisProps(hvze);
  SetAxisProps(hvzw);


  // x-y vertex distributions
  double x0 = -0.5, y0 = -0.5, x1 = +0.5, y1 = +0.5;
  TH2D *hve = new TH2D(Form("hve_%d_%d", prod, subit),
                       Form("East vertex - prod %d step %d;"
                            "x [cm];y [cm]", prod, subit),
                       500, x0, x1, 500, y0, y1);
  SetAxisProps(hve, 205, 205);


  TH2D *hvw = new TH2D(Form("hvw_%d_%d", prod, subit),
                       Form("West vertex - prod %d step %d;"
                            "x [cm];y [cm]", prod, subit),
                       500, x0, x1, 500, y0, y1);
  SetAxisProps(hvw, 205, 205);

  // West - East offset histograms
  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k = 0; k < 3; k++)
  {
    double lim = 0.41;
    if (k == 2)
      lim *= 2;

    hdv[k] = new TH1D(Form("hdv_%d__%d_%d", prod, subit, k),
                      Form("West - East vertex difference #Delta%s "
                           "- prod %d step %d; (W-E) #Delta%s [cm];tracks",
                           xyzstr[k], prod, subit, xyzstr[k]),
                      200, -lim, lim);
    SetAxisProps(hdv[k]);
  }

  //Fill histograms
  Info("", "--Filling histograms--");

  for (unsigned int ev = 0; ev < events.size(); ev++)
  {
    TVectorD ve = RetrieveVertex(events[ev], "east");
    TVectorD vw = RetrieveVertex(events[ev], "west");

    //count the number of tracks in each arm
    int ne = 0;
    int nw = 0;
    for (unsigned int trk = 0; trk < events[ev].size(); trk++)
    {
      int arm = (events[ev][trk].hits[0].x < 0.) ? 0 : 1; // 0 = East, 1 = West.

      if (arm == 0) ne++;
      if (arm == 1) nw++;
    }

    if (ne > 1 && nw > 1)
    {
      for (int k = 0; k < 3; k++)
        hdv[k]->Fill(vw(k) - ve(k));
    }

    if (ne > 1)
    {
      hve->Fill(ve(0), ve(1));
      hvze->Fill(ve(2));
    }

    if (nw > 1)
    {
      hvw->Fill(vw(0), vw(1));
      hvzw->Fill(vw(2));
    }

  }

  // Plot--------------------------------------------------------
  Info("", "--Plotting histograms--");
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);
  TCanvas *c = 0;

  // East and West xy vertex distributions
  hve->GetXaxis()->SetRangeUser(hve->GetMean(1) - 0.1, hve->GetMean(1) + 0.1);
  hve->GetYaxis()->SetRangeUser(hve->GetMean(2) - 0.1, hve->GetMean(2) + 0.1);
  hvw->GetXaxis()->SetRangeUser(hvw->GetMean(1) - 0.1, hvw->GetMean(1) + 0.1);
  hvw->GetYaxis()->SetRangeUser(hvw->GetMean(2) - 0.1, hvw->GetMean(2) + 0.1);
  c = new TCanvas("xy_vertex",
                  "xy_vertex_ew", 1000, 600);
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.15); // L, R, B, T
  hve->Draw("col");
  ltx.DrawLatex(0.25, 0.95, "East");
  ltx.DrawLatex(0.45, 0.95, Form("mean (%.0f, %.0f) #mum",
                                 1e4 * hve->GetMean(1), 1e4 * hve->GetMean(2)));
  ltx.DrawLatex(0.45, 0.90, Form("std dev (%.0f, %.0f) #mum",
                                 1e4 * hve->GetRMS(1), 1e4 * hve->GetRMS(2)));
  c->cd(2);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.15); // L, R, B, T
  hvw->Draw("col");
  ltx.DrawLatex(0.25, 0.95, "West");
  ltx.DrawLatex(0.45, 0.95, Form("mean (%.0f, %.0f) #mum",
                                 1e4 * hvw->GetMean(1), 1e4 * hvw->GetMean(2)));
  ltx.DrawLatex(0.45, 0.90, Form("std dev (%.0f, %.0f) #mum",
                                 1e4 * hvw->GetRMS(1), 1e4 * hvw->GetRMS(2)));
  cList->Add(c);

  // East and West z vertex distributions
  c = new TCanvas("z_vertex",
                  "z_vertex_ew", 1000, 500);
  c->Divide(2, 1);
  c->cd(1);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  hvze->GetYaxis()->SetRangeUser(0, 1.2 * hvze->GetMaximum());
  hvze->Draw("");
  ltx.DrawLatex(0.25, 0.92, "East");
  ltx.DrawLatex(0.45, 0.92, Form("mean %.3f cm", hvze->GetMean(1)));
  ltx.DrawLatex(0.45, 0.87, Form("std dev %.3f cm", hvze->GetRMS(1)));
  c->cd(2);
  gPad->SetMargin(0.12, 0.02, 0.12, 0.02); // L, R, B, T
  hvzw->GetYaxis()->SetRangeUser(0, 1.2 * hvzw->GetMaximum());
  hvzw->Draw("");
  ltx.DrawLatex(0.25, 0.92, "West");
  ltx.DrawLatex(0.45, 0.92, Form("mean %.3f cm", hvzw->GetMean(1)));
  ltx.DrawLatex(0.45, 0.87, Form("std dev %.3f cm", hvzw->GetRMS(1)));
  cList->Add(c);

  // East-west offsets in x, y, and z
  c = new TCanvas("ew_offsets_%d",
                  "ew_offsets_%d", 1000, 500);
  c->Divide(3, 1, 0.001, 0.001);
  for (int k = 0; k < 3; k++)
  {
    c->cd(k + 1);
    gPad->SetMargin(0.15, 0.02, 0.12, 0.02); // L, R, B, T
    hdv[k]->GetYaxis()->SetRangeUser(0, 1.2 * hdv[k]->GetMaximum());
    hdv[k]->Draw();
    ltx.DrawLatex(0.25, 0.95, Form("Mean %.3f", hdv[k]->GetMean()));
    ltx.DrawLatex(0.25, 0.90, Form("Std dev %.3f", hdv[k]->GetRMS()));
  }
  cList->Add(c);


  Info("","--Printing PDFS--");
  PrintPDFs(cList, Form("pdfs/run%d-pro%dsub%d-%s",
                        run, prod, subit, treename), "");
  const char *pdfName = Form("pdfs/ewoffset-run%d-pro%dsub%d-%s",
                             run, prod, subit, treename);
  PrintPDF(cList, pdfName, "");
  if (openPDF)
    gSystem->Exec(Form("open %s.pdf", pdfName));


  return;
}

void
SetAxisProps(TH1 *h, int xdivs, int ydivs,
             double labelSize, double titleSize,
             double xTitleOffset, double yTitleOffset,
             bool centerX, bool centerY)
{
  TAxis *ax = h->GetXaxis();
  TAxis *ay = h->GetYaxis();
  ax->SetNdivisions(xdivs);
  ay->SetNdivisions(ydivs);

  ax->SetLabelSize(labelSize);
  ay->SetLabelSize(labelSize);

  ax->SetTitleSize(titleSize);
  ay->SetTitleSize(titleSize);

  ax->SetTitleOffset(xTitleOffset);
  ay->SetTitleOffset(yTitleOffset);

  ax->CenterTitle(centerX);
  ay->CenterTitle(centerY);

  return;
}

