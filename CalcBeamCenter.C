#include "VtxAlignBase.h"
#include "GLSFitter.h"
#include "VertexFinder.h"
#include "DcaFunctions.h"
#include "BeamCenterFuncs.h"
#include "ParameterDefs.h"
#include "VtxIO.h"
#include "VtxVis.h"

void WriteConfigFile(const char *filename,
                     int run,
                     int prod,
                     int subiter,
                     TVectorD &bc,
                     TVectorD &ewdiff,
                     TVectorD &cntdiff,
                     TString notes = "");
TH1D *
ProfileSummary(TProfile *prof, const char *name);

void CalcBeamCenter(int run = 123456,
                    int prod = 0,
                    int subiter = 0)
{
  bool write = true;

  // Inputs:
  TString rootFileIn = Form("rootfiles/%d-%d-%d.root", run, prod, subiter);
  TString pisaFileIn = Form("geom/%d-%d-%d.par", run, prod, subiter);

  // Outputs:
  TString bcFileOut  = Form("rootfiles/bc-%d-%d-%d.root", run, prod, subiter);
  TString configFile = Form("production/config/config-%d-%d-%d.txt",
                            run, prod, subiter);

  TFile *f = new TFile(rootFileIn.Data());
  TTree *t = (TTree *)f->Get("vtxtrks");

  gStyle->SetOptStat(0);
  gStyle->SetPalette(56, 0, 0.5); // Inverted "radiator", 50% transparency
  TObjArray *cList = new TObjArray();
  TLatex ltx;
  ltx.SetNDC();
  ltx.SetTextFont(42);

  SvxTGeo *geo = VTXModel(pisaFileIn.Data());

  geoEvents events;
  GetEventsFromTree(t, geo, events, -1);
  FitTracks(events);

  TH1D *hdcae = new TH1D("hdcae", "", 200, -0.05, 0.05);
  TH1D *hdcaw = new TH1D("hdcaw", "", 200, -0.05, 0.05);
  TH2D *hdca2d   = new TH2D("hdca2d", ";#phi [rad];DCA [cm]",
                            100, -TMath::PiOver2(), 3*TMath::PiOver2(),
                            100, -0.05, +0.05);

  // Event multiplicity histograms
  TH1D *hne = new TH1D("hne", Form("VTXE multiplicity - prod %d step %d;"
                                   "tracks/event;events",
                                   prod, subiter), 100, 0, 100);
  TH1D *hnw = new TH1D("hnw", Form("VTXW multiplicity - prod %d step %d;"
                                   "tracks/event;events",
                                   prod, subiter), 100, 0, 100);

  // x-y vertex distributions
  double x0 = -1.0, y0 = -1.0, x1 = +1.0, y1 = +1.0;
  TH2D *hve = new TH2D("hve", Form("East vertex - prod %d step %d;"
                                   "x [cm];y [cm]", prod, subiter),
                       500, x0, x1, 500, y0, y1);
  TH2D *hvw = new TH2D("hvw", Form("West vertex - prod %d step %d;"
                                   "x [cm];y [cm]", prod, subiter),
                       500, x0, x1, 500, y0, y1);

  // West - East offset histograms
  TH1D *hdv[3];
  const char *xyzstr[3] = {"x", "y", "z"};
  for (int k=0; k<3; k++)
    hdv[k] = new TH1D(Form("hdv%d",k),
                      Form("West - East vertex difference #Delta%s "
                           "- prod %d step %d; #Delta%s [cm];events",
                           xyzstr[k], prod, subiter, xyzstr[k]),
                      200, -1, 1);

  Printf("Filling mult/vertex/DCA distributions...");
  int minmult = 10;
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ne=0, nw=0;
    int mult = events[ev].size();
    for (int t=0; t<mult; t++)
    {
      SvxGeoTrack trk = events[ev][t];
      if (East(trk.phi0))
      {
        ne++;
        hdcae->Fill(trk.xydca);
      }
      else
      {
        nw++;
        hdcaw->Fill(trk.xydca);
      }

      double p = trk.phi0;
      if (p > 3*TMath::PiOver2())
        p -= TMath::TwoPi();
      hdca2d->Fill(p, trk.xydca);

    }

    hne->Fill(ne);
    hnw->Fill(nw);

    if (mult >= minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east", "");
      TVectorD vw = Vertex(events.at(ev), "west", "");

      hve->Fill(ve(0), ve(1));
      hvw->Fill(vw(0), vw(1));

      for (int k=0; k<3; k++)
        hdv[k]->Fill(vw(k) - ve(k));
    }
  }

  Printf("Computing quantiles for {e,w}{x,y}");
  const int nq = 5;
  double probs[nq] = {0.01, 0.32, 0.50, 0.68, 0.99};
  double qxe[nq] = {0}, qxw[nq] = {0}, qye[nq] = {0}, qyw[nq] = {0};
  hve->ProjectionX()->GetQuantiles(nq, qxe, probs);
  hvw->ProjectionX()->GetQuantiles(nq, qxw, probs);
  hve->ProjectionY()->GetQuantiles(nq, qye, probs);
  hvw->ProjectionY()->GetQuantiles(nq, qyw, probs);

  TVectorD bce(2); bce(0) = qxe[2]; bce(1) = qye[2];
  TVectorD bcw(2); bcw(0) = qxw[2]; bcw(1) = qyw[2];
  TGraphErrors *gbc = new TGraphErrors();
  gbc->SetTitle("Beamcenter from east arm (point 0) and west arm (point 1)");
  gbc->SetPoint(0, bce(0), bce(1));
  gbc->SetPoint(1, bcw(0), bcw(1));
  gbc->SetPointError(0, 0.5*(qxe[3]-qxe[1]), 0.5*(qye[3]-qye[1]));
  gbc->SetPointError(1, 0.5*(qxw[3]-qxw[1]), 0.5*(qyw[3]-qyw[1]));

  // Zoom in on beam spot, keeping 98% of the data in view.
  // Maintain a 1:1 aspect ratio.
  double eastrange = TMath::Max(qxe[nq-1]-qxe[0], qye[nq-1]-qye[0]);
  double westrange = TMath::Max(qxw[nq-1]-qxw[0], qyw[nq-1]-qyw[0]);
  if (eastrange > 1.)
    eastrange = TMath::Min(westrange, 1.);
  if (westrange > 1.)
    westrange = TMath::Min(eastrange, 1.);
  x0 = qxe[2] - eastrange/2;
  y0 = qye[2] - eastrange/2;
  x1 = qxe[2] + eastrange/2;
  y1 = qye[2] + eastrange/2;
  hve->GetXaxis()->SetRangeUser(x0, x1);
  hve->GetYaxis()->SetRangeUser(y0, y1);
  x0 = qxw[2] - westrange/2;
  y0 = qyw[2] - westrange/2;
  x1 = qxw[2] + westrange/2;
  y1 = qyw[2] + westrange/2;
  hvw->GetXaxis()->SetRangeUser(x0, x1);
  hvw->GetYaxis()->SetRangeUser(y0, y1);

  Printf("Refitting tracks using beam center...");
  FitTracks(events, gbc);

  Printf("DCA x-y distribution...");
  x0 = -0.12;
  y0 = -0.12;
  x1 = +0.12;
  y1 = +0.12;
  TH2D *hxye = new TH2D("hxye", "DCA (east);x [cm];y [cm]",
                        60, x0, x1, 60, y0, y1);
  TH2D *hxyw = new TH2D("hxyw", "DCA (west);x [cm];y [cm]",
                        60, x0, x1, 60, y0, y1);
  TH1D *hre = new TH1D("hre", "hre", 100, 0, 0.5);
  TH1D *hrw = new TH1D("hrw", "hrw", 100, 0, 0.5);
  hxye->SetLineColorAlpha(kRed+2, 0.4);
  hxyw->SetLineColorAlpha(kRed+2, 0.4);
  DcaDist(events, "east", hxye, hre); // DCA to primary vertex
  DcaDist(events, "west", hxyw, hrw);

  Printf("DCA to beam center vs phi");
  int nphibins = 200;
  double dmin = 0.0, dmax = 0.05; // 500 um
  TH2D *hbp = new TH2D("hbp", "", nphibins, 0, TMath::TwoPi(), 100, dmin, dmax);
  TProfile *bprof = new TProfile("bprof", "", nphibins, 0, TMath::TwoPi(), dmin, dmax);
  bprof->SetMarkerStyle(kFullCircle);

  // Error calculation options:
  // "" : error = sigma/sqrt(N) = std. dev of mean.
  // "s": error = sigma = std. dev of distribution in dmin, dmax.
  // "g": error = 1/sqrt(sum(w)) where w is (weighted) sum of entries.
  bprof->BuildOptions(dmin, dmax, "g");
  DcaVsPhi(events, bce, bcw, hbp, bprof);

  Printf("DCA to primary vertex (in each arm) vs phi...");
  TH2D *hdp = new TH2D("hdp", "", nphibins, 0, TMath::TwoPi(), 100, dmin, dmax);
  TProfile *dprof = new TProfile("dprof", "", nphibins, 0, TMath::TwoPi(), dmin, dmax);
  dprof->SetMarkerStyle(kFullCircle);
  dprof->BuildOptions(dmin, dmax, "g");

  // This is inefficient, repeating all this work...improve later
  for (unsigned int ev=0; ev<events.size(); ev++)
  {
    int ntrk = events[ev].size();
    if (ntrk > minmult)
    {
      TVectorD ve = Vertex(events.at(ev), "east");
      TVectorD vw = Vertex(events.at(ev), "west");
      DcaVsPhi(events[ev], ve, vw, hdp, dprof);
    }
  }

  // TVectorD avgbc = 0.5*(bce + bcw);
  TVectorD ewdiff(3);
  for (int k=0; k<3; k++)
    ewdiff(k) = hdv[k]->GetMean();
  TVectorD cntdiff(3);

  if (write)
  {
    WriteConfigFile(configFile.Data(),
                    run,
                    prod,
                    subiter,
                    bcw,
                    ewdiff,
                    cntdiff,
                    "Using west beam center.");

    TFile *bcf = new TFile(bcFileOut.Data(), "recreate");
    gbc->Write("gbc");
    hve->Write();
    hvw->Write();
    bcf->Close();
  }


  // ---- ---- ---- ---- ---- ---- ----
  //                Plot
  // ---- ---- ---- ---- ---- ---- ----
  Printf("Plotting...");
  DrawObject(hne, "", "east_mult", cList);
  gPad->SetLogx();
  gPad->SetLogy();
  DrawObject(hnw, "", "west_mult", cList);
  gPad->SetLogx();
  gPad->SetLogy();


  TCanvas *cdca = new TCanvas("cdca", "dca_east_west", 1000, 500);
  SetYMax(hdcae, hdcaw);
  cdca->Divide(2,1,0,0);
  cdca->cd(1);
  hdcae->Draw();
  ltx.DrawLatex(0.15, 0.90, Form("East: (Mean, Std Dev) = (%.0f, %.0f) #mum",
                                 1e4*hdcae->GetMean(), 1e4*hdcae->GetRMS()));
  cdca->cd(2);
  hdcaw->Draw();
  ltx.DrawLatex(0.15, 0.90, Form("West: (Mean, Std Dev) = (%.0f, %.0f) #mum",
                                 1e4*hdcaw->GetMean(), 1e4*hdcaw->GetRMS()));
  cList->Add(cdca);

  DrawObject(hdca2d, "col", "dca2d", cList);
  hdca2d->GetXaxis()->CenterTitle();
  hdca2d->GetYaxis()->CenterTitle();
  TProfile* dca2dprof = hdca2d->ProfileX("dca2dprof", 1, -1, "d,same");
  dca2dprof->SetMarkerStyle(kFullCircle);
  ltx.DrawLatex(0.25, 0.80, Form("West"));
  ltx.DrawLatex(0.65, 0.80, Form("East"));

  TCanvas *cve = new TCanvas("cve", "cve", 500, 500);
  hve->Draw("col");
  ltx.DrawLatex(0.15, 0.85, Form("BC = %.3f, %.3f", bce(0), bce(1)));
  ltx.DrawLatex(0.15, 0.80, Form("#sigma_{qx} = %.3f, #sigma_{qy} = %.3f",
                                 0.5*(qxe[3]-qxe[1]),
                                 0.5*(qye[3]-qye[1])));
  cList->Add(cve);
  TCanvas *cvw = new TCanvas("cvw", "cvw", 500, 500);
  hvw->Draw("col");
  ltx.DrawLatex(0.15, 0.85, Form("BC = %.3f, %.3f", bcw(0), bcw(1)));
  ltx.DrawLatex(0.15, 0.80, Form("#sigma_{qx} = %.3f, #sigma_{qy} = %.3f",
                                 0.5*(qxw[3]-qxw[1]),
                                 0.5*(qyw[3]-qyw[1])));
  cList->Add(cvw);

  for (int k=0; k<3; k++)
  {
    DrawObject(hdv[k], "", Form("cdv%d",k), cList);
    ltx.DrawLatex(0.15, 0.85, Form("Mean %.3f", hdv[k]->GetMean()));
    ltx.DrawLatex(0.15, 0.80, Form("Std dev %.3f", hdv[k]->GetRMS()));
  }

  // DCA contour plot
  TCanvas *c = new TCanvas("c", "dca2d", 1000, 500);
  c->Divide(2,1,0,0);
  c->cd(1);
  TAxis *ax = hxye->GetXaxis();
  TAxis *ay = hxye->GetYaxis();
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(208);
  ay->SetNdivisions(208);
  hxye->Draw("cont3");
  gPad->SetLogz();
  c->cd(2);
  ax = hxyw->GetXaxis();
  ay = hxyw->GetYaxis();
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(208);
  ay->SetNdivisions(208);
  hxyw->Draw("cont3");
  gPad->SetLogz();
  cList->Add(c);

  TCanvas *cr = new TCanvas("cr", "dcar", 500, 500);
  SetHistProps(hre,kAzure+1,kNone,kAzure+1);
  SetHistProps(hrw,kOrange,kNone,kOrange);
  hre->SetLineWidth(2);
  hrw->SetLineWidth(2);
  hrw->SetTitle("Zero-field DCA #color[861]{ east}, #color[800]{ west};DCA [cm];");
  hrw->Draw("");
  hre->Draw("same");
  ltx.DrawLatex(0.2, 0.85, Form("East mean %.3f, std dev %.3f",
                                hre->GetMean(), hre->GetRMS()));
  ltx.DrawLatex(0.2, 0.80, Form("West mean %.3f, std dev %.3f",
                                hrw->GetMean(), hrw->GetRMS()));
  cr->SetLogy();
  cList->Add(cr);

  TH1D *bpsummary = ProfileSummary(bprof, "BC DCA vs phi summary");
  TH1D *dpsummary = ProfileSummary(dprof, "DCA vs phi summary");
  string bps = Form("mean %.4f, std dev %.4f",
                    bpsummary->GetMean(), bpsummary->GetRMS());
  string dps = Form("mean %.4f, std dev %.4f",
                    dpsummary->GetMean(), dpsummary->GetRMS());

  DrawObject(hbp, "colz", "bcdcavsphi", cList, 900, 500);
  hbp->SetTitle("DCA to beam center");
  ax = hbp->GetXaxis();
  ay = hbp->GetYaxis();
  ax->SetTitle("#Leftarrow  west      #phi + #pi/2      east  #Rightarrow");
  ay->SetTitle("DCA [cm]");
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(408);
  ay->SetNdivisions(205);
  gPad->SetGridy();
  bprof->Draw("same");
  ltx.DrawLatex(0.2, 0.85, bps.c_str());

  DrawObject(hdp, "colz", "dcavsphi", cList, 900, 500);
  hdp->SetTitle("Zero-field DCA");
  ax = hdp->GetXaxis();
  ay = hdp->GetYaxis();
  ax->SetTitle("#Leftarrow  west      #phi + #pi/2      east  #Rightarrow");
  ay->SetTitle("DCA [cm]");
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetNdivisions(408);
  ay->SetNdivisions(205);
  gPad->SetGridy();
  dprof->Draw("same");
  ltx.DrawLatex(0.2, 0.85, dps.c_str());

  DrawObject(bpsummary, "", "bp_summary", cList, 500, 500);
  ltx.DrawLatex(0.2, 0.85, bps.c_str());

  DrawObject(dpsummary, "", "dp_summary", cList, 500, 500);
  ltx.DrawLatex(0.2, 0.85, dps.c_str());

  Printf("Fraction outside 1000um: %.3f (e) %.3f (w)",
         1.0 - hre->Integral(1,hre->FindBin(0.099))/hre->Integral(),
         1.0 - hrw->Integral(1,hrw->FindBin(0.099))/hrw->Integral());

  if (write)
  {
    PrintPDFs(cList, Form("pdfs/run%d-prod%d-subiter%d", run, prod, subiter), "");
    PrintPDF(cList, Form("pdfs/beam-center-run%d-pro%d-sub%d", run, prod, subiter));
  }

  return;
}

void
WriteConfigFile(const char *filename,
                int run,
                int prod,
                int subiter,
                TVectorD &bc,
                TVectorD &ewdiff,
                TVectorD &cntdiff,
                TString notes)
{
  cout << "Writing " << filename << "..." << flush;

  ofstream fs(filename);

  // Metadata - parsed by humans - start lines with #
  fs << "# This is " << gSystem->GetFromPipe("pwd").Data()
     << filename << endl;
  fs << "# Created by user " << gSystem->GetFromPipe("whoami").Data() << endl;
  fs << "# On " << gSystem->GetFromPipe("date").Data() << endl;
  fs << "# At " << gSystem->GetFromPipe("hostname").Data() << endl;
  fs << "# Run number: " << run << endl;
  fs << "# Production step: " << prod << endl;
  fs << "# Sub-iteration step: " << subiter << endl;
  fs << "# Source: "
     << gSystem->GetFromPipe("git config --get remote.origin.url").Data()
     << endl;
  fs << "# Latest commit: "
     << gSystem->GetFromPipe("git rev-parse HEAD").Data()
     << endl;
  fs << "# Notes: " << notes.Data() << endl;
  fs << endl;

  // Production configuration data below
  fs << "beamcenter: " << bc(0) << " " << bc(1) << endl;
  fs << "east-to-west: " << ewdiff(0) << " " << ewdiff(1) << " " << ewdiff(2)
     << endl;
  fs << "vtx-to-cnt: " << cntdiff(0) << " " << cntdiff(1) << " " << cntdiff(2)
     << endl;
  fs << "geomfile: " << Form("%d-%d-%d.par", run, prod, subiter) << endl;

  fs.close();

  Printf("done.");

  return;
}

TH1D *
ProfileSummary(TProfile *prof, const char *name)
{
  TH1D *h = new TH1D(name, name, 100, prof->GetMinimum(), prof->GetMaximum());
  for (int i=1; i<=prof->GetNbinsX(); ++i)
  {
    double e = prof->GetBinError(i);
    h->Fill(prof->GetBinContent(i), e>0 ? 1./e : 0.);
  }

  return h;
}
