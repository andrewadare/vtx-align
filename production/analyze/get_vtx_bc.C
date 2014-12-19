////////////////////////////////////////////////////////////////////////////////
//
// Get the VTX Beamcenter for a given run
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 12-18-2014
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>

#include <iostream>

using namespace std;

void get_vtx_bc()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //===================================================//
  // SET RUNNING CONDITIONS
  //===================================================//

  const char *fileName =
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-405869-22-11/testvtxproduction_fieldon-405869-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-22-11_1/testvtxproduction_fieldon-407951-22-11_1.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-408219-22-11/testvtxproduction_fieldon-408219-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-408408-22-11/testvtxproduction_fieldon-408408-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-409702-22-11/testvtxproduction_fieldon-409702-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-410112-22-11/testvtxproduction_fieldon-410112-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-410682-22-11/testvtxproduction_fieldon-410682-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-410942-22-11/testvtxproduction_fieldon-410942-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-411644-22-11/testvtxproduction_fieldon-411644-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-411914-22-11/testvtxproduction_fieldon-411914-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-413721-22-11/testvtxproduction_fieldon-413721-22-11.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-414431-22-11/testvtxproduction_fieldon-414431-22-11.root";
    "/direct/phenix+prod01/phnxreco/millepede/ZF_production/zf-405836-22-11/testvtxproduction_zf-405836-22-11.root";

  bool print_plots = true;
  const char *fileLabel = "zf-405836-22-11";

  const char *eventCuts = "!tickCut && "
                          "TMath::Abs(vtx[2])<10 && "
                          // "(bbcq[0]+bbcq[1])>400 && " //central
                          // "(bbcq[0]+bbcq[1])<181 && " //peripheral
                          "which_vtx==\"SVX_PRECISE\"";

  //===================================================//
  // DELCARE VARIABLES
  //===================================================//

  TTree *ntp_event;

  TH1F *hvtx[2];

  TF1 *fgaus = new TF1("fgaus", "gaus", -0.4, 0.4);
  fgaus->SetLineColor(kRed);

  //===================================================//
  // GET TREE & HISTOGRAMS FROM FILE
  //===================================================//
  cout << endl;
  cout << "--> Getting data from " << fileName << endl;

  TFile *fin = TFile::Open(fileName);
  if (!fin)
  {
    cout << "ERROR!!! Unable to open " << fileName << endl;
    return;
  }

  ntp_event = (TTree *) fin->Get("ntp_event");
  if (!ntp_event)
  {
    cout << "ERROR!!! Unable to get ntp_event from " << fileName << endl;
    return;
  }

  ntp_event->Draw("vtx[0]>>htmp(1000,-0.5,0.5)", eventCuts, "goff");
  hvtx[0] = (TH1F *) gDirectory->FindObject("htmp");
  hvtx[0]->SetDirectory(0);
  hvtx[0]->SetName("vtx_0");
  hvtx[0]->SetTitle(";X_{vtx}");


  ntp_event->Draw("vtx[1]>>htmp(1000,-0.5,0.5)", eventCuts, "goff");
  hvtx[1] = (TH1F *) gDirectory->FindObject("htmp");
  hvtx[1]->SetDirectory(0);
  hvtx[1]->SetName("vtx_1");
  hvtx[1]->SetTitle(";Y_{vtx}");

  //===================================================//
  // FIT
  //===================================================//
  cout << endl;
  cout << "--> Fitting bc" << endl;

  for (int i = 0; i < 2; i++)
  {
    double mean = hvtx[i]->GetMean();
    double rms = hvtx[i]->GetRMS();

    fgaus->SetRange(mean - rms, mean + rms);
    fgaus->SetParameters(hvtx[i]->GetMaximum(), mean, rms);
    hvtx[i]->Fit(fgaus, "RQ0N");

    //we can get some weird shapes, so refit using the output of the first fit
    mean = fgaus->GetParameter(1);
    rms = fgaus->GetParameter(2);
    fgaus->SetRange(mean - 0.75 * rms, mean + 0.75 * rms);
    fgaus->SetParameters(hvtx[i]->GetMaximum(), mean, rms);
    hvtx[i]->Fit(fgaus, "RQ0");

    cout << " vtx[" << i << "]:"
         << " mean=" << fgaus->GetParameter(1) << " +/- " << fgaus->GetParError(1)
         << " sig=" << fgaus->GetParameter(2) << " +/- " << fgaus->GetParError(2)
         << endl;
  }

  //===================================================//
  // PLOT OBJECT
  //===================================================//
  cout << endl;
  cout << "--> Plotting" << endl;

  TLatex ltitle;
  ltitle.SetNDC();
  ltitle.SetTextAlign(22);

  TLatex lstat;
  lstat.SetNDC();

  //===================================================//
  // PLOT
  //===================================================//

  TCanvas *cvtx = new TCanvas("cvtx", "vtx", 1200, 700);
  cvtx->Divide(2, 1);

  for (int i = 0; i < 2; ++i)
  {
    cvtx->cd(i + 1);
    hvtx[i]->Draw();

    TF1 *ftmp = hvtx[i]->GetFunction("fgaus");
    ftmp->DrawCopy("same");

    ltitle.DrawLatex(0.5, 0.95, fileLabel);
    lstat.DrawLatex(0.6, 0.8, Form("mean=%.4f", ftmp->GetParameter(1)));
    lstat.DrawLatex(0.6, 0.7, Form("#sigma=%.4f", ftmp->GetParameter(2)));
  }

  //===================================================//
  // PRINT PLOTS
  //===================================================//
  if (print_plots)
  {
    cout << endl;
    cout << "--> Printing plots" << endl;

    cvtx->Print(Form("pdfs/bc/bc_%s.pdf", fileLabel));
  }

}