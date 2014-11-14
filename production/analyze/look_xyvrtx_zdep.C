////////////////////////////////////////////////////////////////////////////////
//
// Look at the x & y, E & W, primary vertex as a function of the z vertex
// using output of TestVTXProduction module
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 11-14-2014
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>

#include <iostream>

using namespace std;

void look_xyvrtx_zdep()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //=============================================================//
  // SET RUNNING CONDITIONS
  //=============================================================//

  const char *inFile =
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-1-1/testvtxproduction_fieldon-407951-1-1.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-taebong-p2-v8_newfixes/testvtxproduction_407951-taebong-p2-v8-newfixes.root";
    "/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-0-3/testvtxproduction_407951-0-3.root";
  const bool print_plots = true;

  // beamcenter from config-fieldon-407951-1-1.txt
  // const float bc[2] = { -0.120756, 0.0741091};
  // beamcenter from config-taebong-p2-v8.txt
  // const float bc[2] = { -0.120756, 0.0741091};
  // beamcenter from config-zf-411768-0-3.txt
  const float bc[2] = {0.3532, 0.0528};

  //=============================================================//
  // DECLARE VARIABLES
  //=============================================================//

  TTree *ntp_event;

  TH2D *hvrtx_zvrtx[2][2]; //[E/W][X/Y]
  TProfile *pvrtx_zvrtx[2][2]; //[E/W][X/Y]

  char draw[500];


  TF1 *fline = new TF1("fline", "pol1", -8, 8);
  fline->SetLineColor(kOrange + 2);
  float slope[2][2] = {{0}, {0}};
  float slope_e[2][2] = {{0}, {0}};

  //=============================================================//
  // GET TREE FROM FILE
  //=============================================================//
  cout << endl;
  cout << "--> Getting TTree from " << inFile << endl;

  TFile *fin = TFile::Open(inFile);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << inFile << endl;
    return;
  }

  ntp_event = (TTree *) fin->Get("ntp_event");
  if (!ntp_event)
  {
    cout << "ERROR!! Unable to get ntp_event from " << inFile << endl;
    return;
  }

  //=============================================================//
  // MAKE HISTOGRAMS
  //=============================================================//
  cout << endl;
  cout << "--> Extract histograms" << endl;

  const char *armLabel[2] = {"W", "E"};
  const char *vrtxLabel[2] = {"X", "Y"};
  for (int iarm = 0; iarm < 2; iarm++)
  {
    for (int ivrtx = 0; ivrtx < 2; ivrtx++)
    {
      sprintf(draw, "vtx%s[%i]%+f:vtx[2]>>htmp(500,-10,10,1000,-0.05,0.05)",
              armLabel[iarm], ivrtx, -1 * bc[ivrtx]);
      cout << " " << draw << endl;

      ntp_event->Draw(
        // Form("vtx%s[%i]%+f:vtx[2]>>htmp(500,-10,10,1000,-0.1,0.1",
        //      armLabel[iarm], ivrtx, bc[ivrtx]),
        draw,
        "which_vtx==\"SVX_PRECISE\"",
        "goff");
      hvrtx_zvrtx[iarm][ivrtx] = (TH2D *) gDirectory->FindObject("htmp");
      hvrtx_zvrtx[iarm][ivrtx]->SetDirectory(0);
      hvrtx_zvrtx[iarm][ivrtx]->SetName(Form("hvrtx_zvrtx_%i_%i", iarm, ivrtx));
      hvrtx_zvrtx[iarm][ivrtx]->SetTitle(
        Form(";vtx Z [cm];%s vtx %s %+0.4f",
             armLabel[iarm], vrtxLabel[ivrtx], bc[ivrtx]));

      pvrtx_zvrtx[iarm][ivrtx] = new TProfile(
        Form("pvrtx_zvrtx_%i_%i", iarm, ivrtx),
        Form(";vtx Z [cm];%s vtx %s %+0.4f",
             armLabel[iarm], vrtxLabel[ivrtx], bc[ivrtx]),
        500, -10, 10,
        -0.1, 0.1);

      sprintf(draw, "vtx%s[%i]%+f:vtx[2]>>pvrtx_zvrtx_%i_%i",
              armLabel[iarm], ivrtx, -1 * bc[ivrtx], iarm, ivrtx);
      cout << " " << draw << endl;

      ntp_event->Draw(
        draw,
        "which_vtx==\"SVX_PRECISE\"",
        "goff");
      pvrtx_zvrtx[iarm][ivrtx]->SetLineColor(kRed);

      pvrtx_zvrtx[iarm][ivrtx]->Fit(fline, "RQ0");

      slope[iarm][ivrtx] = fline->GetParameter(1);
      slope_e[iarm][ivrtx] = fline->GetParError(1);

    }
  }

  //=============================================================//
  // PLOT OBJECTS
  //=============================================================//
  cout << endl;
  cout << "--> Plotting" << endl;

  TLatex lslope;
  lslope.SetNDC();
  lslope.SetTextAlign(22);

  //=============================================================//
  // PLOT
  //=============================================================//

  TCanvas *cvrtx = new TCanvas("cvrtx", "vrtx", 1200, 800);
  cvrtx->Divide(2, 2);
  for (int iarm = 0; iarm < 2; iarm++)
  {
    for (int ivrtx = 0; ivrtx < 2; ivrtx++)
    {
      cvrtx->cd(iarm * 2 + ivrtx + 1);
      hvrtx_zvrtx[iarm][ivrtx]->Draw("colz");
      pvrtx_zvrtx[iarm][ivrtx]->Draw("same");
      lslope.DrawLatex(0.5, 0.95,
                       Form("slope = %.4f #pm %.4f", slope[iarm][ivrtx], slope_e[iarm][ivrtx]));
    }
  }

  //=============================================================//
  // PRINT PLOTS
  //=============================================================//
  if (print_plots)
  {
    cout << endl;
    cout << "--> Printing Plots" << endl;

    cvrtx->Print("pdfs/look_xyvrtx_zdep_vrtx.pdf");
  }



}