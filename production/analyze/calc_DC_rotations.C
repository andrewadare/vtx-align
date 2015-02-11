////////////////////////////////////////////////////////////////////////////////
//
// Calculate the rotations of the PHCentralTrack phi0, the0 using the
// s and z residuals calculated from SvxCentralTrackReco
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 12-10-2014
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TCut.h>
#include <TLine.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TEventList.h>
#include <TMath.h>

#include <iostream>

using namespace std;

void calc_DC_rotations()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //==========================================================================//
  // SET RUNNING CONDITIONS
  //==========================================================================//

  //should really use zero field for this ...
  const char *fileName =
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-3-2/testvtxproduction_fieldon-407951-3-2.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-3-2_2/testvtxproduction_fieldon-407951-3-2_2.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-3-2_4/testvtxproduction_fieldon-407951-3-2_4.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-405869-22-11_1/testvtxproduction_fieldon-405869-22-11_1.root";
    // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-414431-22-11_2/testvtxproduction_fieldon-414431-22-11_2.root";
    "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-22-11_1/testvtxproduction_fieldon-407951-22-11_1.root";

  bool print_plots = true;
  const char *plot_tag = "fieldon-407951-22-11_1";

  // float plot_rot_min = 0.025; //for _1
  float plot_rot_min = 0.005; //for _2

  const int NLAYER = 4;
  const int NARM = 2;

  int color[2] =
  {
    kBlue,
    kRed,
  };

  TCut trackCuts = "(dchquality == 31 || dchquality == 63) && "
                   "Nclus > 2 && Nclus_layer[0]>0 && Nclus_layer[1]>0 && "
                   "chisq/ndf < 3 && "
                   "pT>1";
  TCut eventCuts = "!tickCut && "
                   "TMath::Abs(vtx[2])<10 && "
                   //"(bbcq[0]+bbcq[1])>400 && " //central
                   // "(bbcq[0]+bbcq[1])<181 && " //peripheral
                   "which_vtx==\"SVX_PRECISE\"";

  TCut armCut[2] = {"dcarm==0", "dcarm==1"};
  TCut layerCut[4] = {"layer==0", "layer==1", "layer==2", "layer==3"};


  //==========================================================================//
  // DECLARE VARIABLES
  //==========================================================================//

  //-- Tree & variables
  TTree *ntp_SVXCNT;
  float dchquality;
  int Nclus;
  int Nclus_layer[4];
  float chisq;
  float ndf;
  float pT;
  bool tickCut;
  float vtx[3];
  int dcarm;
  int layer[7];
  float clus_global_x[7];
  float clus_global_y[7];
  float clus_global_z[7];
  float res_s[7];
  float res_z[7];
  float phi0;
  float the0;

  //-- Histograms
  TH3D *hdphiLayerArm;
  TH3D *hdtheLayerArm;

  TH1D *hdphi_arm_layer[NARM][NLAYER];
  TH1D *hdthe_arm_layer[NARM][NLAYER];

  TGraphErrors *gdphi_arm[NARM];
  TGraphErrors *gdthe_arm[NARM];

  TF1 *fl = new TF1("fl", "pol0", 0, 4);
  fl->SetLineStyle(2);

  double dphi[NARM] = {0};
  double dthe[NARM] = {0};

  double dphi_e[NARM] = {0};
  double dthe_e[NARM] = {0};

  TF1 *fgaus = new TF1("fgaus", "gaus", -1, 1);
  fgaus->SetLineColor(kBlack);


  //==========================================================================//
  // DEFINE HISTOGRAMS
  //==========================================================================//
  cout << endl;
  cout << "--> Defining histograms" << endl;

  hdphiLayerArm = new TH3D("hdphiLayerArm",
                           ";d#phi;Layer;#phi",
                           1000, -0.5, 0.5,
                           4, -0.5, 3.5,
                           2, -0.5, 1.5);
  hdtheLayerArm = new TH3D("hdtheLayerArm",
                           ";d#theta;Layer;#phi",
                           1000, -0.5, 0.5,
                           4, -0.5, 3.5,
                           2, -0.5, 1.5);


  //==========================================================================//
  // GET DATA
  //==========================================================================//
  cout << endl;
  cout << "--> Getting data from " << fileName << endl;

  TFile *fin = TFile::Open(fileName);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << fileName << endl;
    return;
  }

  //-- Set up Tree
  ntp_SVXCNT = (TTree *) fin->Get("ntp_SVXCNT");
  if (!ntp_SVXCNT)
  {
    cout << "ERROR!! Unable to find ntp_SVXCNT in " << fileName << endl;
    return;
  }

  ntp_SVXCNT->SetBranchAddress("dchquality", &dchquality);
  ntp_SVXCNT->SetBranchAddress("Nclus", &Nclus);
  ntp_SVXCNT->SetBranchAddress("Nclus_layer", &Nclus_layer);
  ntp_SVXCNT->SetBranchAddress("chisq", &chisq);
  ntp_SVXCNT->SetBranchAddress("ndf", &ndf);
  ntp_SVXCNT->SetBranchAddress("pT", &pT);
  ntp_SVXCNT->SetBranchAddress("tickCut", &tickCut);
  ntp_SVXCNT->SetBranchAddress("vtx", &vtx);
  ntp_SVXCNT->SetBranchAddress("dcarm", &dcarm);
  ntp_SVXCNT->SetBranchAddress("layer", &layer);
  ntp_SVXCNT->SetBranchAddress("clus_global_x", &clus_global_x);
  ntp_SVXCNT->SetBranchAddress("clus_global_y", &clus_global_y);
  ntp_SVXCNT->SetBranchAddress("clus_global_z", &clus_global_z);
  ntp_SVXCNT->SetBranchAddress("res_s", &res_s);
  ntp_SVXCNT->SetBranchAddress("res_z", &res_z);
  ntp_SVXCNT->SetBranchAddress("phi0", &phi0);
  ntp_SVXCNT->SetBranchAddress("the0", &the0);


  //-- Loop over entries
  unsigned int NENTRIES = ntp_SVXCNT->GetEntries();
  cout << "--> Will loop over " << NENTRIES << " entries" << endl;
  for (int ientry = 0; ientry < NENTRIES; ientry++)
  {
    ntp_SVXCNT->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "----> Entry " << ientry << endl;

    if ((dchquality == 31 || dchquality == 63) &&
        Nclus >= 3 &&
        Nclus_layer[0] > 0 && Nclus_layer[1] > 0 &&
        chisq / ndf < 3 &&
        pT > 1 &&
        !tickCut &&
        TMath::Abs(vtx[2]) < 10)
    {

      for (int iclus = 0; iclus < Nclus; iclus++)
      {
        float r = TMath::Sqrt(TMath::Power(clus_global_x[iclus], 2) +
                              TMath::Power(clus_global_y[iclus], 2));
        hdphiLayerArm->Fill(TMath::ASin(res_s[iclus] / r),
                            layer[iclus],
                            dcarm);
        hdtheLayerArm->Fill(TMath::ASin(res_z[iclus] / r),
                            layer[iclus],
                            dcarm);

      }
    }
  }

  //==========================================================================//
  // CALCULATE ROTATIONS
  //==========================================================================//
  cout << endl;
  cout << "--> Projecting and finding rotations" << endl;

  for (int iarm = 0; iarm < NARM; iarm++)
  {
    gdphi_arm[iarm] = new TGraphErrors();
    gdphi_arm[iarm]->SetLineColor(color[iarm]);
    gdphi_arm[iarm]->SetMarkerColor(color[iarm]);
    gdphi_arm[iarm]->SetMarkerStyle(20 + iarm);
    gdphi_arm[iarm]->SetTitle(";Layer;d#phi");

    gdthe_arm[iarm] = new TGraphErrors();
    gdthe_arm[iarm]->SetLineColor(color[iarm]);
    gdthe_arm[iarm]->SetMarkerColor(color[iarm]);
    gdthe_arm[iarm]->SetMarkerStyle(20 + iarm);
    gdthe_arm[iarm]->SetTitle(";Layer;d#theta");



    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      //-- dphi
      hdphi_arm_layer[iarm][ilayer] = (TH1D *)
                                      hdphiLayerArm->ProjectionX(
                                        Form("hdphi_%i_%i", iarm, ilayer),
                                        ilayer + 1, ilayer + 1,
                                        iarm + 1, iarm + 1);
      hdphi_arm_layer[iarm][ilayer]->SetLineColor(color[iarm]);

      float mean = hdphi_arm_layer[iarm][ilayer]->GetMean();
      float rms = hdphi_arm_layer[iarm][ilayer]->GetRMS();

      float ll = mean - 0.75 * rms;
      float hl = mean + 0.75 * rms;

      fgaus->SetRange(ll, hl);
      fgaus->SetParameters(hdphi_arm_layer[iarm][ilayer]->GetMaximum(),
                           mean,
                           rms);

      hdphi_arm_layer[iarm][ilayer]->Fit(fgaus, "RQ0");

      cout << "   mdphi[" << iarm << "][" << ilayer << "] = "
           << fgaus->GetParameter(1)
           << " +/- " << fgaus->GetParError(1)
           << endl;


      gdphi_arm[iarm]->SetPoint(ilayer,
                                ilayer,
                                fgaus->GetParameter(1));
      gdphi_arm[iarm]->SetPointError(ilayer,
                                     0,
                                     fgaus->GetParError(1));



      //-- dtheta
      hdthe_arm_layer[iarm][ilayer] = (TH1D *)
                                      hdtheLayerArm->ProjectionX(
                                        Form("hdthe_%i_%i", iarm, ilayer),
                                        ilayer + 1, ilayer + 1,
                                        iarm + 1, iarm + 1);
      hdthe_arm_layer[iarm][ilayer]->SetLineColor(color[iarm]);

      mean = hdthe_arm_layer[iarm][ilayer]->GetMean();
      rms = hdthe_arm_layer[iarm][ilayer]->GetRMS();

      ll = mean - 0.75 * rms;
      hl = mean + 0.75 * rms;

      fgaus->SetRange(ll, hl);
      fgaus->SetParameters(hdthe_arm_layer[iarm][ilayer]->GetMaximum(),
                           mean,
                           rms);

      hdthe_arm_layer[iarm][ilayer]->Fit(fgaus, "RQ0");

      cout << "   mdthe[" << iarm << "][" << ilayer << "] = "
           << fgaus->GetParameter(1)
           << " +/- " << fgaus->GetParError(1)
           << endl;


      gdthe_arm[iarm]->SetPoint(ilayer,
                                ilayer,
                                fgaus->GetParameter(1));
      gdthe_arm[iarm]->SetPointError(ilayer,
                                     0,
                                     fgaus->GetParError(1));
    }

    //-- Fit phi
    gdphi_arm[iarm]->Fit("fl", "RQN");

    dphi[iarm] = fl->GetParameter(0);
    dphi_e[iarm] = fl->GetParError(0);

    cout << " Arm" << iarm
         << " dphi = " << dphi[iarm]
         << " +/- " << dphi_e[iarm]
         << endl;

    //-- Fit theta
    gdthe_arm[iarm]->Fit("fl", "RQN");

    dthe[iarm] = fl->GetParameter(0);
    dthe_e[iarm] = fl->GetParError(0);

    cout << " Arm" << iarm
         << " dthe = " << dthe[iarm]
         << " +/- " << dthe_e[iarm]
         << endl;

  }


  //==========================================================================//
  // PLOT OBJECTS
  //==========================================================================//
  cout << endl;
  cout << "--> Plotting" << endl;

  TLatex ltitle;
  ltitle.SetNDC();
  ltitle.SetTextAlign(22);

  TLatex lmean;
  lmean.SetNDC();

  TLine lm;
  lm.SetLineColor(kBlack);
  lm.SetLineStyle(2);

  TH1D *haxis = new TH1D("haxis", "", 4, -0.5, 3.5);

  TF1 *ftmp;

  //==========================================================================//
  // PLOT
  //==========================================================================//

  TCanvas *cdphi = new TCanvas("cdphi", "dphi", 1200, 700);
  cdphi->Divide(NLAYER, NARM);
  for (int iarm = 0; iarm < NARM; iarm++)
  {
    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      cdphi->cd(iarm * NLAYER + ilayer + 1);
      hdphi_arm_layer[iarm][ilayer]->GetXaxis()->SetRangeUser(-0.1, 0.1);
      hdphi_arm_layer[iarm][ilayer]->Draw();

      ftmp = hdphi_arm_layer[iarm][ilayer]->GetFunction("fgaus");
      ftmp->DrawCopy("same");

      ltitle.DrawLatex(0.5, 0.95, Form("%s Arm%i B%i", plot_tag, iarm, ilayer));
      lmean.DrawLatex(0.6, 0.8, Form("m = %.4f", ftmp->GetParameter(1)));

      lm.DrawLine(
        ftmp->GetParameter(1),
        0,
        ftmp->GetParameter(1),
        hdphi_arm_layer[iarm][ilayer]->GetMaximum());

    }
  }

  TCanvas *cdthe = new TCanvas("cdthe", "dthe", 1200, 700);
  cdthe->Divide(NLAYER, NARM);
  for (int iarm = 0; iarm < NARM; iarm++)
  {
    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      cdthe->cd(iarm * NLAYER + ilayer + 1);
      hdthe_arm_layer[iarm][ilayer]->GetXaxis()->SetRangeUser(-0.1, 0.1);
      hdthe_arm_layer[iarm][ilayer]->Draw();

      ftmp = hdthe_arm_layer[iarm][ilayer]->GetFunction("fgaus");
      ftmp->DrawCopy("same");

      ltitle.DrawLatex(0.5, 0.95, Form("%s Arm%i B%i", plot_tag, iarm, ilayer));
      lmean.DrawLatex(0.6, 0.8, Form("m = %.4f", ftmp->GetParameter(1)));

      lm.DrawLine(
        ftmp->GetParameter(1),
        0,
        ftmp->GetParameter(1),
        hdthe_arm_layer[iarm][ilayer]->GetMaximum());

    }
  }

  TCanvas *cres = new TCanvas("cres", "results", 1200, 700);
  cres->Divide(2, 1);

  cres->cd(1);
  haxis->SetTitle(";Layer;d#phi");
  haxis->SetMinimum(-1 * plot_rot_min);
  haxis->SetMaximum( plot_rot_min);
  haxis->DrawCopy();
  for (int i = 0; i < NARM; ++i)
  {
    gdphi_arm[i]->Draw("P");

    fl->SetLineColor(color[i]);
    fl->SetParameter(0, dphi[i]);
    fl->DrawCopy("same");

    lmean.SetTextColor(color[i]);
    lmean.DrawLatex(0.2, 0.85 - 0.08 * i,
                    Form("Arm%i d#phi=%.4f",
                         i, dphi[i]));
    ltitle.DrawLatex(0.5, 0.95, plot_tag);

  }

  cres->cd(2);
  haxis->SetTitle(";Layer;d#theta");
  haxis->SetMinimum(-1 * plot_rot_min);
  haxis->SetMaximum( plot_rot_min);
  haxis->DrawCopy();
  for (int i = 0; i < NARM; ++i)
  {
    gdthe_arm[i]->Draw("P");

    fl->SetLineColor(color[i]);
    fl->SetParameter(0, dthe[i]);
    fl->DrawCopy("same");

    lmean.SetTextColor(color[i]);
    lmean.DrawLatex(0.2, 0.85 - 0.08 * i,
                    Form("Arm%i d#theta=%.4f",
                         i, dthe[i]));
    ltitle.DrawLatex(0.5, 0.95, plot_tag);
  }



  //==========================================================================//
  // PRINT PLOTS
  //==========================================================================//
  if (print_plots)
  {
    cout << endl;
    cout << "--> Printing plots with tag " << plot_tag << endl;

    cdphi->Print(Form("pdfs/dcrot/dphi_layer_arm_%s.pdf", plot_tag));
    cdthe->Print(Form("pdfs/dcrot/dtheta_layer_arm_%s.pdf", plot_tag));
    cres->Print(Form("pdfs/dcrot/phcnt_phi0the0_rotations_%s.pdf", plot_tag));

  }


}