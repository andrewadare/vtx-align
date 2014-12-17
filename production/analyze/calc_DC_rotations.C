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
    "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-3-2/testvtxproduction_fieldon-407951-3-2.root";
  // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-3-2_2/testvtxproduction_fieldon-407951-3-2_2.root";
  // "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-3-2_4/testvtxproduction_fieldon-407951-3-2_4.root";
  bool print_plots = false;
  const char *plot_tag = "fieldon-407951-3-2_phi";

  const int NLAYER = 4;
  const int NARM = 2;
  const int NPHI = 4;
  float phil[NPHI] = { -1.00, 0.25, 1.50, 3.25};
  float phih[NPHI] = { 0.25, 1.50, 3.25, 4.00};

  int color[4] =
  {
    kBlack,
    kBlue,
    kRed,
    kGreen + 2,
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
  TH3D *hdphiLayerPhi;
  TH3D *hdtheLayerPhi;

  TH3D *hsphiLayer[NPHI];
  TH3D *hzcotthetaLayer[NPHI];

  TProfile *ps_phi[NPHI][NLAYER];
  TProfile *pz_cottheta[NPHI][NLAYER];

  TH2D *hsphi_layer[NPHI][NLAYER];
  TH2D *hzcottheta_layer[NPHI][NLAYER];

  TH1D *hdphi_arm_layer[NPHI][NLAYER];
  TH1D *hdthe_arm_layer[NPHI][NLAYER];

  TGraphErrors *gdphi_arm[NPHI];
  TGraphErrors *gdthe_arm[NPHI];

  TF1 *fl = new TF1("fl", "pol0", 0, 4);
  fl->SetLineStyle(2);

  double dphi[NPHI] = {0};
  double dthe[NPHI] = {0};

  double dphi_e[NPHI] = {0};
  double dthe_e[NPHI] = {0};

  TF1 *fgaus = new TF1("fgaus", "gaus", -1, 1);
  fgaus->SetLineColor(kBlack);

  TF1 *fress = new TF1("fress",
                       "[0] + [1]*TMath::Sin(x) + [2]*TMath::Cos(x)",
                       -1, 4);
  fress->SetLineColor(kBlack);
  fress->SetLineStyle(2);
  fress->SetParLimits(0, -0.02, 0.02);
  fress->SetParLimits(1, -0.02, 0.02);
  fress->SetParLimits(2, -0.02, 0.02);

  //==========================================================================//
  // DEFINE HISTOGRAMS
  //==========================================================================//
  cout << endl;
  cout << "--> Defining histograms" << endl;

  hdphiLayerPhi = new TH3D("hdphiLayerPhi",
                           ";d#phi;Layer;#phi",
                           1000, -0.5, 0.5,
                           4, -0.5, 3.5,
                           75, -1, 4);
  hdtheLayerPhi = new TH3D("hdtheLayerPhi",
                           ";d#theta;Layer;#phi",
                           1000, -0.5, 0.5,
                           4, -0.5, 3.5,
                           75, -1, 4);

  for (int iarm = 0; iarm < NARM; iarm++)
  {
    hsphiLayer[iarm] = new TH3D(Form("hsphiLayer_%i", iarm),
                                ";PHCNT #phi^{0};#Deltas/r;Layer",
                                500, -1, 4,
                                500, -0.2, 0.2,
                                4, -0.5, 3.5);
    hzcotthetaLayer[iarm] = new TH3D(Form("hzcotthetaLayer_%i", iarm),
                                     ";PHCNT cot(#theta^{0});#Deltaz/r;Layer",
                                     100, -0.5, 0.5,
                                     500, -0.3, 0.3,
                                     4, -0.5, 3.5);

    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      ps_phi[iarm][ilayer] =
        new TProfile(Form("ps_phi_%i_%i", iarm, ilayer),
                     ";PHCNT #phi^{0};#Deltas/r",
                     500, -4, 7,
                     -0.3, 0.3);
      ps_phi[iarm][ilayer]->SetLineColor(kRed);

      pz_cottheta[iarm][ilayer] =
        new TProfile(Form("pz_cottheta_%i_%i", iarm, ilayer),
                     ";PHCNT cot(#theta^{0});#Deltaz/r;",
                     500, -4, 7,
                     -0.3, 0.3);
      pz_cottheta[iarm][ilayer]->SetLineColor(kRed);

    }
  }

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
        hdphiLayerPhi->Fill(TMath::ASin(res_s[iclus] / r),
                            layer[iclus],
                            phi0);
        hdtheLayerPhi->Fill(TMath::ASin(res_z[iclus] / r),
                            layer[iclus],
                            phi0);

        hsphiLayer[dcarm]->Fill(phi0,
                                res_s[iclus] / r,
                                layer[iclus]);

        ps_phi[dcarm][layer[iclus]]->Fill(phi0,
                                          res_s[iclus] / r);

        hzcotthetaLayer[dcarm]->Fill(1. / TMath::Tan(the0),
                                     res_z[iclus] / r,
                                     layer[iclus]);

        pz_cottheta[dcarm][layer[iclus]]->Fill(1. / TMath::Tan(the0),
                                               res_z[iclus] / r);

      }
    }
  }

  //-- Projections
  for (int iarm = 0; iarm < NARM; iarm++)
  {
    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      //-- project to 2D
      cout << "    project s arm" << iarm << " B" << ilayer << endl;
      hsphiLayer[iarm]->GetZaxis()->SetRange(ilayer + 1, ilayer + 1);
      hsphi_layer[iarm][ilayer] = (TH2D *) hsphiLayer[iarm]->Project3D("yx");
      hsphi_layer[iarm][ilayer]->SetName(Form("hsphi_%i_%i", iarm, ilayer));

      hzcotthetaLayer[iarm]->GetZaxis()->SetRange(ilayer + 1, ilayer + 1);
      hzcottheta_layer[iarm][ilayer] = (TH2D *)
                                       hzcotthetaLayer[iarm]->Project3D("yx");
      hzcottheta_layer[iarm][ilayer]->SetName(Form("hzcottheta_%i_%i", iarm, ilayer));
    }
  }

  //==========================================================================//
  // CALCULATE ROTATIONS
  //==========================================================================//
  cout << endl;
  cout << "--> Projecting and finding rotations" << endl;

  for (int iphi = 0; iphi < NPHI; iphi++)
  {
    gdphi_arm[iphi] = new TGraphErrors();
    gdphi_arm[iphi]->SetLineColor(color[iphi]);
    gdphi_arm[iphi]->SetMarkerColor(color[iphi]);
    gdphi_arm[iphi]->SetMarkerStyle(20 + iphi);
    gdphi_arm[iphi]->SetTitle(";Layer;d#phi");

    gdthe_arm[iphi] = new TGraphErrors();
    gdthe_arm[iphi]->SetLineColor(color[iphi]);
    gdthe_arm[iphi]->SetMarkerColor(color[iphi]);
    gdthe_arm[iphi]->SetMarkerStyle(20 + iphi);
    gdthe_arm[iphi]->SetTitle(";Layer;d#theta");



    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      //-- dphi
      int bl = hdphiLayerPhi->GetZaxis()->FindBin(phil[iphi]);
      int bh = hdphiLayerPhi->GetZaxis()->FindBin(phih[iphi]) - 1;
      hdphi_arm_layer[iphi][ilayer] = (TH1D *)
                                      hdphiLayerPhi->ProjectionX(
                                        Form("hdphi_%i_%i", iphi, ilayer),
                                        ilayer + 1, ilayer + 1,
                                        bl, bh);
      hdphi_arm_layer[iphi][ilayer]->SetLineColor(color[iphi]);

      float mean = hdphi_arm_layer[iphi][ilayer]->GetMean();
      float rms = hdphi_arm_layer[iphi][ilayer]->GetRMS();

      float ll = mean - 0.75 * rms;
      float hl = mean + 0.75 * rms;

      fgaus->SetRange(ll, hl);
      fgaus->SetParameters(hdphi_arm_layer[iphi][ilayer]->GetMaximum(),
                           mean,
                           rms);

      hdphi_arm_layer[iphi][ilayer]->Fit(fgaus, "RQ0");

      cout << "   mdphi[" << iphi << "][" << ilayer << "] = "
           << fgaus->GetParameter(1)
           << " +/- " << fgaus->GetParError(1)
           << endl;


      gdphi_arm[iphi]->SetPoint(ilayer,
                                ilayer,
                                fgaus->GetParameter(1));
      gdphi_arm[iphi]->SetPointError(ilayer,
                                     0,
                                     fgaus->GetParError(1));



      //-- dtheta
      bl = hdtheLayerPhi->GetZaxis()->FindBin(phil[iphi]);
      bh = hdtheLayerPhi->GetZaxis()->FindBin(phih[iphi]) - 1;
      hdthe_arm_layer[iphi][ilayer] = (TH1D *)
                                      hdtheLayerPhi->ProjectionX(
                                        Form("hdthe_%i_%i", iphi, ilayer),
                                        ilayer + 1, ilayer + 1,
                                        bl, bh);
      hdthe_arm_layer[iphi][ilayer]->SetLineColor(color[iphi]);

      mean = hdthe_arm_layer[iphi][ilayer]->GetMean();
      rms = hdthe_arm_layer[iphi][ilayer]->GetRMS();

      ll = mean - 0.75 * rms;
      hl = mean + 0.75 * rms;

      fgaus->SetRange(ll, hl);
      fgaus->SetParameters(hdthe_arm_layer[iphi][ilayer]->GetMaximum(),
                           mean,
                           rms);

      hdthe_arm_layer[iphi][ilayer]->Fit(fgaus, "RQ0");

      cout << "   mdthe[" << iphi << "][" << ilayer << "] = "
           << fgaus->GetParameter(1)
           << " +/- " << fgaus->GetParError(1)
           << endl;


      gdthe_arm[iphi]->SetPoint(ilayer,
                                ilayer,
                                fgaus->GetParameter(1));
      gdthe_arm[iphi]->SetPointError(ilayer,
                                     0,
                                     fgaus->GetParError(1));
    }

    //-- Fit phi
    gdphi_arm[iphi]->Fit("fl", "RQN");

    dphi[iphi] = fl->GetParameter(0);
    dphi_e[iphi] = fl->GetParError(0);

    cout << " " << phil[iphi] << " < phi < " << phih[iphi] << " dphi = "
         << dphi[iphi]
         << " +/- " << dphi_e[iphi]
         << endl;

    //-- Fit theta
    gdthe_arm[iphi]->Fit("fl", "RQN");

    dthe[iphi] = fl->GetParameter(0);
    dthe_e[iphi] = fl->GetParError(0);

    cout << " " << phil[iphi] << " < phi < " << phih[iphi] << " dthe = "
         << dthe[iphi]
         << " +/- " << dthe_e[iphi]
         << endl;

  }


  //==========================================================================//
  // FIT S/R VS PHI PLOTS
  //==========================================================================//
  cout << endl;
  cout << "--> Fitting s/r vs phi plots" << endl;

  for (int iarm = 0; iarm < NARM; iarm++)
  {
    cout << endl;
    cout << "----------- Arm=" << iarm << " -----------" << endl;
    if (iarm == 0)
      fress->SetRange(2.5, 3.25);
      // fress->SetRange(2, 4);
    else
      fress->SetRange(-1, 1.5);

    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      fress->SetParameters(dphi[NPHI - 2 * iarm - 1],
                           0.0,
                           0.0);
      ps_phi[iarm][ilayer]->Fit(fress, "R0");

    }
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

      ltitle.DrawLatex(0.5, 0.95, Form("Arm%i B%i", iarm, ilayer));
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

      ltitle.DrawLatex(0.5, 0.95, Form("Arm%i B%i", iarm, ilayer));
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
  haxis->SetMinimum(-0.1);
  haxis->SetMaximum( 0.1);
  haxis->DrawCopy();
  for (int i = 0; i < NPHI; ++i)
  {
    gdphi_arm[i]->Draw("P");

    fl->SetLineColor(color[i]);
    fl->SetParameter(0, dphi[i]);
    fl->DrawCopy("same");

    lmean.SetTextColor(color[i]);
    lmean.DrawLatex(0.2, 0.85 - 0.08 * i,
                    Form("%.2f<#phi<%.2f d#phi=%.4f",
                         phil[i], phih[i], dphi[i]));
  }

  cres->cd(2);
  haxis->SetTitle(";Layer;d#theta");
  haxis->SetMinimum(-0.1);
  haxis->SetMaximum( 0.1);
  haxis->DrawCopy();
  for (int i = 0; i < NPHI; ++i)
  {
    gdthe_arm[i]->Draw("P");

    fl->SetLineColor(color[i]);
    fl->SetParameter(0, dthe[i]);
    fl->DrawCopy("same");

    lmean.SetTextColor(color[i]);
    lmean.DrawLatex(0.2, 0.85 - 0.08 * i,
                    Form("%.2f<#phi<%.2f d#theta=%.4f",
                         phil[i], phih[i], dthe[i]));
  }


  TCanvas *csphi[NARM];
  for (int iarm = 0; iarm < NARM; iarm++)
  {
    csphi[iarm] = new TCanvas(Form("csphi_%i", iarm),
                              Form("s vs phi Arm%i", iarm),
                              1200, 700);
    csphi[iarm]->Divide(2, 2);
    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      csphi[iarm]->cd(ilayer + 1);
      if (iarm == 0)
        hsphi_layer[iarm][ilayer]->GetXaxis()->SetRangeUser(2, 4);
      else
        hsphi_layer[iarm][ilayer]->GetXaxis()->SetRangeUser(-1, 1.5);
      hsphi_layer[iarm][ilayer]->GetYaxis()->SetRangeUser(-0.05, 0.05);
      hsphi_layer[iarm][ilayer]->Draw("colz");
      ps_phi[iarm][ilayer]->Draw("same");

      ftmp = (TF1 *) ps_phi[iarm][ilayer]->GetFunction("fress");
      ftmp->DrawCopy("same");

      ltitle.DrawLatex(0.5, 0.95, Form("Arm%i B%i", iarm, ilayer));
    }
  }


  TCanvas *czthe[NARM];
  for (int iarm = 0; iarm < NARM; iarm++)
  {
    czthe[iarm] = new TCanvas(Form("czthe_%i", iarm),
                              Form("z vs cot(the0) Arm%i", iarm),
                              1200, 700);
    czthe[iarm]->Divide(2, 2);
    for (int ilayer = 0; ilayer < NLAYER; ilayer++)
    {
      czthe[iarm]->cd(ilayer + 1);
      hzcottheta_layer[iarm][ilayer]->GetYaxis()->SetRangeUser(-0.05, 0.05);
      hzcottheta_layer[iarm][ilayer]->Draw("colz");
      pz_cottheta[iarm][ilayer]->Draw("same");

      ltitle.DrawLatex(0.5, 0.95, Form("Arm%i B%i", iarm, ilayer));
    }
  }


  //==========================================================================//
  // PRINT PLOTS
  //==========================================================================//
  if (print_plots)
  {
    cout << endl;
    cout << "--> Printing plots with tag " << plot_tag << endl;

    cdphi->Print(Form("pdfs/dphi_layer_arm_%s.pdf", plot_tag));
    cdthe->Print(Form("pdfs/dtheta_layer_arm_%s.pdf", plot_tag));
    cres->Print(Form("pdfs/phcnt_phi0the0_rotations_%s.pdf", plot_tag));

  }


}