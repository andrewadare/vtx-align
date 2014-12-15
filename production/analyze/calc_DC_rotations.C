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
  bool print_plots = true;
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


  const char *plot;

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

  ntp_SVXCNT = (TTree *) fin->Get("ntp_SVXCNT");
  if (!ntp_SVXCNT)
  {
    cout << "ERROR!! Unable to find ntp_SVXCNT in " << fileName << endl;
    return;
  }

  //-- First make an entry list to save time
  ntp_SVXCNT->Draw(">>cutentries", eventCuts && trackCuts);
  TEventList *cutentries = (TEventList *) gDirectory->FindObject("cutentries");
  ntp_SVXCNT->SetEventList(cutentries);


  plot = "phi0:layer:TMath::ASin(res_s/"
         "TMath::Sqrt(TMath::Power(clus_global_x,2) + TMath::Power(clus_global_y,2)))"
         " >> htmp(1000,-0.5,0.5, 4,-0.5,3.5, 75,-1,4)";
  ntp_SVXCNT->Draw(plot, "", "goff");
  hdphiLayerPhi = (TH3D *) gDirectory->FindObject("htmp");
  hdphiLayerPhi->SetDirectory(0);
  hdphiLayerPhi->SetName("hdphiLayerPhi");
  hdphiLayerPhi->SetTitle(";d#phi;Layer;#phi");

  plot = "phi0:layer:TMath::ASin(res_z/"
         "TMath::Sqrt(TMath::Power(clus_global_x,2) + TMath::Power(clus_global_y,2)))"
         " >> htmp(1000,-0.5,0.5, 4,-0.5,3.5, 75,-1,4)";
  ntp_SVXCNT->Draw(plot, "", "goff");
  hdtheLayerPhi = (TH3D *) gDirectory->FindObject("htmp");
  hdtheLayerPhi->SetDirectory(0);
  hdtheLayerPhi->SetName("hdtheLayerPhi");
  hdtheLayerPhi->SetTitle(";d#theta;Layer;Arm");


  //--residual vs phi/cottheta
  for (int iarm = 0; iarm < NARM; iarm++)
  {
    cout << " s arm" << iarm << endl;
    plot = "layer:res_s:phi0"
           " >> htmp(500,-1,4, 500,-0.2,0.2, 4,-0.5,3.5)";
    ntp_SVXCNT->Draw(plot, armCut[iarm], "goff");
    hsphiLayer[iarm] = (TH3D *) gDirectory->FindObject("htmp");
    hsphiLayer[iarm]->SetDirectory(0);
    hsphiLayer[iarm]->SetName(Form("hsphiLayer_%i", iarm));
    hsphiLayer[iarm]->SetTitle(";PHCNT #phi^{0};#Deltas;Layer");

    cout << " z arm" << iarm << endl;
    plot = "layer:res_z:1./TMath::Tan(the0)"
           " >> htmp(100,-0.5,0.5, 500,-0.3,0.3, 4,-0.5,3.5)";
    ntp_SVXCNT->Draw(plot, armCut[iarm], "goff");
    hzcotthetaLayer[iarm] = (TH3D *) gDirectory->FindObject("htmp");
    hzcotthetaLayer[iarm]->SetDirectory(0);
    hzcotthetaLayer[iarm]->SetName(Form("hzcotthetaLayer_%i", iarm));
    hzcotthetaLayer[iarm]->SetTitle(";PHCNT cot(#theta^{0});#Deltaz;Layer");

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

      //-- TProfile
      /*
      cout << "    TProfile s arm" << iarm << " B" << ilayer << endl;

      ps_phi[iarm][ilayer] = new TProfile(Form("ps_phi_%i_%i", iarm, ilayer),
                                          ";PHCNT #phi^{0};#Deltas",
                                          500, -4, 7,
                                          -0.3, 0.3);
      ntp_SVXCNT->Draw(Form("res_s;phi0 >> ps_phi_%i_%i", iarm, ilayer),
                       "",//eventCuts && trackCuts && armCut[iarm] && layerCut[ilayer],
                       "goff");
      cout << "   " << ps_phi[iarm][ilayer]->GetEntries() << endl;
      ps_phi[iarm][ilayer]->SetLineColor(kRed);



      pz_cottheta[iarm][ilayer] = new TProfile(Form("pz_cottheta_%i_%i", iarm, ilayer),
                                          ";PHCNT cot(#theta^{0});#Deltaz;",
                                          500, -4, 7,
                                          -0.3, 0.3);
      ntp_SVXCNT->Draw(Form("res_z:1./TMath::Tan(the0) >> pz_cottheta_%i_%i", iarm, ilayer),
                       eventCuts && trackCuts && armCut[iarm] && layerCut[ilayer],
                       "goff");
      pz_cottheta[iarm][ilayer]->SetLineColor(kRed);
      */
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
      hsphi_layer[iarm][ilayer]->Draw("colz");
      // ps_phi[iarm][ilayer]->Draw("same");

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
      hzcottheta_layer[iarm][ilayer]->Draw("colz");
      // pz_cottheta[iarm][ilayer]->Draw("same");

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