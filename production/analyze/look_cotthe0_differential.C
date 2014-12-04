////////////////////////////////////////////////////////////////////////////////
//
// Look at PHCentralTrack & VTX variables vs cot(the0)
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 12-1-2014
//
////////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLegend.h>

#include <iostream>

using namespace std;

void look_cotthe0_differential()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //==========================================================//
  // SET RUNNING CONDITIONS
  //==========================================================//

  const char *inFile = "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-1-3_7/testvtxproduction_fieldon-407951-1-3_7.root";
  const char *fileLabel = "407951-1-3_7";
  // const char *inFile = "/direct/phenix+prod01/phnxreco/millepede/ZF_production/zf-411768-1-3/testvtxproduction_zf-411768-1-3.root";
  // const char *fileLabel = "zf-407951-1-3";
  // const char *inFile = "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-1-3/testvtxproduction_fieldon-407951-1-3.root";
  // const char *fileLabel = "407951-1-3";

  const int printPlots = true;

  const int NCOT = 2;
  int colorCot[NCOT] =
  {
    kBlue,
    kRed,
  };

  const float cotthe0_cut = 0.04;

  const int NLAYERS = 4;
  int NLADDERS[NLAYERS] = {10, 20, 16, 24};

  const int NARM = 2; //separate E/W arms
  const char *armLabel[] = {"East", "West"};

  const int NVTX = 3;
  const char *vtxLabel[3] = {"x", "y", "z"};

  TCut trackCuts = "(dchquality == 31 || dchquality == 63) && "
                   "Nclus > 2 && Nclus_layer[0]>0 && Nclus_layer[1]>0 && "
                   // "!(the0>1.5 && the0<1.65) &&"
                   //"charge == 1 &&"
                   "pT>1 &&"
                   "chisq/ndf < 3";
  TCut eventCuts = "!tickCut && "
                   "TMath::Abs(vtx[2])<10 && "
                   //"(bbcq[0]+bbcq[1])>400 && " //central
                   // "(bbcq[0]+bbcq[1])<181 && " //peripheral
                   "which_vtx==\"SVX_PRECISE\"";

  TCut cotthe0Cut[NCOT] =
  {
    Form("TMath::Abs(1./TMath::Tan(the0))>%.2f", cotthe0_cut),
    Form("TMath::Abs(1./TMath::Tan(the0))<%.2f", cotthe0_cut),
  };

  TCut layerCut = "layer==3";
  TCut ladderCut = "ladder==7";

  cout << " eventCuts      = " << eventCuts.GetTitle() << endl;
  cout << " trackCuts      = " << trackCuts.GetTitle() << endl;
  for (int icot = 0; icot < NCOT; icot++)
    cout << " cotthe0Cut[" << icot << "]  = " << cotthe0Cut[icot].GetTitle() << endl;

  cout << " layerCut       = " << layerCut.GetTitle() << endl;
  cout << " ladderCut      = " << ladderCut.GetTitle() << endl;


  //==========================================================//
  // DECLARE VARIABLES
  //==========================================================//

  TTree *ntp_SVXCNT;

  TH2F *hreszcotthe0;

  TH1F *hpT_cotthe0[NCOT];
  TH1F *hphi0_cotthe0[NCOT];
  TH1F *hchisq_cotthe0[NCOT];
  TH1F *hvtx_cotthe0[NCOT][NVTX];
  TH2F *hclusxz_cotthe0[NCOT];
  TH2F *hclusrz_cotthe0[NCOT];

  float vtxLim[NVTX][2] =
  {
    {0.2, 0.4},
    {0, 0.1},
    { -12, 12},
  };

  //==========================================================//
  // GET TREE FROM FILE
  //==========================================================//
  cout << endl;
  cout << "--> Getting tree from " << inFile << endl;

  TFile *fin = TFile::Open(inFile);
  if (!fin)
  {
    cout << "ERROR!! Unable to open " << inFile << endl;
    return;
  }

  ntp_SVXCNT = (TTree *) fin->Get("ntp_SVXCNT");
  if (!ntp_SVXCNT)
  {
    cout << "ERROR!! Unable to get ntp_SVXCNT from " << inFile << endl;
    return;
  }

  //==========================================================//
  // PERFORM PROJECTIONS FROM THE TREE
  //==========================================================//
  cout << endl;
  cout << "--> Projecting information from ntp_SVXCNT" << endl;

  ntp_SVXCNT->Draw("res_z:1./TMath::Tan(the0)>>htmp(100,-0.5,0.5,500,-0.5,0.5)",
                   trackCuts && eventCuts && layerCut,
                   "goff");
  hreszcotthe0 = (TH2F *) gDirectory->FindObject("htmp");
  hreszcotthe0->SetDirectory(0);
  hreszcotthe0->SetName("hreszcotthe0");
  hreszcotthe0->SetTitle(";cot(the0); SvxCentralTrack z residual [cm]");

  for (int icot = 0; icot < NCOT; icot++)
  {

    ntp_SVXCNT->Draw("pT>>htmp(50,0,5)",
                     trackCuts && eventCuts && cotthe0Cut[icot] && layerCut,
                     "goff");
    hpT_cotthe0[icot] = (TH1F *) gDirectory->FindObject("htmp");
    hpT_cotthe0[icot]->SetDirectory(0);
    hpT_cotthe0[icot]->SetName(Form("hpT_cotthe0_%i", icot));
    hpT_cotthe0[icot]->SetLineColor(colorCot[icot]);
    hpT_cotthe0[icot]->SetTitle(";p_{T} [GeV/c]");

    ntp_SVXCNT->Draw("phi0>>htmp(500,-1,4)",
                     trackCuts && eventCuts && cotthe0Cut[icot] && layerCut,
                     "goff");
    hphi0_cotthe0[icot] = (TH1F *) gDirectory->FindObject("htmp");
    hphi0_cotthe0[icot]->SetDirectory(0);
    hphi0_cotthe0[icot]->SetName(Form("hphi0_cotthe0_%i", icot));
    hphi0_cotthe0[icot]->SetLineColor(colorCot[icot]);
    hphi0_cotthe0[icot]->SetTitle(";phi0");

    ntp_SVXCNT->Draw("chisq/ndf>>htmp(100,0,4)",
                     trackCuts && eventCuts && cotthe0Cut[icot] && layerCut,
                     "goff");
    hchisq_cotthe0[icot] = (TH1F *) gDirectory->FindObject("htmp");
    hchisq_cotthe0[icot]->SetDirectory(0);
    hchisq_cotthe0[icot]->SetName(Form("hchisq_cotthe0_%i", icot));
    hchisq_cotthe0[icot]->SetLineColor(colorCot[icot]);
    hchisq_cotthe0[icot]->SetTitle(";chisq/ndf");

    for (int ivtx = 0; ivtx < NVTX; ivtx++)
    {
      ntp_SVXCNT->Draw(Form("vtx[%i]>>htmp(120,%f,%f)", ivtx, vtxLim[ivtx][0], vtxLim[ivtx][1]),
                       trackCuts && eventCuts && cotthe0Cut[icot] && layerCut,
                       "goff");
      hvtx_cotthe0[icot][ivtx] = (TH1F *) gDirectory->FindObject("htmp");
      hvtx_cotthe0[icot][ivtx]->SetDirectory(0);
      hvtx_cotthe0[icot][ivtx]->SetName(Form("hvtx_cotthe0_%i_%i", icot, ivtx));
      hvtx_cotthe0[icot][ivtx]->SetLineColor(colorCot[icot]);
      hvtx_cotthe0[icot][ivtx]->SetTitle(Form(";%s vtx", vtxLabel[ivtx]));
    }

    // ntp_SVXCNT->Draw("clus_global_x:clus_global_z>>htmp(300,-15,15, 100,15.,16.5)",
    ntp_SVXCNT->Draw("clus_global_x:clus_global_z>>htmp(300,-15,15, 500,-50.,50.)",
                     trackCuts && eventCuts && cotthe0Cut[icot]
                     && layerCut && ladderCut,
                     "goff");
    hclusxz_cotthe0[icot] = (TH2F *) gDirectory->FindObject("htmp");
    hclusxz_cotthe0[icot]->SetDirectory(0);
    hclusxz_cotthe0[icot]->SetName(Form("hclusxz_cotthe0_%i", icot));
    hclusxz_cotthe0[icot]->SetTitle(";cluster global z position [cm]; cluster global x position [cm]");


    ntp_SVXCNT->Draw("TMath::Sqrt(TMath::Power(clus_global_x, 2)+TMath::Power(clus_global_y, 2)):clus_global_z>>htmp(300,-15,15, 500,0.,50.)",
                     trackCuts && eventCuts && cotthe0Cut[icot]
                     && layerCut && ladderCut,
                     "goff");
    hclusrz_cotthe0[icot] = (TH2F *) gDirectory->FindObject("htmp");
    hclusrz_cotthe0[icot]->SetDirectory(0);
    hclusrz_cotthe0[icot]->SetName(Form("hclusrz_cotthe0_%i", icot));
    hclusrz_cotthe0[icot]->SetTitle(";cluster global z position [cm]; cluster global r position [cm]");



  }

  //==========================================================//
  // PLOT OBJECTS
  //==========================================================//
  cout << endl;
  cout << "--> Plotting" << endl;


  TLine lcut;
  lcut.SetLineStyle(2);
  lcut.SetLineColor(kRed);
  lcut.SetLineWidth(2);

  TLatex ltitle;
  ltitle.SetTextAlign(22);
  ltitle.SetNDC();

  TLegend *leg = new TLegend(0.3,0.7,0.9,0.95);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  for (int icot = 0; icot < NCOT; icot++)
    leg->AddEntry(hphi0_cotthe0[icot],cotthe0Cut[icot].GetTitle(),"L");

  //==========================================================//
  // PLOT
  //==========================================================//

  TCanvas *ccotthe0 = new TCanvas("ccotthe0", "cot(the0)", 700, 700);
  ccotthe0->cd(1);
  hreszcotthe0->Draw("colz");
  ltitle.DrawLatex(05, 0.95, fileLabel);
  lcut.DrawLine(-1 * cotthe0_cut, -0.5, -1 * cotthe0_cut, 0.5);
  lcut.DrawLine(cotthe0_cut, -0.5, cotthe0_cut, 0.5);


  TCanvas *cdc = new TCanvas("cdc", "DC", 1400, 600);
  cdc->Divide(3, 1);

  cdc->cd(1);
  gPad->SetLogy();
  hpT_cotthe0[0]->Draw();
  for (int icot = 1; icot < NCOT; icot++)
    hpT_cotthe0[icot]->Draw("same");
  leg->Draw("same");

  cdc->cd(2);
  gPad->SetLogy();
  hphi0_cotthe0[0]->Draw();
  for (int icot = 1; icot < NCOT; icot++)
    hphi0_cotthe0[icot]->Draw("same");

  cdc->cd(3);
  gPad->SetLogy();
  hchisq_cotthe0[0]->Draw();
  for (int icot = 1; icot < NCOT; icot++)
    hchisq_cotthe0[icot]->Draw("same");



  TCanvas *cvtx = new TCanvas("cvtx", "vtx", 1200, 600);
  cvtx->Divide(NVTX, 1);
  for (int ivtx = 0; ivtx < NVTX; ivtx++)
  {
    cvtx->cd(ivtx + 1);
    hvtx_cotthe0[0][ivtx]->Draw();
    for (int icot = 0; icot < NCOT; icot++)
      hvtx_cotthe0[icot][ivtx]->Draw("same");
  }

  TCanvas *cclus = new TCanvas("cclus", "clus pos", 1200, 600);
  cclus->Divide(NCOT, 1);
  for (int icot = 0; icot < NCOT; icot++)
  {
    cclus->cd(icot + 1);
    hclusxz_cotthe0[icot]->Draw("colz");
    ltitle.DrawLatex(0.5, 0.95,
                     Form("B%cL%c %s",
                          layerCut.GetTitle()[strlen(layerCut.GetTitle())-1], 
                          ladderCut.GetTitle()[strlen(ladderCut.GetTitle())-1], 
                          cotthe0Cut[icot].GetTitle()));
  }

  TCanvas *cclusr = new TCanvas("cclusr", "clus pos", 1200, 600);
  cclusr->Divide(NCOT, 1);
  for (int icot = 0; icot < NCOT; icot++)
  {
    cclusr->cd(icot + 1);
    hclusrz_cotthe0[icot]->Draw("colz");
    ltitle.DrawLatex(0.5, 0.95,
                     Form("B%cL%c %s",
                          layerCut.GetTitle()[strlen(layerCut.GetTitle())-1], 
                          ladderCut.GetTitle()[strlen(ladderCut.GetTitle())-1], 
                          cotthe0Cut[icot].GetTitle()));
  }


  //==========================================================//
  // PRINT PLOTS
  //==========================================================//
  if (printPlots)
  {
    cout << endl;
    cout << "--> Printing plots" << endl;

    ccotthe0->Print("pdfs/look_cotthe0_resz_vs_cotthe0.pdf");
    cdc->Print("pdfs/look_cotthe0_trackvar.pdf");
    cvtx->Print("pdfs/look_cotthe0_primvtx.pdf");
    cclus->Print("pdfs/look_cotthe0_clusxz.pdf");
    cclusr->Print("pdfs/look_cotthe0_clusrz.pdf");
  }

}