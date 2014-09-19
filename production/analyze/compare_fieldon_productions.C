///////////////////////////////////////////////////////////
//
// Compare the output from different productions
// output from TestVTXProduction module
// (different productions use different par files)
//
///////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 6-30-2014
//
///////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TMath.h>

#include <iostream>

using namespace std;

void compare_fieldon_productions()
{

    gStyle->SetOptStat(0);
    //gStyle->SetPalette(56, 0, 0.5);

    //================================================//
    // SET RUNNING CONDITIONS
    //================================================//

    bool print_plots = true;

    const int NFILES = 5;
    const char *fileName[] =
    {

        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/taebong-407951-p2-v8/ntuple/testvtxproduction_taebong-407951-p2-v8.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-ideal-407951-0-0/ntuple/testvtxproduction_ideal-407951-0-0.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-0-0/ntuple/testvtxproduction_407951-0-0.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-ideal-407951-0-1/ntuple/testvtxproduction_ideal-407951-0-1.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-0-1/ntuple/testvtxproduction_407951-0-1.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-ideal-407951-0-2/ntuple/testvtxproduction_ideal-407951-0-2.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-0-2/ntuple/testvtxproduction_407951-0-2.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-0-3/ntuple/testvtxproduction_407951-0-3.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-0-4/ntuple/testvtxproduction_407951-0-4.root",
        //"/direct/phenix+prod01/phnxreco/millepede/fieldon/fieldon-407951-ideal/ntuple/testvtxproduction_407951-ideal.root",
        "/phenix/prod01/phnxreco/millepede/fieldon/taebong-407951-p2-v8-Wpvtx/ntuple/testvtxproduction_407951-taebong-p2-v8-Wpvtx.root",
        //"/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-0-4-Wpvtx/ntuple/testvtxproduction_407951-0-4-Wpvtx.root",
        //"/phenix/prod01/phnxreco/millepede/fieldon/taebong-407951-p2-v8-Wpvtx-WW/ntuple/testvtxproduction_taebong-407951-p2-v8-Wpvtx-WW.root",
        //"/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-0-2-Wpvtx-WW/ntuple/testvtxproduction_407951-0-2-Wpvtx-WW.root",
        //"/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-0-4-Wpvtx-WW/ntuple/testvtxproduction_407951-0-2-Wpvtx-WW.root",
        "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-0-0-Wpvtx/ntuple/testvtxproduction_407951-0-0-Wpvtx.root",
        "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-1-0-Wpvtx/ntuple/testvtxproduction_407951-1-0-Wpvtx.root",
        "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-4-0-Wpvtx/ntuple/testvtxproduction_407951-4-0-Wpvtx.root",
        "/phenix/prod01/phnxreco/millepede/fieldon/fieldon-407951-5-0-Wpvtx/ntuple/testvtxproduction_407951-5-0-Wpvtx.root",
    };
    const char *fileLabel[] =
    {
        //"Taebong p2 v8",
        //"Andrew ideal-407951-0-0",
        //"407951-0-0",
        //"Andrew ideal-407951-0-1",
        //"407951-0-1",
        //"Andrew ideal-407951-0-2",
        //"407951-0-2",
        //"407951-0-3",
        //"407951-0-4",
        //"Ideal",
        "Taebong p2 v8 Wpvtx",
        //"407951-0-4 Wpvtx",
        //"Taebong Wpvtx-WW",
        //"407951-0-2 Wpvtx-WW",
        //"407951-0-4 Wpvtx-WW",
        "407951-0-0 Wpvtx",
        "407951-1-0 Wpvtx",
        "407951-4-0 Wpvtx",
        "407951-5-0 Wpvtx",
    };
    int color[] =
    {
        kBlack,
        kBlue,
        kRed,
        kGreen + 2,
        kOrange,
    };
    int mstyle[] =
    {
        20,
        21,
        22,
        23,
        24,
    };

    const int NPTSVXCNT = 6;
    float ptl_svxcnt[] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0};
    float pth_svxcnt[] = {1.5, 2.0, 2.5, 3.0, 4.0, 5.0};

    const int NLAYERS = 4;
    int NLADDERS[NLAYERS] = {10, 20, 16, 24};

    const int NARM = 2; //separate E/W arms
    const char *armLabel[] = {"East", "West"};

    const char *trackCuts = "(dchquality == 31 || dchquality == 63) && "
                            "Nclus > 2 && Nclus_layer[0]>0 && Nclus_layer[1]>0 && "
                            "chisq/ndf < 3 && "
                            //"charge == 1 &&"
                            "pT>1";
    const char *eventCuts = "!tickCut && "
                            "TMath::Abs(vtx[2])<10 && "
                            //"(bbcq[0]+bbcq[1])>400 && " //central
                            //"(bbcq[0]+bbcq[1])<181 && " //peripheral
                            "which_vtx==\"SVX_PRECISE\"";
    const char *standaloneCuts = "Nclus_layer[0]>0 && Nclus_layer[1]>0 && pT>0.5";

    cout << " eventCuts      = " << eventCuts << endl;
    cout << " trackCuts      = " << trackCuts << endl;
    cout << " standaloneCuts = " << standaloneCuts << endl;

    //for phcentraltracks
    const char *dchtrackCuts = "(dchquality == 31 || dchquality == 63) && "
                               "pT>1";
    const char *dcheventCuts = "TMath::Abs(vtx[2])<10";


    //================================================//
    // DELCARE VARIABLES
    //================================================//

    TFile *fin;
    TTree *ntp_SVXCNT;
    TTree *ntp_SEG;
    TTree *ntp_event;
    TTree *ntp_CNT;

    TH2F *hdcapt_svxcnt[NARM][NFILES];
    TH2F *hdcaphi_svxcnt[NFILES];
    TProfile *pdcaphi_svxcnt[NFILES];

    TH2F *hdcazpt_svxcnt[NARM][NFILES];
    TH2F *hdcazphi_svxcnt[NFILES];
    TProfile *pdcazphi_svxcnt[NFILES];

    TH2F *hdcaphi_standalone[NFILES];
    TProfile *pdcaphi_standalone[NFILES];

    TH1F *hdca_svxcnt[NARM][NFILES];
    TH1F *hdca_standalone[NFILES];

    TH1F *hntracksphi_svxcnt[NFILES];
    TH1F *hntracksphi_cnt[NFILES];
    TH1F *hntracksphi_eff[NFILES];
    TH1F *hntracksphi_standalone[NFILES];
    float nevents[NFILES] = {0};

    TH1F *hvtx[NFILES][3];

    TH1F *hvtx_EW[NFILES][3];
    TH1F *hvtx_E[NFILES][3];
    TH1F *hvtx_W[NFILES][3];
    TH1F *hvtxzbbcz[NFILES];

    TH2F *hdcapt_standalone[NARM][NFILES];
    TH1F *hdca_pt_standalone[NARM][NFILES][NPTSVXCNT];

    TH1F *hdca_pt_svxcnt[NARM][NFILES][NPTSVXCNT];
    TH1F *hdcaz_pt_svxcnt[NARM][NFILES][NPTSVXCNT];

    TGraphErrors *gdca_ptres_svxcnt[NARM][NFILES];

    TH1F *hchisqndf_svxcnt[NARM][NFILES];
    TH1F *hchisqndf_standalone[NARM][NFILES];

    TH3F *hres_z_laylad[NFILES];
    TH3F *hres_s_laylad[NFILES];

    TH1F *hres_z[NFILES][NLAYERS][24];
    TH1F *hres_s[NFILES][NLAYERS][24];

    TGraphErrors *gres_z[NFILES];
    TGraphErrors *gres_s[NFILES];

    TH1F *hNclus_svxcnt[NARM][NFILES];
    TH1F *hNclus_standalone[NARM][NFILES];

    TF1 *fgaus = new TF1("fgaus", "gaus", -0.02, 0.02);
    fgaus->SetLineColor(kRed);

    float vtx_EW_mean[NFILES][3];

    //================================================//
    // READ FROM ROOT FILE
    //================================================//
    cout << endl;
    cout << "--> Reading root files " << endl;

    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        cout << "---> Reading " << fileName[ifile] << endl;
        fin = TFile::Open(fileName[ifile]);
        if (!fin)
        {
            cout << "ERROR!! Unable to open " << fileName[ifile] << endl;
            return;
        }

        //--- SvxCentralTrack ntuple ---//
        ntp_SVXCNT = (TTree *) fin->Get("ntp_SVXCNT");
        if (!ntp_SVXCNT)
        {
            cout << "ERROR!! Unable to get ntp_SVXCNT from " << fileName[ifile] << endl;
            return;
        }

        ntp_SVXCNT->Draw("pT:dca2D*charge>>htmp(200,-0.1,0.1, 50,0,5)",
                         Form("%s && %s && dcarm==0", eventCuts, trackCuts),
                         "goff");
        hdcapt_svxcnt[0][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcapt_svxcnt[0][ifile]->SetDirectory(0);
        hdcapt_svxcnt[0][ifile]->SetName(Form("hdcapt_svxcnt_0_%i", ifile));
        hdcapt_svxcnt[0][ifile]->SetTitle(";DCA 2D * chg [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("pT:dca2D*charge>>htmp(200,-0.1,0.1, 50,0,5)",
                         Form("%s && %s && dcarm==1", eventCuts, trackCuts),
                         "goff");
        hdcapt_svxcnt[1][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcapt_svxcnt[1][ifile]->SetDirectory(0);
        hdcapt_svxcnt[1][ifile]->SetName(Form("hdcapt_svxcnt_1_%i", ifile));
        hdcapt_svxcnt[1][ifile]->SetTitle(";DCA 2D * chg [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("dca2D*charge:phi0>>htmp(250,-1,4, 400,-0.06,0.06)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hdcaphi_svxcnt[ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcaphi_svxcnt[ifile]->SetDirectory(0);
        hdcaphi_svxcnt[ifile]->SetName(Form("hdcaphi_svxcnt_%i", ifile));
        hdcaphi_svxcnt[ifile]->SetTitle(";phi;DCA 2D * chg [cm]");
        hdcaphi_svxcnt[ifile]->GetYaxis()->SetTitleOffset(1.5);

        pdcaphi_svxcnt[ifile] = new TProfile(Form("pdcaphi_svxcnt_%i", ifile), "", 125, -1, 4);
        pdcaphi_svxcnt[ifile]->SetLineColor(kRed);
        ntp_SVXCNT->Draw(Form("dca2D*charge:phi0>>pdcaphi_svxcnt_%i", ifile),
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        pdcaphi_svxcnt[ifile]->SetDirectory(0);




        ntp_SVXCNT->Draw("pT:dcaz>>htmp(700,-0.7,0.7, 50,0,5)",
                         Form("%s && %s && dcarm==0", eventCuts, trackCuts),
                         "goff");
        hdcazpt_svxcnt[0][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcazpt_svxcnt[0][ifile]->SetDirectory(0);
        hdcazpt_svxcnt[0][ifile]->SetName(Form("hdcazpt_svxcnt_0_%i", ifile));
        hdcazpt_svxcnt[0][ifile]->SetTitle(";DCA Z [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("pT:dcaz>>htmp(700,-0.7,0.7, 50,0,5)",
                         Form("%s && %s && dcarm==1", eventCuts, trackCuts),
                         "goff");
        hdcazpt_svxcnt[1][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcazpt_svxcnt[1][ifile]->SetDirectory(0);
        hdcazpt_svxcnt[1][ifile]->SetName(Form("hdcazpt_svxcnt_1_%i", ifile));
        hdcazpt_svxcnt[1][ifile]->SetTitle(";DCA Z [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("dcaz:phi0>>htmp(250,-1,4, 200,-0.2,0.2)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hdcazphi_svxcnt[ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcazphi_svxcnt[ifile]->SetDirectory(0);
        hdcazphi_svxcnt[ifile]->SetName(Form("hdcazphi_svxcnt_%i", ifile));
        hdcazphi_svxcnt[ifile]->SetTitle(";phi;DCA Z [cm]");
        hdcazphi_svxcnt[ifile]->GetYaxis()->SetTitleOffset(1.5);

        pdcazphi_svxcnt[ifile] = new TProfile(Form("pdcazphi_svxcnt_%i", ifile), "", 125, -1, 4);
        pdcazphi_svxcnt[ifile]->SetLineColor(kRed);
        ntp_SVXCNT->Draw(Form("dcaz:phi0>>pdcazphi_svxcnt_%i", ifile),
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        pdcazphi_svxcnt[ifile]->SetDirectory(0);




        ntp_SVXCNT->Draw("phi0>>htmp(75,-1,4)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hntracksphi_svxcnt[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hntracksphi_svxcnt[ifile]->SetDirectory(0);
        hntracksphi_svxcnt[ifile]->SetName(Form("hntracksphi_svxcnt_%i", ifile));
        hntracksphi_svxcnt[ifile]->SetTitle(";phi;N SvxCentralTracks / event");
        hntracksphi_svxcnt[ifile]->GetYaxis()->SetTitleOffset(1.5);
        hntracksphi_svxcnt[ifile]->SetLineColor(color[ifile]);
        hntracksphi_svxcnt[ifile]->Sumw2();



        ntp_SVXCNT->Draw("chisq/ndf>>htmp(100,0,5)",
                         Form("%s && %s && dcarm==0", eventCuts, trackCuts),
                         "goff");
        hchisqndf_svxcnt[0][ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hchisqndf_svxcnt[0][ifile]->SetDirectory(0);
        hchisqndf_svxcnt[0][ifile]->SetName(Form("hchisqndf_svxcnt_0_%i", ifile));
        hchisqndf_svxcnt[0][ifile]->SetTitle(";#Chi^{2}/ndf");
        hchisqndf_svxcnt[0][ifile]->GetYaxis()->SetTitleOffset(1.5);
        hchisqndf_svxcnt[0][ifile]->SetLineColor(color[ifile]);


        ntp_SVXCNT->Draw("chisq/ndf>>htmp(100,0,5)",
                         Form("%s && %s && dcarm==1", eventCuts, trackCuts),
                         "goff");
        hchisqndf_svxcnt[1][ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hchisqndf_svxcnt[1][ifile]->SetDirectory(0);
        hchisqndf_svxcnt[1][ifile]->SetName(Form("hchisqndf_svxcnt_1_%i", ifile));
        hchisqndf_svxcnt[1][ifile]->SetTitle(";#Chi^{2}/ndf");
        hchisqndf_svxcnt[1][ifile]->GetYaxis()->SetTitleOffset(1.5);
        hchisqndf_svxcnt[1][ifile]->SetLineColor(color[ifile]);



        ntp_SVXCNT->Draw("layer:ladder:res_z>>htmp(500,-0.5,0.5,25,-0.5,24.5,4,-0.5,3.5)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hres_z_laylad[ifile] = (TH3F *) gDirectory->FindObject("htmp");
        hres_z_laylad[ifile]->SetDirectory(0);
        hres_z_laylad[ifile]->SetName(Form("hres_z_laylad_%i", ifile));
        hres_z_laylad[ifile]->SetTitle(";z residual [cm];ladder;layer");
        hres_z_laylad[ifile]->SetLineColor(color[ifile]);
        hres_z_laylad[ifile]->Sumw2();

        ntp_SVXCNT->Draw("layer:ladder:res_s>>htmp(800,-4,4,25,-0.5,24.5,4,-0.5,3.5)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hres_s_laylad[ifile] = (TH3F *) gDirectory->FindObject("htmp");
        hres_s_laylad[ifile]->SetDirectory(0);
        hres_s_laylad[ifile]->SetName(Form("hres_s_laylad_%i", ifile));
        hres_s_laylad[ifile]->SetTitle(";s residual [cm];ladder;layer");
        hres_s_laylad[ifile]->SetLineColor(color[ifile]);
        hres_s_laylad[ifile]->Sumw2();


        for (int iarm = 0; iarm < 2; iarm++)
        {
            ntp_SVXCNT->Draw("Nclus>>htmp(8,-0.5,7.5)",
                             Form("%s && %s && dcarm==%i", eventCuts, trackCuts, iarm),
                             "goff");
            hNclus_svxcnt[iarm][ifile] = (TH1F *) gDirectory->FindObject("htmp");
            hNclus_svxcnt[iarm][ifile]->SetDirectory(0);
            hNclus_svxcnt[iarm][ifile]->SetName(Form("hNclus_%i_%i", iarm, ifile));
            hNclus_svxcnt[iarm][ifile]->SetTitle(Form("DC Arm=%i;SvxCentralTracks N_{clusters}", iarm));
            hNclus_svxcnt[iarm][ifile]->SetLineColor(color[ifile]);
            hNclus_svxcnt[iarm][ifile]->Scale(1. /
                                              hNclus_svxcnt[iarm][ifile]->Integral(1,
                                                      hNclus_svxcnt[iarm][ifile]->GetNbinsX()));
        }

        //--- CNT ntuple ---//
        ntp_CNT = (TTree *) fin->Get("ntp_CNT");
        if (!ntp_CNT)
        {
            cout << "ERROR!! Unable to find ntp_CNT in " << fileName[ifile] << endl;
            return;
        }

        ntp_CNT->Draw("phi0>>htmp(75,-1,4)",
                      Form("%s && %s", dcheventCuts, dchtrackCuts),
                      "goff");
        hntracksphi_cnt[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hntracksphi_cnt[ifile]->SetDirectory(0);
        hntracksphi_cnt[ifile]->SetName(Form("hntracksphi_cnt_%i", ifile));
        hntracksphi_cnt[ifile]->SetTitle(";phi;N PHCentralTracks / event");
        hntracksphi_cnt[ifile]->GetYaxis()->SetTitleOffset(1.5);
        hntracksphi_cnt[ifile]->SetLineColor(color[ifile]);
        hntracksphi_cnt[ifile]->Sumw2();

        //calculate the SVXCentralTrack finding efficiency
        //assume binomial uncertainties
        //see http://home.fnal.gov/~paterno/images/effic.pdf
        //for reference
        hntracksphi_eff[ifile] = (TH1F *) hntracksphi_svxcnt[ifile]->Clone(Form("hntracksphi_eff_%i", ifile));
        hntracksphi_eff[ifile]->SetDirectory(0);
        hntracksphi_eff[ifile]->SetTitle(";phi;N SvxCentralTracks / N PHCentralTracks");
        hntracksphi_eff[ifile]->SetMarkerStyle(mstyle[ifile]);
        hntracksphi_eff[ifile]->SetMarkerColor(color[ifile]);
        for (int ix = 1; ix <= hntracksphi_eff[ifile]->GetNbinsX(); ix++)
        {
            float ncnt = hntracksphi_cnt[ifile]->GetBinContent(ix);
            float nsvxcnt = hntracksphi_svxcnt[ifile]->GetBinContent(ix);

            if (ncnt <= 0)
            {
                hntracksphi_eff[ifile]->SetBinContent(ix, 0);
                hntracksphi_eff[ifile]->SetBinError(ix, 0);
            }
            else if (nsvxcnt > ncnt)
            {
                cout << "  Huh... " << nsvxcnt << " > " << ncnt << " phi=" << hntracksphi_cnt[ifile]->GetBinCenter(ix) << endl;
                hntracksphi_eff[ifile]->SetBinContent(ix, 0);
                hntracksphi_eff[ifile]->SetBinError(ix, 0);

            }
            else
            {
                float eff = nsvxcnt / ncnt;
                float eff_e = 1. / ncnt * TMath::Sqrt(nsvxcnt * (1 - nsvxcnt / ncnt));

                hntracksphi_eff[ifile]->SetBinContent(ix, eff);
                hntracksphi_eff[ifile]->SetBinError(ix, eff_e);
            }
        }

        cout << " nsvxcnt = " << hntracksphi_svxcnt[ifile]->Integral(1, hntracksphi_svxcnt[ifile]->GetNbinsX()) << endl;
        cout << " cnt     = " << hntracksphi_cnt[ifile]->Integral(1, hntracksphi_cnt[ifile]->GetNbinsX()) << endl;
        cout << " eff     = " << hntracksphi_eff[ifile]->Integral(1, hntracksphi_eff[ifile]->GetNbinsX()) << endl;

        //--- Standalone ntuple ---//
        ntp_SEG = (TTree *) fin->Get("ntp_SEG");
        if (!ntp_SEG)
        {
            cout << "ERROR!! Unable to get ntp_SEG from " << fileName[ifile] << endl;
            return;
        }

        ntp_SEG->Draw("pT:dca2D*charge>>htmp(300,-0.3,0.3, 50,0,5)",
                      Form("%s && %s && TMath::Abs(TMath::ATan2(p[1],p[0]))>1.5", eventCuts, standaloneCuts),
                      "goff");
        hdcapt_standalone[0][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcapt_standalone[0][ifile]->SetDirectory(0);
        hdcapt_standalone[0][ifile]->SetName(Form("hdcapt_standalone_0_%i", ifile));
        hdcapt_standalone[0][ifile]->SetTitle(";SvxSegment DCA 2D * chg [cm]; p_{T} [GeV/c]");

        ntp_SEG->Draw("pT:dca2D*charge>>htmp(300,-0.3,0.3, 50,0,5)",
                      Form("%s && %s && TMath::Abs(TMath::ATan2(p[1],p[0]))<1.5", eventCuts, standaloneCuts),
                      "goff");
        hdcapt_standalone[1][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcapt_standalone[1][ifile]->SetDirectory(0);
        hdcapt_standalone[1][ifile]->SetName(Form("hdcapt_standalone_1_%i", ifile));
        hdcapt_standalone[1][ifile]->SetTitle(";SvxSegment DCA 2D * chg [cm]; p_{T} [GeV/c]");

        ntp_SEG->Draw("dca2D*charge>>htmp(350,-0.3,0.3)",
                      Form("%s && %s", eventCuts, standaloneCuts),
                      "goff");
        hdca_standalone[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hdca_standalone[ifile]->SetDirectory(0);
        hdca_standalone[ifile]->SetName(Form("hdca_standalone_%i", ifile));
        hdca_standalone[ifile]->SetTitle(";SvxSegment DCA 2D * chg [cm]");
        hdca_standalone[ifile]->SetLineColor(color[ifile]);
        hdca_standalone[ifile]->SetLineWidth(2);
        //hdca_standalone[ifile]->Scale(1. / hdca_standalone[ifile]->Integral(1, hdca_standalone[ifile]->GetNbinsX()));

        ntp_SEG->Draw("TMath::ATan2(p[1],p[0])>>htmp(200,-4,4)",
                      Form("%s && %s", eventCuts, standaloneCuts),
                      "goff");
        hntracksphi_standalone[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hntracksphi_standalone[ifile]->SetDirectory(0);
        hntracksphi_standalone[ifile]->SetName(Form("hntracksphi_standalone_%i", ifile));
        hntracksphi_standalone[ifile]->SetTitle(";phi;N SvxSegments / event");
        hntracksphi_standalone[ifile]->GetYaxis()->SetTitleOffset(1.5);
        hntracksphi_standalone[ifile]->SetLineColor(color[ifile]);
        hntracksphi_standalone[ifile]->Sumw2();
        cout << " nseg = " << hntracksphi_standalone[ifile]->Integral(1, hntracksphi_standalone[ifile]->GetNbinsX()) << endl;



        ntp_SEG->Draw("dca2D*charge:TMath::ATan2(p[1],p[0])>>htmp(400,-4,4, 200,-0.2,0.2)",
                      Form("%s && %s", eventCuts, standaloneCuts),
                      "goff");
        hdcaphi_standalone[ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcaphi_standalone[ifile]->SetDirectory(0);
        hdcaphi_standalone[ifile]->SetName(Form("hdcaphi_standalone_%i", ifile));
        hdcaphi_standalone[ifile]->SetTitle(";phi;SvxSegment DCA 2D * chg [cm]");
        hdcaphi_standalone[ifile]->GetYaxis()->SetTitleOffset(1.5);

        pdcaphi_standalone[ifile] = new TProfile(Form("pdcaphi_standalone_%i", ifile), "", 200, -4, 4);
        pdcaphi_standalone[ifile]->SetLineColor(kRed);
        ntp_SEG->Draw(Form("dca2D*charge:TMath::ATan2(p[1],p[0])>>pdcaphi_standalone_%i", ifile),
                      Form("%s && %s", eventCuts, standaloneCuts),
                      "goff");
        pdcaphi_standalone[ifile]->SetDirectory(0);




        ntp_SEG->Draw("chisq/ndf>>htmp(100,0,10)",
                      Form("%s && %s && TMath::Abs(TMath::ATan2(p[1],p[0]))>1.5", eventCuts, standaloneCuts),
                      "goff");
        hchisqndf_standalone[0][ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hchisqndf_standalone[0][ifile]->SetDirectory(0);
        hchisqndf_standalone[0][ifile]->SetName(Form("hchisqndf_standalone_0_%i", ifile));
        hchisqndf_standalone[0][ifile]->SetTitle(";SvxSegment #Chi^{2}/NDF");
        hchisqndf_standalone[0][ifile]->SetLineColor(color[ifile]);

        ntp_SEG->Draw("chisq/ndf>>htmp(100,0,10)",
                      Form("%s && %s && TMath::Abs(TMath::ATan2(p[1],p[0]))<1.5", eventCuts, standaloneCuts),
                      "goff");
        hchisqndf_standalone[1][ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hchisqndf_standalone[1][ifile]->SetDirectory(0);
        hchisqndf_standalone[1][ifile]->SetName(Form("hchisqndf_standalone_1_%i", ifile));
        hchisqndf_standalone[1][ifile]->SetTitle(";SvxSegment #Chi^{2}/NDF");
        hchisqndf_standalone[1][ifile]->SetLineColor(color[ifile]);


        ntp_SEG->Draw("Nclus>>htmp(8,-0.5,7.5)",
                      Form("%s && %s && TMath::Abs(TMath::ATan2(p[1],p[0]))<1.5", eventCuts, standaloneCuts),
                      "goff");
        hNclus_standalone[0][ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hNclus_standalone[0][ifile]->SetDirectory(0);
        hNclus_standalone[0][ifile]->SetName(Form("hNclus_%i_%i", 0, ifile));
        hNclus_standalone[0][ifile]->SetTitle(Form("DC Arm=%i;SvxSegments N_{clusters}", 0));
        hNclus_standalone[0][ifile]->SetLineColor(color[ifile]);
        hNclus_standalone[0][ifile]->Scale(1. /
                                           hNclus_standalone[0][ifile]->Integral(1,
                                                   hNclus_standalone[0][ifile]->GetNbinsX()));


        ntp_SEG->Draw("Nclus>>htmp(8,-0.5,7.5)",
                      Form("%s && %s && TMath::Abs(TMath::ATan2(p[1],p[0]))>1.5", eventCuts, standaloneCuts),
                      "goff");
        hNclus_standalone[1][ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hNclus_standalone[1][ifile]->SetDirectory(0);
        hNclus_standalone[1][ifile]->SetName(Form("hNclus_%i_%i", 1, ifile));
        hNclus_standalone[1][ifile]->SetTitle(Form("DC Arm=%i;SvxSegments N_{clusters}", 1));
        hNclus_standalone[1][ifile]->SetLineColor(color[ifile]);
        hNclus_standalone[1][ifile]->Scale(1. /
                                           hNclus_standalone[1][ifile]->Integral(1,
                                                   hNclus_standalone[1][ifile]->GetNbinsX()));



        //--- Event ntuple ---//
        ntp_event = (TTree *) fin->Get("ntp_event");
        if (!ntp_event)
        {
            cout << "ERROR!! Unable to get ntp_event from " << fileName[ifile] << endl;
            return;
        }

        ntp_event->Draw("run>>htmp", eventCuts, "goff");
        nevents[ifile] = ((TH1F *) gDirectory->FindObject("htmp"))->GetEntries();
        cout << "  Found " << ntp_event->GetEntries() << " events in " << fileName[ifile] << endl;
        cout << "  Found " << nevents[ifile] << " events which pass cuts in " << fileName[ifile] << endl;

        if (nevents[ifile] > 0)
        {
            hntracksphi_svxcnt[ifile]->Scale(1. / nevents[ifile]);
        }

        ntp_event->Draw("vtx[0]>>htmp(200,-0.2,0.0)", eventCuts, "goff");
        hvtx[ifile][0] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx[ifile][0]->SetDirectory(0);
        hvtx[ifile][0]->SetName(Form("vtx_%i_%i", ifile, 0));
        hvtx[ifile][0]->SetTitle(";X_{vtx}");
        hvtx[ifile][0]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtx[1]>>htmp(200,0.0,0.2)", eventCuts, "goff");
        hvtx[ifile][1] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx[ifile][1]->SetDirectory(0);
        hvtx[ifile][1]->SetName(Form("vtx_%i_%i", ifile, 1));
        hvtx[ifile][1]->SetTitle(";Y_{vtx}");
        hvtx[ifile][1]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtx[2]>>htmp(150,-15,15)", eventCuts, "goff");
        hvtx[ifile][2] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx[ifile][2]->SetDirectory(0);
        hvtx[ifile][2]->SetName(Form("vtx_%i_%i", ifile, 2));
        hvtx[ifile][2]->SetTitle(";Z_{vtx}");
        hvtx[ifile][2]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxE[0]>>htmp(200,-0.2,0.0)", eventCuts, "goff");
        hvtx_E[ifile][0] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_E[ifile][0]->SetDirectory(0);
        hvtx_E[ifile][0]->SetName(Form("vtx_%i_%i", ifile, 0));
        hvtx_E[ifile][0]->SetTitle(";E X_{vtx}");
        hvtx_E[ifile][0]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxE[1]>>htmp(200,0.0,0.2)", eventCuts, "goff");
        hvtx_E[ifile][1] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_E[ifile][1]->SetDirectory(0);
        hvtx_E[ifile][1]->SetName(Form("vtx_%i_%i", ifile, 1));
        hvtx_E[ifile][1]->SetTitle(";E Y_{vtx}");
        hvtx_E[ifile][1]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxE[2]>>htmp(150,-15,15)", eventCuts, "goff");
        hvtx_E[ifile][2] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_E[ifile][2]->SetDirectory(0);
        hvtx_E[ifile][2]->SetName(Form("vtx_%i_%i", ifile, 2));
        hvtx_E[ifile][2]->SetTitle(";E Z_{vtx}");
        hvtx_E[ifile][2]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxW[0]>>htmp(200,-0.2,0.0)", eventCuts, "goff");
        hvtx_W[ifile][0] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_W[ifile][0]->SetDirectory(0);
        hvtx_W[ifile][0]->SetName(Form("vtx_%i_%i", ifile, 0));
        hvtx_W[ifile][0]->SetTitle(";W X_{vtx}");
        hvtx_W[ifile][0]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxW[1]>>htmp(200,0.0,0.2)", eventCuts, "goff");
        hvtx_W[ifile][1] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_W[ifile][1]->SetDirectory(0);
        hvtx_W[ifile][1]->SetName(Form("vtx_%i_%i", ifile, 1));
        hvtx_W[ifile][1]->SetTitle(";W Y_{vtx}");
        hvtx_W[ifile][1]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxW[2]>>htmp(150,-15,15)", eventCuts, "goff");
        hvtx_W[ifile][2] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_W[ifile][2]->SetDirectory(0);
        hvtx_W[ifile][2]->SetName(Form("vtx_%i_%i", ifile, 2));
        hvtx_W[ifile][2]->SetTitle(";W Z_{vtx}");
        hvtx_W[ifile][2]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxW[0]-vtxE[0]>>htmp(200,-0.2,0.2)", eventCuts, "goff");
        hvtx_EW[ifile][0] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_EW[ifile][0]->SetDirectory(0);
        hvtx_EW[ifile][0]->SetName(Form("vtx_EW_%i_%i", ifile, 0));
        hvtx_EW[ifile][0]->SetTitle(";W-E X_{vtx}");
        hvtx_EW[ifile][0]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxW[1]-vtxE[1]>>htmp(200,-0.2,0.2)", eventCuts, "goff");
        hvtx_EW[ifile][1] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_EW[ifile][1]->SetDirectory(0);
        hvtx_EW[ifile][1]->SetName(Form("vtx_EW_%i_%i", ifile, 1));
        hvtx_EW[ifile][1]->SetTitle(";W-E Y_{vtx}");
        hvtx_EW[ifile][1]->SetLineColor(color[ifile]);


        ntp_event->Draw("vtxW[2]-vtxE[2]>>htmp(200,-0.2,0.2)", eventCuts, "goff");
        hvtx_EW[ifile][2] = (TH1F *) gDirectory->FindObject("htmp");
        hvtx_EW[ifile][2]->SetDirectory(0);
        hvtx_EW[ifile][2]->SetName(Form("vtx_EW_%i_%i", ifile, 2));
        hvtx_EW[ifile][2]->SetTitle(";W-E Z_{vtx}");
        hvtx_EW[ifile][2]->SetLineColor(color[ifile]);


        vtx_EW_mean[ifile][0] = hvtx_EW[ifile][0]->GetMean();
        vtx_EW_mean[ifile][1] = hvtx_EW[ifile][1]->GetMean();
        vtx_EW_mean[ifile][2] = hvtx_EW[ifile][2]->GetMean();


        ntp_event->Draw("vtx[2]-bbcz>>htmp(200,-10,10)", eventCuts, "goff");
        hvtxzbbcz[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hvtxzbbcz[ifile]->SetDirectory(0);
        hvtxzbbcz[ifile]->SetName(Form("hvtxzbbcz_%i", ifile));
        hvtxzbbcz[ifile]->SetTitle(";VTX_{Z} - BBC_{Z} [cm]");
        hvtxzbbcz[ifile]->SetLineColor(color[ifile]);

        //cleanup
        delete ntp_SVXCNT;
        delete ntp_SEG;
        delete ntp_event;
        fin->Close();
        delete fin;

    }


    //================================================//
    // PROJECT OF PT
    //================================================//
    cout << endl;
    cout << "--> Projecting DCA over pT bins" << endl;

    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        for (int iarm = 0; iarm < NARM; iarm++)
        {
            int bl = hdcapt_svxcnt[iarm][ifile]->GetYaxis()->FindBin(ptl_svxcnt[0]);
            int bh = hdcapt_svxcnt[iarm][ifile]->GetYaxis()->FindBin(pth_svxcnt[NPTSVXCNT - 1]);
            hdca_svxcnt[iarm][ifile] = (TH1F *) hdcapt_svxcnt[iarm][ifile]->ProjectionX(Form("hdca_svxcnt_%i", ifile), bl, bh);
            hdca_svxcnt[iarm][ifile]->SetTitle(Form("%.1f<p_{T}<%.1f;DCA 2D * chg [cm]", ptl_svxcnt[0], pth_svxcnt[NPTSVXCNT - 1]));
            hdca_svxcnt[iarm][ifile]->SetLineColor(color[ifile]);
            hdca_svxcnt[iarm][ifile]->Scale(1. / hdca_svxcnt[iarm][ifile]->Integral(1, hdca_svxcnt[iarm][ifile]->GetNbinsX()));

            gdca_ptres_svxcnt[iarm][ifile] = new TGraphErrors();
            gdca_ptres_svxcnt[iarm][ifile]->SetTitle(";p_{T} [GeV/c];DCA2D resolution [#mum]");
            gdca_ptres_svxcnt[iarm][ifile]->SetMarkerStyle(mstyle[ifile]);
            gdca_ptres_svxcnt[iarm][ifile]->SetLineColor(color[ifile]);
            gdca_ptres_svxcnt[iarm][ifile]->SetMarkerColor(color[ifile]);
            gdca_ptres_svxcnt[iarm][ifile]->GetYaxis()->SetTitleOffset(1.5);

            for (int ipt = 0; ipt < NPTSVXCNT; ipt++)
            {
                bl = hdcapt_svxcnt[iarm][ifile]->GetYaxis()->FindBin(ptl_svxcnt[ipt]);
                bh = hdcapt_svxcnt[iarm][ifile]->GetYaxis()->FindBin(pth_svxcnt[ipt]);
                hdca_pt_svxcnt[iarm][ifile][ipt] =
                    (TH1F *) hdcapt_svxcnt[iarm][ifile]->ProjectionX(
                        Form("hdca_pt_svxcnt_%i_%i_%i", iarm, ifile, ipt), bl, bh);
                hdca_pt_svxcnt[iarm][ifile][ipt]->SetTitle(
                    Form("DCArm%i %.1f<p_{T}<%.1f;SvxCNT DCA 2D * chg [cm];rescaled", iarm, ptl_svxcnt[ipt], pth_svxcnt[ipt]));
                hdca_pt_svxcnt[iarm][ifile][ipt]->SetLineColor(color[ifile]);
                hdca_pt_svxcnt[iarm][ifile][ipt]->Sumw2();
                hdca_pt_svxcnt[iarm][ifile][ipt]->Scale(
                    1. / hdca_pt_svxcnt[iarm][ifile][ipt]->GetMaximum());

                fgaus->SetParameters(hdca_pt_svxcnt[iarm][ifile][ipt]->GetMaximum(),
                                     hdca_pt_svxcnt[iarm][ifile][ipt]->GetMean(),
                                     hdca_pt_svxcnt[iarm][ifile][ipt]->GetRMS() );

                //fgaus->SetRange(hdca_pt_svxcnt[iarm][ifile][ipt]->GetMean() - 0.75 * hdca_pt_svxcnt[iarm][ifile][ipt]->GetRMS(),
                //                hdca_pt_svxcnt[iarm][ifile][ipt]->GetMean() + 0.75 * hdca_pt_svxcnt[iarm][ifile][ipt]->GetRMS() );

                hdca_pt_svxcnt[iarm][ifile][ipt]->Fit(fgaus, "R0Q");

                gdca_ptres_svxcnt[iarm][ifile]->SetPoint(ipt,
                        0.5 * (ptl_svxcnt[ipt] + pth_svxcnt[ipt]),
                        fgaus->GetParameter(2) * 1e4);
                gdca_ptres_svxcnt[iarm][ifile]->SetPointError(ipt,
                        0,
                        fgaus->GetParError(2) * 1e4);



                //-- dcaz
                bl = hdcazpt_svxcnt[iarm][ifile]->GetYaxis()->FindBin(ptl_svxcnt[ipt]);
                bh = hdcazpt_svxcnt[iarm][ifile]->GetYaxis()->FindBin(pth_svxcnt[ipt]);
                hdcaz_pt_svxcnt[iarm][ifile][ipt] =
                    (TH1F *) hdcazpt_svxcnt[iarm][ifile]->ProjectionX(
                        Form("hdcaz_pt_svxcnt_%i_%i_%i", iarm, ifile, ipt), bl, bh);
                hdcaz_pt_svxcnt[iarm][ifile][ipt]->SetTitle(
                    Form("DC Arm=%i %.1f<p_{T}<%.1f;SvxCNT DCA Z [cm];rescaled", iarm, ptl_svxcnt[ipt], pth_svxcnt[ipt]));
                hdcaz_pt_svxcnt[iarm][ifile][ipt]->SetLineColor(color[ifile]);
                hdcaz_pt_svxcnt[iarm][ifile][ipt]->Scale(
                    1. / hdcaz_pt_svxcnt[iarm][ifile][ipt]->GetMaximum());

                //-- dca2d standalone
                bl = hdcapt_standalone[iarm][ifile]->GetYaxis()->FindBin(ptl_svxcnt[ipt]);
                bh = hdcapt_standalone[iarm][ifile]->GetYaxis()->FindBin(pth_svxcnt[ipt]);
                hdca_pt_standalone[iarm][ifile][ipt] =
                    (TH1F *) hdcapt_standalone[iarm][ifile]->ProjectionX(
                        Form("hdca_pt_standalone_%i_%i_%i", iarm, ifile, ipt), bl, bh);
                hdca_pt_standalone[iarm][ifile][ipt]->SetTitle(
                    Form("DCArm%i %.1f<p_{T}<%.1f;SvxSegment DCA 2D * chg [cm]", iarm, ptl_svxcnt[ipt], pth_svxcnt[ipt]));
                hdca_pt_standalone[iarm][ifile][ipt]->SetLineColor(color[ifile]);
                hdca_pt_standalone[iarm][ifile][ipt]->Scale(
                    1. / hdca_pt_standalone[iarm][ifile][ipt]->GetMaximum());


            }
        }

    }


    //================================================//
    // CALCULATE RESIDUALS FOR EACH LADDER
    //================================================//
    cout << endl;
    cout << "--> Projecting residuals for each layer/ladder" << endl;

    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        gres_z[ifile] = new TGraphErrors();
        gres_z[ifile]->SetName(Form("gres_z_%i", ifile));
        gres_z[ifile]->SetLineColor(color[ifile]);
        gres_z[ifile]->SetMarkerColor(color[ifile]);
        gres_z[ifile]->SetMarkerStyle(mstyle[ifile]);

        gres_s[ifile] = new TGraphErrors();
        gres_s[ifile]->SetName(Form("gres_s_%i", ifile));
        gres_s[ifile]->SetLineColor(color[ifile]);
        gres_s[ifile]->SetMarkerColor(color[ifile]);
        gres_s[ifile]->SetMarkerStyle(mstyle[ifile]);

        int idx = 0;
        for (int ilayer = 0; ilayer < NLAYERS; ilayer++)
        {
            for (int iladder = 0; iladder < NLADDERS[ilayer]; iladder++)
            {
                //-- Z residual
                hres_z_laylad[ifile]->GetYaxis()->SetRange(iladder + 1, iladder + 1);
                hres_z_laylad[ifile]->GetZaxis()->SetRange(ilayer + 1, ilayer + 1);

                hres_z[ifile][ilayer][iladder] =
                    (TH1F *) hres_z_laylad[ifile]->Project3D("x");
                hres_z[ifile][ilayer][iladder]->SetName(Form("hres_z_%i_%i_%i", ifile, ilayer, iladder));
                hres_z[ifile][ilayer][iladder]->SetTitle(Form("B%iL%02i;z residual [cm]", ilayer, iladder));
                hres_z[ifile][ilayer][iladder]->Scale(
                    1. / hres_z[ifile][ilayer][iladder]->GetMaximum());

                if (hres_z[ifile][ilayer][iladder]->GetEntries() > 100)
                {
                    fgaus->SetRange(-0.5, 0.5);
                    fgaus->SetParameters(
                        hres_z[ifile][ilayer][iladder]->GetMaximum(),
                        hres_z[ifile][ilayer][iladder]->GetMean(),
                        hres_z[ifile][ilayer][iladder]->GetRMS());

                    hres_z[ifile][ilayer][iladder]->Fit(fgaus, "R0Q");

                    gres_z[ifile]->SetPoint(
                        idx,
                        idx,
                        fgaus->GetParameter(1));

                    gres_z[ifile]->SetPointError(
                        idx,
                        0,
                        fgaus->GetParameter(2));
                }
                else
                {
                    gres_z[ifile]->SetPoint(
                        idx,
                        idx,
                        -9999.);

                    gres_z[ifile]->SetPointError(
                        idx,
                        0,
                        0);

                }


                //-- S residual
                hres_s_laylad[ifile]->GetYaxis()->SetRange(iladder + 1, iladder + 1);
                hres_s_laylad[ifile]->GetZaxis()->SetRange(ilayer + 1, ilayer + 1);

                hres_s[ifile][ilayer][iladder] =
                    (TH1F *) hres_s_laylad[ifile]->Project3D("x");
                hres_s[ifile][ilayer][iladder]->SetName(Form("hres_s_%i_%i_%i", ifile, ilayer, iladder));
                hres_s[ifile][ilayer][iladder]->SetTitle(Form("B%iL%02i;s residual [cm]", ilayer, iladder));
                hres_s[ifile][ilayer][iladder]->Scale(
                    1. / hres_s[ifile][ilayer][iladder]->GetMaximum());


                if (hres_s[ifile][ilayer][iladder]->GetEntries() > 100)
                {
                    fgaus->SetRange(-4, 4);
                    fgaus->SetParameters(
                        hres_s[ifile][ilayer][iladder]->GetMaximum(),
                        hres_s[ifile][ilayer][iladder]->GetMean(),
                        hres_s[ifile][ilayer][iladder]->GetRMS());

                    hres_s[ifile][ilayer][iladder]->Fit(fgaus, "R0Q");

                    gres_s[ifile]->SetPoint(
                        idx,
                        idx,
                        fgaus->GetParameter(1));

                    gres_s[ifile]->SetPointError(
                        idx,
                        0,
                        fgaus->GetParameter(2));
                }
                else
                {
                    gres_s[ifile]->SetPoint(
                        idx,
                        idx,
                        -9999.);

                    gres_s[ifile]->SetPointError(
                        idx,
                        0,
                        0);
                }

                idx++;
            }
        }
    }
    //================================================//
    // PLOT OBJECTS
    //================================================//
    cout << endl;
    cout << "--> Plotting" << endl;


    TLegend *leg = new TLegend(0.6, 0.7, 0.98, 0.85);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        leg->AddEntry(hdca_svxcnt[0][ifile], fileLabel[ifile], "LP");
    }

    TLegend *legwevents = new TLegend(0.3, 0.6, 0.7, 0.9);
    legwevents->SetFillStyle(0);
    legwevents->SetBorderSize(0);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        legwevents->AddEntry(hdca_svxcnt[0][ifile], fileLabel[ifile], "L");
        legwevents->AddEntry((TObject *)0, Form("   %.0f events", nevents[ifile]), "");
    }

    TLegend *legvtxEW[3];
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        legvtxEW[ivtx] = new TLegend(0.5, 0.75, 1., 0.9);
        legvtxEW[ivtx]->SetFillStyle(0);
        legvtxEW[ivtx]->SetBorderSize(0);
        legvtxEW[ivtx]->SetTextSize(0.03);
        for (int ifile = 0; ifile < NFILES; ifile++)
        {
            legvtxEW[ivtx]->AddEntry(hvtx_EW[ifile][ivtx], fileLabel[ifile], "L");
            legvtxEW[ivtx]->AddEntry((TObject *)0, Form("    mean=%.4f cm", vtx_EW_mean[ifile][ivtx]), "");
        }
    }


    TLegend *legvtx[3];
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        if (ivtx == 2)
            legvtx[ivtx] = new TLegend(0.3, 0.25, 0.7, 0.5);
        else
            legvtx[ivtx] = new TLegend(0.53, 0.75, 1., 0.9);
        legvtx[ivtx]->SetFillStyle(0);
        legvtx[ivtx]->SetBorderSize(0);
        legvtx[ivtx]->SetTextSize(0.03);
        for (int ifile = 0; ifile < NFILES; ifile++)
        {
            legvtx[ivtx]->AddEntry(hvtx[ifile][ivtx], fileLabel[ifile], "L");
            legvtx[ivtx]->AddEntry((TObject *)0, Form("    mean=%.4f cm", hvtx[ifile][ivtx]->GetMean()), "");
        }
    }


    TLatex ltitle;
    ltitle.SetNDC();
    ltitle.SetTextAlign(22);

    TLatex llabel;
    llabel.SetNDC();
    llabel.SetTextAlign(22);

    //================================================//
    // PLOT
    //================================================//

    //--- SvxCentralTracks ---//
    TCanvas *cdca_pt_svxcnt[NARM];
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdca_pt_svxcnt[iarm] = new TCanvas(Form("cdca_pt_svxcnt_%i", iarm),
                                           Form("svxcnt dca v pt arm%i", iarm),
                                           1200, 700);
        cdca_pt_svxcnt[iarm]->Divide(3, 2);
        for (int ipt = 0; ipt < NPTSVXCNT; ipt++)
        {
            cdca_pt_svxcnt[iarm]->GetPad(ipt + 1)->SetRightMargin(0.02);
            cdca_pt_svxcnt[iarm]->GetPad(ipt + 1)->SetTicks(1, 1);

            cdca_pt_svxcnt[iarm]->cd(ipt + 1);
            //gPad->SetLogy();
            hdca_pt_svxcnt[iarm][0][ipt]->GetXaxis()->SetRangeUser(-0.06, 0.06);
            hdca_pt_svxcnt[iarm][0][ipt]->Draw("HIST");
            for (int ifile = 1; ifile < NFILES; ifile++)
                hdca_pt_svxcnt[iarm][ifile][ipt]->Draw("HIST same");
            leg->Draw("same");

        }
    }


    TCanvas *cdca_svxcnt_fits[NARM];
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdca_svxcnt_fits[iarm] = new TCanvas(Form("cdca_svxcnt_fits_%i", iarm),
                                             Form("svxcnt dca fits arm%i", iarm),
                                             1200, 700);
        cdca_svxcnt_fits[iarm]->Divide(NPTSVXCNT, NFILES);
        TF1 *ftmp;
        for (int ifile = 0; ifile < NFILES; ifile++)
        {
            for (int ipt = 0; ipt < NPTSVXCNT; ipt++)
            {
                cdca_svxcnt_fits[iarm]->cd(ifile * NPTSVXCNT + ipt + 1);
                gPad->SetLogy();
                hdca_pt_svxcnt[iarm][ifile][ipt]->Draw();
                ftmp = (TF1 *) hdca_pt_svxcnt[iarm][ifile][ipt]->GetFunction("fgaus");
                if (!ftmp)
                {
                    cout << "WARNING!! No fit function for ifile=" << ifile << " ipt=" << ipt << endl;
                    continue;
                }

                ftmp->DrawCopy("same");
            }
        }
    }



    TCanvas *cdcares_svxcnt = new TCanvas("cdcares_svxcnt",
                                          "dca res svxcnt",
                                          1200, 600);
    cdcares_svxcnt->Divide(NARM, 1, 0, 0);
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdcares_svxcnt->GetPad(iarm + 1)->SetTopMargin(0.08);
        if (iarm == 0) cdcares_svxcnt->GetPad(iarm + 1)->SetRightMargin(0.00);
        else cdcares_svxcnt->GetPad(iarm + 1)->SetRightMargin(0.06);
        cdcares_svxcnt->GetPad(iarm + 1)->SetBottomMargin(0.12);
        if (iarm == 1) cdcares_svxcnt->GetPad(iarm + 1)->SetLeftMargin(0.0);
        else cdcares_svxcnt->GetPad(iarm + 1)->SetLeftMargin(0.12);
        cdcares_svxcnt->GetPad(iarm + 1)->SetTicks(1, 1);

        cdcares_svxcnt->cd(iarm + 1);
        gdca_ptres_svxcnt[iarm][0]->GetYaxis()->SetRangeUser(0, 175);
        gdca_ptres_svxcnt[iarm][0]->Draw("AP");
        for (int ifile = 1; ifile < NFILES; ifile++)
            gdca_ptres_svxcnt[iarm][ifile]->Draw("P");
        leg->Draw("same");

        ltitle.DrawLatex(0.5, 0.95, Form("DC Arm=%i (%s)", iarm, armLabel[iarm]));
    }




    TCanvas *cdcaphi_svxcnt = new TCanvas("cdcaphi_svxcnt", "svxcnt dca v phi", 1200, 700);
    cdcaphi_svxcnt->Divide(NFILES, 1);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        cdcaphi_svxcnt->cd(ifile + 1);
        hdcaphi_svxcnt[ifile]->Draw("colz");
        pdcaphi_svxcnt[ifile]->Draw("P same");
        ltitle.DrawLatex(0.5, 0.95, fileLabel[ifile]);
        llabel.SetTextColor(kGray);
        llabel.DrawLatex(0.25, 0.05, "WEST");
        llabel.DrawLatex(0.75, 0.05, "EAST");

    }



    TCanvas *cdcaz_pt_svxcnt[NARM];
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdcaz_pt_svxcnt[iarm] = new TCanvas(Form("cdcaz_pt_svxcnt_%i", iarm),
                                            Form("svxcnt dcaz v pt arm%i", iarm),
                                            1200, 700);
        cdcaz_pt_svxcnt[iarm]->Divide(3, 2);
        for (int ipt = 0; ipt < NPTSVXCNT; ipt++)
        {
            cdcaz_pt_svxcnt[iarm]->GetPad(ipt + 1)->SetRightMargin(0.02);
            cdcaz_pt_svxcnt[iarm]->GetPad(ipt + 1)->SetTicks(1, 1);

            cdcaz_pt_svxcnt[iarm]->cd(ipt + 1);
            //gPad->SetLogy();
            hdcaz_pt_svxcnt[iarm][0][ipt]->GetXaxis()->SetRangeUser(-0.2, 0.2);
            hdcaz_pt_svxcnt[iarm][0][ipt]->Draw();
            for (int ifile = 1; ifile < NFILES; ifile++)
                hdcaz_pt_svxcnt[iarm][ifile][ipt]->Draw("same");
            leg->Draw("same");

        }
    }


    TCanvas *cdcazphi_svxcnt = new TCanvas("cdcazphi_svxcnt", "svxcnt dcaz v phi", 1200, 700);
    cdcazphi_svxcnt->Divide(NFILES, 1);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        cdcazphi_svxcnt->cd(ifile + 1);
        hdcazphi_svxcnt[ifile]->Draw("colz");
        pdcazphi_svxcnt[ifile]->Draw("P same");
        ltitle.DrawLatex(0.5, 0.95, fileLabel[ifile]);
        llabel.SetTextColor(kGray);
        llabel.DrawLatex(0.25, 0.05, "WEST");
        llabel.DrawLatex(0.75, 0.05, "EAST");
    }



    TCanvas *cntracksphi_svxcnt = new TCanvas("cntracksphi_svxcnt", "n tracks v phi svxcnt", 900, 600);
    cntracksphi_svxcnt->cd(1);
    hntracksphi_svxcnt[0]->Draw("hist");
    for (int ifile = 1; ifile < NFILES; ifile++)
        hntracksphi_svxcnt[ifile]->Draw("hist same");
    legwevents->Draw("same");
    llabel.SetTextColor(kGray);
    llabel.DrawLatex(0.25, 0.05, "WEST");
    llabel.DrawLatex(0.75, 0.05, "EAST");



    TCanvas *cntracksphi_eff = new TCanvas("cntracksphi_eff", "n tracks v phi eff", 900, 600);
    cntracksphi_eff->cd(1);
    hntracksphi_eff[0]->GetYaxis()->SetRangeUser(0.0, 0.4);
    hntracksphi_eff[0]->Draw("hist E");
    for (int ifile = 1; ifile < NFILES; ifile++)
        hntracksphi_eff[ifile]->Draw("hist E same");
    legwevents->Draw("same");
    llabel.SetTextColor(kGray);
    llabel.DrawLatex(0.25, 0.03, "WEST");
    llabel.DrawLatex(0.75, 0.03, "EAST");



    TCanvas *cchisq_svxcnt = new TCanvas("cchisq_svxcnt", "svxcnt chisq/ndf", 1200, 600);
    cchisq_svxcnt->Divide(NARM, 1);
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cchisq_svxcnt->cd(iarm + 1);
        hchisqndf_svxcnt[iarm][0]->Draw();
        for (int ifile = 1; ifile < NFILES; ifile++)
            hchisqndf_svxcnt[iarm][ifile]->Draw("same");
        if (iarm == 0) leg->Draw("same");

        ltitle.DrawLatex(0.5, 0.95, Form("DC Arm=%i (%s)", iarm, armLabel[iarm]));
    }


    TCanvas *cresz[NLAYERS];
    float rangez[4] = {0.2, 0.3, 0.4, 0.5};
    int cx[4] = {5, 5, 4, 6};
    int cy[4] = {2, 4, 4, 4};
    for (int ilayer = 0; ilayer < NLAYERS; ilayer++)
    {
        cresz[ilayer] = new TCanvas(Form("cresz_%i", ilayer),
                                    Form("res z B%i", ilayer),
                                    cx[ilayer] * 200,
                                    cy[ilayer] * 150);
        cresz[ilayer]->Divide(cx[ilayer], cy[ilayer]);
        for (int iladder = 0; iladder < NLADDERS[ilayer]; iladder++)
        {
            cresz[ilayer]->cd(iladder + 1);
            if (hres_z[0][ilayer][iladder]->GetEntries() < 100) continue;
            hres_z[0][ilayer][iladder]->GetXaxis()->SetRangeUser(-1 * rangez[ilayer], rangez[ilayer]);
            hres_z[0][ilayer][iladder]->Draw("HIST");
            for (int ifile = 1; ifile < NFILES; ifile++)
                hres_z[ifile][ilayer][iladder]->Draw("hist same");
            //ltitle.DrawLatex(0.5, 0.95, Form("B%iL%02i", ilayer, iladder));
        }
    }

    TCanvas *cress[NLAYERS];
    float ranges[4] = {0.2, 0.4, 2.0, 4.0};
    for (int ilayer = 0; ilayer < NLAYERS; ilayer++)
    {
        cress[ilayer] = new TCanvas(Form("cress_%i", ilayer),
                                    Form("res s B%i", ilayer),
                                    cx[ilayer] * 200,
                                    cy[ilayer] * 150);
        cress[ilayer]->Divide(cx[ilayer], cy[ilayer]);
        for (int iladder = 0; iladder < NLADDERS[ilayer]; iladder++)
        {
            cress[ilayer]->cd(iladder + 1);
            if (hres_z[0][ilayer][iladder]->GetEntries() < 100) continue;
            hres_s[0][ilayer][iladder]->GetXaxis()->SetRangeUser(-1 * ranges[ilayer], ranges[ilayer]);
            hres_s[0][ilayer][iladder]->Draw("hist");
            for (int ifile = 1; ifile < NFILES; ifile++)
                hres_s[ifile][ilayer][iladder]->Draw("hist same");
            //ltitle.DrawLatex(0.5, 0.95, Form("B%iL%02i", ilayer, iladder));
        }
    }

    TCanvas *creszsum = new TCanvas("creszsum", "z res summary", 1200, 600);
    creszsum->cd(1);
    gres_z[0]->GetYaxis()->SetRangeUser(-0.5, 0.5);
    gres_z[0]->Draw("AP");
    for (int ifile = 1; ifile < NFILES; ifile++)
        gres_z[ifile]->Draw("P");


    TCanvas *cresssum = new TCanvas("cresssum", "s res summary", 1200, 600);
    cresssum->cd(1);
    gres_s[0]->GetYaxis()->SetRangeUser(-3., 3.);
    gres_s[0]->Draw("AP");
    for (int ifile = 1; ifile < NFILES; ifile++)
        gres_s[ifile]->Draw("P");


    TCanvas *cnclus_svxcnt = new TCanvas("cnclus_svxcnt", "nclus svxcnt", 1200, 600);
    cnclus_svxcnt->Divide(NARM, 1);
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cnclus_svxcnt->cd(iarm + 1);
        hNclus_svxcnt[iarm][0]->Draw();
        for (int ifile = 1; ifile < NFILES; ifile++)
            hNclus_svxcnt[iarm][ifile]->Draw("same");
        leg->Draw("same");
    }



    //--- Standalone ---//
    TCanvas *cdca_pt_standalone[NARM];
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdca_pt_standalone[iarm] = new TCanvas(Form("cdca_pt_standalone_%i", iarm),
                                               Form("standalone dca v pt arm%i", iarm),
                                               1200, 700);
        cdca_pt_standalone[iarm]->Divide(3, 2);
        for (int ipt = 0; ipt < NPTSVXCNT; ipt++)
        {
            cdca_pt_standalone[iarm]->GetPad(ipt + 1)->SetRightMargin(0.02);
            cdca_pt_standalone[iarm]->GetPad(ipt + 1)->SetTicks(1, 1);

            cdca_pt_standalone[iarm]->cd(ipt + 1);
            //gPad->SetLogy();
            hdca_pt_standalone[iarm][0][ipt]->Draw();
            for (int ifile = 1; ifile < NFILES; ifile++)
                hdca_pt_standalone[iarm][ifile][ipt]->Draw("same");
            leg->Draw("same");

        }
    }

    TCanvas *cdca_standalone = new TCanvas("cdca_standalone", "standalone dca", 600, 600);
    cdca_standalone->cd(1);
    gPad->SetLogy();
    hdca_standalone[0]->Draw();
    for (int ifile = 1; ifile < NFILES; ifile++)
        hdca_standalone[ifile]->Draw("same");
    leg->Draw("same");


    TCanvas *cdcaphi_standalone = new TCanvas("cdcaphi_standalone", "standalone dca v phi", 1200, 700);
    cdcaphi_standalone->Divide(NFILES, 1);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        cdcaphi_standalone->cd(ifile + 1);
        hdcaphi_standalone[ifile]->Draw("colz");
        pdcaphi_standalone[ifile]->Draw("P same");
        ltitle.DrawLatex(0.5, 0.95, fileLabel[ifile]);
    }


    TCanvas *cntracksphi_standalone = new TCanvas("cntracksphi_standalone", "n tracks v phi standalone", 900, 600);
    cntracksphi_standalone->cd(1);
    hntracksphi_standalone[0]->Draw("hist");
    for (int ifile = 1; ifile < NFILES; ifile++)
        hntracksphi_standalone[ifile]->Draw("hist same");
    legwevents->Draw("same");


    TCanvas *cchisq_standalone = new TCanvas("cchisq_standalone", "seg chisq/ndf", 1200, 600);
    cchisq_standalone->Divide(NARM, 1);
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cchisq_standalone->cd(iarm + 1);
        hchisqndf_standalone[iarm][0]->Draw();
        for (int ifile = 1; ifile < NFILES; ifile++)
            hchisqndf_standalone[iarm][ifile]->Draw("same");
        if (iarm == 0) leg->Draw("same");

        ltitle.DrawLatex(0.5, 0.95, Form("DC Arm=%i (%s)", iarm, armLabel[iarm]));
    }


    TCanvas *cnclus_standalone = new TCanvas("cnclus_standalone", "nclus svxcnt", 1200, 600);
    cnclus_standalone->Divide(NARM, 1);
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cnclus_standalone->cd(iarm + 1);
        hNclus_standalone[iarm][0]->Draw();
        for (int ifile = 1; ifile < NFILES; ifile++)
            hNclus_standalone[iarm][ifile]->Draw("same");
        leg->Draw("same");
    }



    //--- Event ---//
    TCanvas *cvtx = new TCanvas("cvtx", "vtx", 1200, 600);
    cvtx->Divide(3, 1);
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        cvtx->GetPad(ivtx + 1)->SetTopMargin(0.02);
        cvtx->GetPad(ivtx + 1)->SetRightMargin(0.02);
        cvtx->GetPad(ivtx + 1)->SetBottomMargin(0.12);
        cvtx->GetPad(ivtx + 1)->SetLeftMargin(0.12);
        cvtx->GetPad(ivtx + 1)->SetTicks(1, 1);

        cvtx->cd(ivtx + 1);
        hvtx[0][ivtx]->Draw();
        for (int ifile = 0; ifile < NFILES; ifile++)
            hvtx[ifile][ivtx]->Draw("same");

        legvtx[ivtx]->Draw("same");
    }

    TCanvas *cvtxbbc = new TCanvas("cvtxbbc", "vtx-bbc", 600, 600);
    cvtxbbc->cd(1);
    hvtxzbbcz[0]->Draw();
    for (int ifile = 1; ifile < NFILES; ifile++)
        hvtxzbbcz[ifile]->Draw("same");
    leg->Draw("same");

    TCanvas *cvtxewdiff = new TCanvas("cvtxewdiff", "W_E vtx", 1200, 500);
    cvtxewdiff->Divide(3, 1);
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        cvtxewdiff->GetPad(ivtx + 1)->SetTopMargin(0.02);
        cvtxewdiff->GetPad(ivtx + 1)->SetTicks(1, 1);

        cvtxewdiff->cd(ivtx + 1);
        hvtx_EW[0][ivtx]->Draw();
        for (int ifile = 0; ifile < NFILES; ifile++)
            hvtx_EW[ifile][ivtx]->Draw("same");

        legvtxEW[ivtx]->Draw("same");
    }

    TCanvas *cvtxew = new TCanvas("cvtxew", "W_E vtx", 1200, 700);
    cvtxew->Divide(3, 2);
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        //E
        cvtxew->GetPad(ivtx + 1)->SetTopMargin(0.02);
        cvtxew->GetPad(ivtx + 1)->SetTicks(1, 1);

        cvtxew->cd(ivtx + 1);
        hvtx_E[0][ivtx]->Draw();
        for (int ifile = 0; ifile < NFILES; ifile++)
            hvtx_E[ifile][ivtx]->Draw("same");

        if (ivtx == 0) leg->Draw("same");

        //W
        cvtxew->GetPad(3 + ivtx + 1)->SetTopMargin(0.02);
        cvtxew->GetPad(3 + ivtx + 1)->SetTicks(1, 1);

        cvtxew->cd(3 + ivtx + 1);
        hvtx_W[0][ivtx]->Draw();
        for (int ifile = 0; ifile < NFILES; ifile++)
            hvtx_W[ifile][ivtx]->Draw("same");

        if (ivtx == 0) leg->Draw("same");

    }

    //================================================//
    // PRINT PLOTS
    //================================================//
    if (print_plots)
    {
        cout << endl;
        cout << "--> Printing plots" << endl;

        //-- Event
        cvtx->Print("pdfs/prim_vertex.pdf");
        cvtxew->Print("pdfs/prim_vertex_ew.pdf");
        cvtxewdiff->Print("pdfs/ew_offset.pdf");
        cvtxbbc->Print("pdfs/vtxz-bbcz.pdf");

        //-- Standalone
        cdca_standalone->Print("pdfs/dca2d_standalone.pdf");
        cntracksphi_standalone->Print("pdfs/nstandalone_phi.pdf");
        cdcaphi_standalone->Print("pdfs/dca2d_vs_phi_standalone.pdf");
        cchisq_standalone->Print("pdfs/chisqndf_standalone.pdf");
        for (int iarm = 0; iarm < NARM; iarm++)
            cdca_pt_standalone[iarm]->Print(Form("pdfs/dca2d_standalone_vs_pt_dcarm%i.pdf", iarm));
        cnclus_standalone->Print("pdfs/nclus_standalone.pdf");

        //-- SvxCentral
        cntracksphi_svxcnt->Print("pdfs/nsvxcnt_phi.pdf");
        cntracksphi_eff->Print("pdfs/nsvxcnteff_phi.pdf");

        cdcaphi_svxcnt->Print("pdfs/dca2d_vs_phi.pdf");
        cdcares_svxcnt->Print("pdfs/dca2d_res_pt.pdf");
        for (int iarm = 0; iarm < NARM; iarm++)
            cdca_pt_svxcnt[iarm]->Print(Form("pdfs/dca2d_vs_pt_dcarm%i.pdf", iarm));

        cdcazphi_svxcnt->Print("pdfs/dcaz_vs_phi.pdf");
        for (int iarm = 0; iarm < NARM; iarm++)
            cdcaz_pt_svxcnt[iarm]->Print(Form("pdfs/dcaz_vs_pt_dcarm%i.pdf", iarm));

        cchisq_svxcnt->Print("pdfs/chisqndf_svxcnt.pdf");

        for (int ilayer = 0; ilayer < NLAYERS; ilayer++)
        {
            cresz[ilayer]->Print(Form("pdfs/svxcnt_res_z_B%i.pdf", ilayer));
            cress[ilayer]->Print(Form("pdfs/svxcnt_res_s_B%i.pdf", ilayer));
        }
        creszsum->Print("pdfs/svxcnt_res_z_summary.pdf");
        cresssum->Print("pdfs/svxcnt_res_s_summary.pdf");
        cnclus_svxcnt->Print("pdfs/nclus_svxcnt.pdf");
    }


}
