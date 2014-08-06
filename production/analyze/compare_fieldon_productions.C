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
    gStyle->SetPalette(56, 0, 0.5);

    //================================================//
    // SET RUNNING CONDITIONS
    //================================================//

    bool print_plots = true;

    const int NFILES = 5;
    const char *fileName[] =
    {
        //        "/direct/phenix+prod01/phnxreco/millepede/fieldon/taebong-407951-p2-v8/ntuple/testvtxproduction_taebong-407951-p2-v8.root",
        //        "/direct/phenix+prod01/phnxreco/millepede/fieldon/ideal-407951-2-3/ntuple/testvtxproduction_407951-2-3.root",
        "testvtxproduction_taebong-407951-p2-v8.root",
        //"testvtxproduction_407951-0-0.root",
        //"testvtxproduction_407951-2-3.root",
        //"testvtxproduction_407951-2-3-1.root",
        //"testvtxproduction_407951-2-3-2.root",
        "testvtxproduction_407951-2-3-3.root",
        "testvtxproduction_407951-2-3-4.root",
        "testvtxproduction_407951-2-3-5.root",
        "testvtxproduction_407951-2-3-6.root",
    };
    const char *fileLabel[] =
    {
        "Taebong p2 v8",
        //"Ideal Geo",
        //"Andrew 2-3",
        //"Andrew 2-3-1 (Z shift)",
        //"Andrew 2-3-2 (0 offsets)",
        "Andrew 2-3-3 (offsets)",
        "Andrew 2-3-4 (no cnt-z)",
        "Andrew 2-3-5 (add W-E)",
        "Andrew 2-3-6 (ideal geo)",
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

    const int NARM = 2; //separate E/W arms
    const char *armLabel[] = {"East", "West"};

    const char *trackCuts = "(dchquality == 31 || dchquality == 63) && "
                            "Nclus > 2 && Nclus_layer[0]>0 && Nclus_layer[1]>0 && "
                            "chisq/ndf < 3 && "
                            //"charge == 1 &&"
                            "pT>1";
    const char *eventCuts = "!tickCut && "
                            "TMath::Abs(vtx[2])<10 && "
                            "which_vtx==\"SVX_PRECISE\"";

    cout << " trackCuts = " << trackCuts << endl;
    cout << " eventCuts = " << eventCuts << endl;

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

    TH1F *hdca_svxcnt[NARM][NFILES];
    TH1F *hdca_standalone[NFILES];

    TH2F *htracks_phized[NFILES];

    TH1F *hntracksphi_svxcnt[NFILES];
    TH1F *hntracksphi_cnt[NFILES];
    TH1F *hntracksphi_eff[NFILES];
    //TH1F *hntracksphi_standalone[NFILES];
    float nevents[NFILES] = {0};

    TH1F *hvtx[NFILES][3];

    TH1F *hvtx_EW[NFILES][3];


    TH1F *hdca_pt_svxcnt[NARM][NFILES][NPTSVXCNT];
    TH1F *hdcaz_pt_svxcnt[NARM][NFILES][NPTSVXCNT];

    TGraphErrors *gdca_ptres_svxcnt[NARM][NFILES];

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

        ntp_SVXCNT->Draw("pT:dca2D*charge>>htmp(100,-0.1,0.1, 50,0,5)",
                         Form("%s && %s && dcarm==0", eventCuts, trackCuts),
                         "goff");
        hdcapt_svxcnt[0][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcapt_svxcnt[0][ifile]->SetDirectory(0);
        hdcapt_svxcnt[0][ifile]->SetName(Form("hdcapt_svxcnt_0_%i", ifile));
        hdcapt_svxcnt[0][ifile]->SetTitle(";DCA 2D * chg [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("pT:dca2D*charge>>htmp(100,-0.1,0.1, 50,0,5)",
                         Form("%s && %s && dcarm==1", eventCuts, trackCuts),
                         "goff");
        hdcapt_svxcnt[1][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcapt_svxcnt[1][ifile]->SetDirectory(0);
        hdcapt_svxcnt[1][ifile]->SetName(Form("hdcapt_svxcnt_1_%i", ifile));
        hdcapt_svxcnt[1][ifile]->SetTitle(";DCA 2D * chg [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("dca2D*charge:phi>>htmp(250,-1,4, 400,-0.06,0.06)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hdcaphi_svxcnt[ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcaphi_svxcnt[ifile]->SetDirectory(0);
        hdcaphi_svxcnt[ifile]->SetName(Form("hdcaphi_svxcnt_%i", ifile));
        hdcaphi_svxcnt[ifile]->SetTitle(";phi;DCA 2D * chg [cm]");
        hdcaphi_svxcnt[ifile]->GetYaxis()->SetTitleOffset(1.5);

        pdcaphi_svxcnt[ifile] = new TProfile(Form("pdcaphi_svxcnt_%i", ifile), "", 125, -1, 4);
        pdcaphi_svxcnt[ifile]->SetLineColor(kBlue);
        ntp_SVXCNT->Draw(Form("dca2D*charge:phi>>pdcaphi_svxcnt_%i", ifile),
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        pdcaphi_svxcnt[ifile]->SetDirectory(0);




        ntp_SVXCNT->Draw("pT:dcaz*charge>>htmp(700,-0.7,0.7, 50,0,5)",
                         Form("%s && %s && dcarm==0", eventCuts, trackCuts),
                         "goff");
        hdcazpt_svxcnt[0][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcazpt_svxcnt[0][ifile]->SetDirectory(0);
        hdcazpt_svxcnt[0][ifile]->SetName(Form("hdcazpt_svxcnt_0_%i", ifile));
        hdcazpt_svxcnt[0][ifile]->SetTitle(";DCA Z * chg [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("pT:dcaz*charge>>htmp(700,-0.7,0.7, 50,0,5)",
                         Form("%s && %s && dcarm==1", eventCuts, trackCuts),
                         "goff");
        hdcazpt_svxcnt[1][ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcazpt_svxcnt[1][ifile]->SetDirectory(0);
        hdcazpt_svxcnt[1][ifile]->SetName(Form("hdcazpt_svxcnt_1_%i", ifile));
        hdcazpt_svxcnt[1][ifile]->SetTitle(";DCA Z * chg [cm]; p_{T} [GeV/c]");

        ntp_SVXCNT->Draw("dcaz*charge:phi>>htmp(250,-1,4, 200,-0.2,0.2)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hdcazphi_svxcnt[ifile] = (TH2F *) gDirectory->FindObject("htmp");
        hdcazphi_svxcnt[ifile]->SetDirectory(0);
        hdcazphi_svxcnt[ifile]->SetName(Form("hdcazphi_svxcnt_%i", ifile));
        hdcazphi_svxcnt[ifile]->SetTitle(";phi;DCA Z * chg [cm]");
        hdcazphi_svxcnt[ifile]->GetYaxis()->SetTitleOffset(1.5);

        pdcazphi_svxcnt[ifile] = new TProfile(Form("pdcazphi_svxcnt_%i", ifile), "", 125, -1, 4);
        pdcazphi_svxcnt[ifile]->SetLineColor(kBlue);
        ntp_SVXCNT->Draw(Form("dcaz*charge:phi>>pdcazphi_svxcnt_%i", ifile),
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        pdcazphi_svxcnt[ifile]->SetDirectory(0);




        ntp_SVXCNT->Draw("phi>>htmp(125,-1,4)",
                         Form("%s && %s", eventCuts, trackCuts),
                         "goff");
        hntracksphi_svxcnt[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hntracksphi_svxcnt[ifile]->SetDirectory(0);
        hntracksphi_svxcnt[ifile]->SetName(Form("hntracksphi_svxcnt_%i", ifile));
        hntracksphi_svxcnt[ifile]->SetTitle(";phi;N SvxCentralTracks / event");
        hntracksphi_svxcnt[ifile]->GetYaxis()->SetTitleOffset(1.5);
        hntracksphi_svxcnt[ifile]->SetLineColor(color[ifile]);
        hntracksphi_svxcnt[ifile]->Sumw2();



        //--- CNT ntuple ---//
        ntp_CNT = (TTree *) fin->Get("ntp_CNT");
        if (!ntp_CNT)
        {
            cout << "ERROR!! Unable to find ntp_CNT in " << fileName[ifile] << endl;
            return;
        }

        ntp_SVXCNT->Draw("phi>>htmp(125,-1,4)",
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

        //ntp_SEG->Draw("dca2D*charge>>htmp(700,-0.7,0.7)", "pT>0.5 && pT<1. && chisq/ndf<3 && Nclus==4", "goff");
        ntp_SEG->Draw("dca2D*charge>>htmp(350,-0.7,0.7)", "Nclus_layer[0]>0 && Nclus_layer[1]>0 && pT>0.5", "goff");
        //ntp_SEG->Draw("dca2D*charge>>htmp(350,-0.7,0.7)", "", "goff", 1000000);
        hdca_standalone[ifile] = (TH1F *) gDirectory->FindObject("htmp");
        hdca_standalone[ifile]->SetDirectory(0);
        hdca_standalone[ifile]->SetName(Form("hdca_standalone_%i", ifile));
        hdca_standalone[ifile]->SetTitle(";DCA 2D * chg [cm]");
        hdca_standalone[ifile]->SetLineColor(color[ifile]);
        hdca_standalone[ifile]->SetLineWidth(2);
        //hdca_standalone[ifile]->Scale(1. / hdca_standalone[ifile]->Integral(1, hdca_standalone[ifile]->GetNbinsX()));


        //histograms
        htracks_phized[ifile] = (TH2F *) fin->Get("htrack_phized");
        if (!htracks_phized[ifile])
        {
            cout << "ERROR!! Unable to find htrack_phized in " << fileName[ifile] << endl;
            return;
        }
        htracks_phized[ifile]->SetName(Form("htracks_phized_%i", ifile));
        htracks_phized[ifile]->SetDirectory(0);
        htracks_phized[ifile]->Scale(1. / htracks_phized[ifile]->Integral());



        //--- Event ntuple ---//
        ntp_event = (TTree *) fin->Get("ntp_event");
        if (!ntp_event)
        {
            cout << "ERROR!! Unable to get ntp_event from " << fileName[ifile] << endl;
            return;
        }

        ntp_event->Draw("run>>htmp", eventCuts, "goff");
        nevents[ifile] = ((TH1F *) gDirectory->FindObject("htmp"))->GetEntries();
        cout << "  Found " << nevents[ifile] << " events in " << fileName[ifile] << endl;

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
                hdca_pt_svxcnt[iarm][ifile][ipt] = (TH1F *) hdcapt_svxcnt[iarm][ifile]->ProjectionX(Form("hdca_pt_svxcnt_%i_%i_%i", iarm, ifile, ipt), bl, bh);
                hdca_pt_svxcnt[iarm][ifile][ipt]->SetTitle(Form("DCArm%i %.1f<p_{T}<%.1f;DCA 2D * chg [cm]", iarm, ptl_svxcnt[ipt], pth_svxcnt[ipt]));
                hdca_pt_svxcnt[iarm][ifile][ipt]->SetLineColor(color[ifile]);
                //hdca_pt_svxcnt[iarm][ifile][ipt]->Scale(1. / hdca_pt_svxcnt[iarm][ifile]->Integral(1, hdca_pt_svxcnt[iarm][ifile]->GetNbinsX()));

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
                hdcaz_pt_svxcnt[iarm][ifile][ipt] = (TH1F *) hdcazpt_svxcnt[iarm][ifile]->ProjectionX(Form("hdcaz_pt_svxcnt_%i_%i_%i", iarm, ifile, ipt), bl, bh);
                hdcaz_pt_svxcnt[iarm][ifile][ipt]->SetTitle(Form("DC Arm=%i %.1f<p_{T}<%.1f;DCA Z * chg [cm]", iarm, ptl_svxcnt[ipt], pth_svxcnt[ipt]));
                hdcaz_pt_svxcnt[iarm][ifile][ipt]->SetLineColor(color[ifile]);
                //hdcaz_pt_svxcnt[iarm][ifile][ipt]->Scale(1. / hdcaz_pt_svxcnt[iarm][ifile]->Integral(1, hdcaz_pt_svxcnt[iarm][ifile]->GetNbinsX()));



            }
        }

    }

    //================================================//
    // PLOT OBJECTS
    //================================================//
    cout << endl;
    cout << "--> Plotting" << endl;


    TLegend *leg = new TLegend(0.55, 0.75, 0.9, 0.9);
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

    TLegend *legvtx[3];
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        legvtx[ivtx] = new TLegend(0.5, 0.75, 1., 0.9);
        legvtx[ivtx]->SetFillStyle(0);
        legvtx[ivtx]->SetBorderSize(0);
        legvtx[ivtx]->SetTextSize(0.03);
        for (int ifile = 0; ifile < NFILES; ifile++)
        {
            legvtx[ivtx]->AddEntry(hvtx_EW[ifile][ivtx], fileLabel[ifile], "L");
            legvtx[ivtx]->AddEntry((TObject *)0, Form("    mean=%.4f cm", vtx_EW_mean[ifile][ivtx]), "");
        }
    }


    TLatex ltitle;
    ltitle.SetNDC();
    ltitle.SetTextAlign(22);

    //================================================//
    // PLOT
    //================================================//

    //--- SvxCentralTracks ---//
    /*
    TCanvas *cdca_svxcnt = new TCanvas("cdca_svxcnt", "svxcnt dca", 600, 600);
    cdca_svxcnt->Divide(2, 1);

    cdca_svxcnt->cd(1);
    gPad->SetLogy();
    hdca_svxcnt[0][0]->Draw();
    for (int ifile = 1; ifile < NFILES; ifile++)
        hdca_svxcnt[0][ifile]->Draw("same");
    leg->Draw("same");

    cdca_svxcnt->cd(2);
    gPad->SetLogy();
    hdca_svxcnt[1][0]->Draw();
    for (int ifile = 1; ifile < NFILES; ifile++)
        hdca_svxcnt[1][ifile]->Draw("same");
    leg->Draw("same");
    */


    TCanvas *cdca_pt_svxcnt[NARM];
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdca_pt_svxcnt[iarm] = new TCanvas(Form("cdca_pt_svxcnt_%i", iarm),
                                           Form("svxcnt dca v pt arm%i", iarm),
                                           1200, 700);
        cdca_pt_svxcnt[iarm]->Divide(3, 2);
        for (int ipt = 0; ipt < NPTSVXCNT; ipt++)
        {
            cdca_pt_svxcnt[iarm]->cd(ipt + 1);
            gPad->SetLogy();
            hdca_pt_svxcnt[iarm][0][ipt]->Draw();
            for (int ifile = 1; ifile < NFILES; ifile++)
                hdca_pt_svxcnt[iarm][ifile][ipt]->Draw("same");
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
    cdcares_svxcnt->Divide(NARM, 1);
    for (int iarm = 0; iarm < NARM; iarm++)
    {
        cdcares_svxcnt->cd(iarm + 1);
        gdca_ptres_svxcnt[iarm][0]->GetYaxis()->SetRangeUser(20, 300);
        gdca_ptres_svxcnt[iarm][0]->Draw("AP");
        for (int ifile = 1; ifile < NFILES; ifile++)
            gdca_ptres_svxcnt[iarm][ifile]->Draw("P");
        leg->Draw("same");

        ltitle.DrawLatex(0.5, 0.95, Form("DC Arm=%i (%s?)", iarm, armLabel[iarm]));
    }




    TCanvas *cdcaphi_svxcnt = new TCanvas("cdcaphi_svxcnt", "svxcnt dca v phi", 1200, 700);
    cdcaphi_svxcnt->Divide(NFILES, 1);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        cdcaphi_svxcnt->cd(ifile + 1);
        hdcaphi_svxcnt[ifile]->Draw("colz");
        pdcaphi_svxcnt[ifile]->Draw("P same");
        ltitle.DrawLatex(0.5, 0.95, fileLabel[ifile]);
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
            cdcaz_pt_svxcnt[iarm]->cd(ipt + 1);
            gPad->SetLogy();
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
    }



    /*
    TCanvas *ctracks_phized = new TCanvas("ctracks_phized", "tracks phi vs zed svxcnt", 1200, 700);
    ctracks_phized->Divide(NFILES, 1);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        ctracks_phized->cd(ifile + 1);
        htracks_phized[ifile]->Draw("colz");
    }
    */


    TCanvas *cntracksphi_svxcnt = new TCanvas("cntracksphi_svxcnt", "n tracks v phi svxcnt", 900, 600);
    cntracksphi_svxcnt->cd(1);
    hntracksphi_svxcnt[0]->Draw("hist");
    for (int ifile = 1; ifile < NFILES; ifile++)
        hntracksphi_svxcnt[ifile]->Draw("hist same");
    legwevents->Draw("same");


    TCanvas *cntracksphi_eff = new TCanvas("cntracksphi_eff", "n tracks v phi eff", 900, 600);
    cntracksphi_eff->cd(1);
    hntracksphi_eff[0]->GetYaxis()->SetRangeUser(0.0, 1.1);
    hntracksphi_eff[0]->Draw();
    for (int ifile = 1; ifile < NFILES; ifile++)
        hntracksphi_eff[ifile]->Draw("same");
    //legwevents->Draw("same");



    //--- Standalone ---//
    TCanvas *cdca_standalone = new TCanvas("cdca_standalone", "standalone dca", 600, 600);
    cdca_standalone->cd(1);
    gPad->SetLogy();
    hdca_standalone[0]->Draw();
    for (int ifile = 1; ifile < NFILES; ifile++)
        hdca_standalone[ifile]->Draw("same");
    leg->Draw("same");

    /*
    TCanvas *cdca_each = new TCanvas("cdca_each", "standalone each", 1200, 700);
    cdca_each->Divide(NFILES, 1);
    for (int ifile = 0; ifile < NFILES; ifile++)
    {
        cdca_each->cd(ifile + 1);
        gPad->SetLogy();
        hdca_standalone[ifile]->Draw();
    }
    */


    //--- Event ---//
    TCanvas *cvtx = new TCanvas("cvtx", "vtx", 1200, 600);
    cvtx->Divide(3, 1);
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        cvtx->cd(ivtx + 1);
        hvtx[0][ivtx]->Draw();
        for (int ifile = 0; ifile < NFILES; ifile++)
            hvtx[ifile][ivtx]->Draw("same");

        if (ivtx == 0) leg->Draw("same");
    }

    TCanvas *cvtxew = new TCanvas("cvtxew", "W_E vtx", 1200, 600);
    cvtxew->Divide(3, 1);
    for (int ivtx = 0; ivtx < 3; ivtx++)
    {
        cvtxew->cd(ivtx + 1);
        hvtx_EW[0][ivtx]->Draw();
        for (int ifile = 0; ifile < NFILES; ifile++)
            hvtx_EW[ifile][ivtx]->Draw("same");

        legvtx[ivtx]->Draw("same");
    }


    //================================================//
    // PRINT PLOTS
    //================================================//
    if (print_plots)
    {
        cout << endl;
        cout << "--> Printing plots" << endl;

        cvtx->Print("pdfs/prim_vertex.pdf");
        cvtxew->Print("pdfs/ew_offset.pdf");

        cdca_standalone->Print("pdfs/dca2d_standalone.pdf");

        cntracksphi_svxcnt->Print("pdfs/nsvxcnt_phi.pdf");
        cntracksphi_eff->Print("pdfs/nsvxcnteff_phi.pdf");

        cdcaphi_svxcnt->Print("pdfs/dca2d_vs_phi.pdf");
        cdcares_svxcnt->Print("pdfs/dca2d_res_pt.pdf");
        for (int iarm = 0; iarm < NARM; iarm++)
            cdca_pt_svxcnt[iarm]->Print(Form("pdfs/dca2d_vs_pt_dcarm%i.pdf", iarm));

        cdcazphi_svxcnt->Print("pdfs/dcaz_vs_phi.pdf");
        for (int iarm = 0; iarm < NARM; iarm++)
            cdcaz_pt_svxcnt[iarm]->Print(Form("pdfs/dcaz_vs_pt_dcarm%i.pdf", iarm));
    }


}
