////////////////////////////////////////////////////
//
// Calculate the (x,y,z) beam center from the
// DC using PHCentralTracks
//
////////////////////////////////////////////////////
//
// Darren McGlinchey
// 8-22-2014
//
////////////////////////////////////////////////////
//
// NOTES:
//
//
////////////////////////////////////////////////////


#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>

#include "GLSFitter.h"

#include <iostream>
#include <vector>

using namespace std;


double Median(const TH1D *h1)
{
    /*
    int n = h1->GetXaxis()->GetNbins();
    std::vector<double>  x(n);
    h1->GetXaxis()->GetCenter( &x[0] );
    const double *y = h1->GetArray();
    // exclude underflow/overflows from bin content array y
    return TMath::Median(n, &x[0], &y[1]);
    */
    int n = 0;
    std::vector<double> x;
    std::vector<double> w;
    for (int i = 1; i <= h1->GetNbinsX(); i++)
    {
        if (h1->GetBinContent(i) > 0)
        {
            x.push_back(h1->GetXaxis()->GetBinCenter(i));
            w.push_back(h1->GetBinContent(i));
            n++;
        }
    }
    return TMath::Median(n, &x[0], &w[0]);
}

void CalcDCBeamCenter()
{

    //==================================================//
    // SET RUNNING CONDITIONS
    //==================================================//

    //-- input file containing ntp_CNT
    const char *inFile = "testvtxproduction_0000411768-0001.root";

    //-- DC beam center [x, y]
    // run 411768 from Taebong
    float dcbeamcent[2] = { -0.101223, 0.0838782};

    //-- PC radii from phnx.par:
    //   pc1rinst = 246.3 cm,
    //   pc2rinst = 415.00 cm,
    //   pc3rinst = 488.05 cm,
    float pcr[3] = {246.3, 415.00, 488.05};

    //-- PC resolution from technical note [cm]
    //   PC1 z res 1.7mm
    //   PC2&3 maintain same angular resolution
    float pcres[3] =
    {
        (float)(1. / TMath::Sqrt(12) * 0.17), //PC1
        (float)(1. / TMath::Sqrt(12) * 0.29), //PC2
        (float)(1. / TMath::Sqrt(12) * 0.34) //PC3
    };

    //==================================================//
    // DECLARE VARIABLES
    //==================================================//

    //-- Tree information
    TTree *ntp_CNT;
    int run;
    int event;
    int dchindex;
    float dchquality;
    float the0;
    float phi0;
    int dcarm;
    float zed;
    float phi;
    float vtx[3];
    float alpha;
    float ppc1[3]; //track projection at pc1 [x,y,z]
    float ppc2[3]; //track projection at pc2 [x,y,z]
    float ppc3[3]; //track projection at pc3 [x,y,z]
    float pc2dphi;
    float pc2dz;
    float pc3dphi;
    float pc3dz;
    float bbcq;
    float bbcz;
    float n0;

    //--
    TH1F *hdcz = new TH1F("hdcz", "; DC z_{vtx} [cm]", 150, -15, 15);
    TH1F *hdcbbcz = new TH1F("hdcbbcz", "; DC - BBC z_{vtx} [cm]", 100, -15, 15);
    TH1F *hdcvtxz = new TH1F("hdcvtxz", "; DC - VTX z_{vtx} [cm]", 100, -15, 15);


    //==================================================//
    // GET TREE FROM FILE
    //==================================================//
    cout << endl;
    cout << "--> Getting Tree's from " << inFile << endl;

    TFile *fin = TFile::Open(inFile);
    if (!fin)
    {
        cout << "ERROR!! Unable to open " << inFile << endl;
        return;
    }

    ntp_CNT = (TTree *) fin->Get("ntp_CNT");
    if (!ntp_CNT)
    {
        cout << "ERROR!! Unable to find ntp_CNT in " << inFile << endl;
        return;
    }

    //-- Set Branches

    ntp_CNT->SetBranchAddress("run",        &run);
    ntp_CNT->SetBranchAddress("event",      &event);
    ntp_CNT->SetBranchAddress("bbcq",       &bbcq);
    ntp_CNT->SetBranchAddress("bbcz",       &bbcz);
    ntp_CNT->SetBranchAddress("vtx",        &vtx);
    ntp_CNT->SetBranchAddress("dchindex",   &dchindex);
    ntp_CNT->SetBranchAddress("dchquality", &dchquality);
    ntp_CNT->SetBranchAddress("alpha",      &alpha);
    ntp_CNT->SetBranchAddress("ppc1",       &ppc1);
    ntp_CNT->SetBranchAddress("ppc2",       &ppc2);
    ntp_CNT->SetBranchAddress("ppc3",       &ppc3);
    ntp_CNT->SetBranchAddress("pc2dphi",    &pc2dphi);
    ntp_CNT->SetBranchAddress("pc2dz",      &pc2dz);
    ntp_CNT->SetBranchAddress("pc3dphi",    &pc3dphi);
    ntp_CNT->SetBranchAddress("pc3dz",      &pc3dz);
    ntp_CNT->SetBranchAddress("dcarm",      &dcarm);
    ntp_CNT->SetBranchAddress("n0",         &n0);



    //==================================================//
    // LOOP OVER TRACKS
    //==================================================//
    cout << endl;
    cout << "--> Looping over CNT Tracks" << endl;

    unsigned int NTRKS = ntp_CNT->GetEntries();
    cout << " Will loop over " << NTRKS << endl;

    //histogram to hold z_bc from tracks in each event
    TH1F *hdczvtx = new TH1F("hdczvtx", ";DC z_{BC} [cm]", 150, -15, 15);

    TH1F *hpczres[3];
    for (int i = 0; i < 3; i++)
    {
        hpczres[i] = new TH1F(
            Form("hpczres_%i", i),
            Form(";PC%i z residual [cm]", i + 1),
            200, -0.2, 0.2);
    }



    TCanvas *cdczvtx = new TCanvas("cdczvtx", "dc zvtx", 1200, 700);
    cdczvtx->Divide(3, 2);
    TLatex ldczvtx;
    ldczvtx.SetNDC();
    ldczvtx.SetTextAlign(22);

    int prevevent = -1;
    int nevents = 0;
    for (unsigned int itrk = 0; itrk < NTRKS; itrk++)
    {
        ntp_CNT->GetEntry(itrk);
        if (itrk % 1000 == 0)
            printf("----> Analyzing entry %i \r", (int)itrk);

        if (prevevent != event && itrk > 0)
        {
            // count the event
            nevents++;
            //if (nevents > 10) break;

            //-- Analyze the results for the previous event
            if (hdczvtx->GetEntries() > 40)
            {
                float zvtx = Median((TH1D *)hdczvtx);
                //float zvtx = hdczvtx->GetMean();


                cdczvtx->cd(2);
                hdczvtx->Draw();
                ldczvtx.DrawLatex(0.5, 0.95, Form("DC z_{vtx} = %.3f cm", zvtx));

                for (int i = 0; i < 3; i++)
                {
                    cdczvtx->cd(i + 4);
                    hpczres[i]->Draw();
                }

                if (nevents % 100 == 0)
                    cdczvtx->Update();
                //gPad->WaitPrimitive();

                //fill event-wise histograms
                hdcz->Fill(zvtx);
                hdcbbcz->Fill(bbcz - zvtx);
                hdcvtxz->Fill(vtx[2] - zvtx);
            }

            //-- Reset persistent information
            hdczvtx->Reset();
            for (int i = 0; i < 3; i++)
                hpczres[i]->Reset();
        }

        //-- Accumulate information for all tracks in the event
        // Make track cuts
        if ( (dchquality == 31 || dchquality == 63)
                && dcarm == 1 //West
                && n0 < 2
                
                && TMath::Abs(pc2dphi) < 0.015
                && TMath::Abs(pc2dz) < 6
                && TMath::Abs(pc3dphi) < 0.02
                && TMath::Abs(pc3dz) < 12
                
                /*
                && TMath::Abs(pc2dphi) < 990
                && TMath::Abs(pc2dz) < 990
                && TMath::Abs(pc3dphi) < 990
                && TMath::Abs(pc3dz) < 990
				*/
           )
        {



            //Fill the matrices and perform the fit
            int m = 3; //PC1, PC2, PC3
            int n = 2;
            TMatrixD X(m, n);
            TMatrixD Cinv(m, m);
            TVectorD y(m);

            for (int i = 0; i < 3; i++)
            {
                X(i, 0) = 1;
                X(i, 1) = pcr[i];
                Cinv(i, i) = 1. / TMath::Power(pcres[i], 2);
                //y(i) = pcr[i];
            }
            //X(0, 1) = ppc1[2];
            //X(1, 1) = ppc2[2];
            //X(2, 1) = ppc3[2];
            y(0) = ppc1[2];
            y(1) = ppc2[2];
            y(2) = ppc3[2];

            TVectorD beta = SolveGLS(X, y, Cinv);

            //calculate the z residuals
            float resid[3];

            resid[0] = beta(1) *
                       //TMath::Sqrt(ppc1[0] * ppc1[0] + ppc1[1] * ppc1[1])
                       pcr[0]
                       + beta(0);
            resid[0] -= ppc1[2];
            hpczres[0]->Fill(resid[0]);


            resid[1] = beta(1) *
                       //TMath::Sqrt(ppc2[0] * ppc2[0] + ppc2[1] * ppc2[1])
                       pcr[1]
                       + beta(0);
            resid[1] -= ppc2[2];
            hpczres[1]->Fill(resid[1]);


            resid[2] = beta(1) *
                       //TMath::Sqrt(ppc3[0] * ppc3[0] + ppc3[1] * ppc3[1])
                       pcr[2]
                       + beta(0);
            resid[2] -= ppc3[2];
            hpczres[2]->Fill(resid[2]);



            //find the z value at the beamcenter
            float r_bc = TMath::Sqrt(
                             TMath::Power(dcbeamcent[0], 2) +
                             TMath::Power(dcbeamcent[1], 2) );
            float z_bc = beta(1) * r_bc + beta(0);

            //fill the histogram with the z_bc
            hdczvtx->Fill(z_bc);


            /*
            cout << "  itrk=" << itrk
                 << " ppc1z=" << ppc1[2]
                 << " ppc2z=" << ppc2[2]
                 << " ppc3z=" << ppc3[2]
                 << " beta[0]=" << beta(0)
                 << " beta[1]=" << beta(1)
                 << " r_bc=" << r_bc
                 << " z_bc=" << z_bc
                 << " rpc1=" << resid[0]
                 << " rpc2=" << resid[1]
                 << " rpc3=" << resid[2]
                 << endl;
            */
            prevevent = event;


        }

    }


    //==================================================//
    // PLOT OBJECTS
    //==================================================//
    cout << endl;
    cout << "--> Plotting" << endl;


    //==================================================//
    // PLOT
    //==================================================//


    TCanvas *cz = new TCanvas("cz", "z vtx", 1200, 600);
    cz->Divide(3, 1);

    cz->cd(1);
    hdcz->Draw();

    cz->cd(2);
    hdcbbcz->Draw();

    cz->cd(3);
    hdcvtxz->Draw();





}