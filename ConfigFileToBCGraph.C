#include "VtxIO.h"

void ConfigFileToBCGraph()
{
  int run = 421822;
  const char *configfilename =
    Form("production/config/run15pp200/config-zf-%d-0-1.txt", run);

  float bc[2] = {0};
  float e2w[3] = {0};
  float v2c[3] = {0};
  string parFileName;
  GetParamFromConfig(configfilename, bc, e2w, v2c, parFileName);

  TGraphErrors *gbc = new TGraphErrors();
  gbc->SetTitle("Beamcenter from east arm (point 0) and west arm (point 1)");
  gbc->SetPoint(0, bc[0], bc[1]);
  gbc->SetPoint(1, bc[0], bc[1]);
  gbc->SetPointError(0, 0.015, 0.010);
  gbc->SetPointError(1, 0.015, 0.010);

  gbc->Draw("aep");

  TString bcFileOut  = Form("rootfiles/bc-%d.root", run);
  TFile *bcf = new TFile(bcFileOut.Data(), "recreate");
  gbc->Write("gbc");
  bcf->Close();

  return;
}