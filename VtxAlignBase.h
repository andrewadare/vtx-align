#ifndef __VTXALIGNBASE_H__
#define __VTXALIGNBASE_H__

#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "SvxProj.h"
#include "Mille.h"
#include "UtilFns.h"

#include <TFile.h>
#include <TSystem.h>
#include <TGeoManager.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

typedef vector<double> vecd;
typedef vector<int>    veci;
typedef vector<string> vecs;
typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

SvxTGeo *VTXModel(const char *pisafile);

SvxTGeo *
VTXModel(const char *pisafile)
{
  SvxTGeo *tgeo = new SvxTGeo;
  tgeo->ReadParFile(pisafile);
  tgeo->MakeTopVolume(100, 100, 100);
  tgeo->AddSensors();

  TGeoManager *mgr = tgeo->GeoManager();
  mgr->CloseGeometry();

  return tgeo;
}

#endif
