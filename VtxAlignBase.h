#ifndef __VTXALIGNBASE_H__
#define __VTXALIGNBASE_H__

#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "SvxProj.h"
#include "Mille.h"
#include "UtilFns.h"

// #include "GLSFitter.h"
// #include "BeamCenterFuncs.h"
// #include "ParameterDefs.h"
// #include "ConstraintBuilder.h"
// #include "BadLadders.h"
// #include "VtxIO.h"
// #include "VtxVis.h"

#include <TFile.h>
#include <TSystem.h>
#include <TGeoManager.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TLatex.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

typedef vector<double> vecd;
typedef vector<int>    veci;
typedef vector<string> vecs;
typedef vector<SvxGeoTrack> geoTracks;
typedef vector<geoTracks> geoEvents;

const char *HITVARS = "layer:ladder:sensor:lx:ly:lz:gx:gy:gz:"
                      "x_size:z_size:res_z:res_s:trkid:event";

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
