#include "VtxAlignBase.h"
#include "VtxIO.h"
#include <utility>
#include <map>

using namespace std;

void test_eventsync()
{

  int ntracks = 1000;

  TFile *inFile = new TFile("rootfiles/411768-0-0_3.root", "read");

  SvxTGeo *tgeo = VTXModel("geom/411768-0-0.par");

  TNtuple *svxseg = (TNtuple *)inFile->Get("vtxhits");
  assert(svxseg);

  TNtuple *svxcnt = (TNtuple *)inFile->Get("cnthits");
  assert(svxcnt);


  multimap<int, int> seg_eventsync;
  geoTracks seg_tracks;

  GetTracksFromTree(svxseg, tgeo, seg_tracks, seg_eventsync , ntracks, "");

  geoTracks cnt_tracks;
  multimap<int, int> cnt_eventsync;

  GetTracksFromTree(svxcnt, tgeo, cnt_tracks, cnt_eventsync , ntracks, "");


  for(unsigned int i = 0; i < seg_tracks.size();i++)//loop thru standalone tracks
    {
      multimap<int,int>::iterator it = seg_eventsync.find(i);
      int event = (*it).second;

      for(unsigned int j = 0; j < cnt_tracks.size();j++)//loop thru cnt tracks
	{
	  multimap<int,int>::iterator it2 = cnt_eventsync.find(j);
	  if(event!=(*it2).second)//check to see if they have the same event number
	    {
	      continue;
	    }
	  else
	    {
	      //do stuff here
	    }
	}
    }


}
