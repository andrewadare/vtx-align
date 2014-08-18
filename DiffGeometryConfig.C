#include "VtxAlignBase.h"

#include "ParameterDefs.h"
#include "BadLadders.h"
#include "VtxVis.h"

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z);
SvxTGeo *GetXYZOffset(const char *parfile, vecd &x, vecd &y, vecd &z, float eastToWest[], float vtxToCNT[]);

/**
 *
 * Compare two alignments based on config files
 *
 */
void DiffGeometryConfig(const char *confa = "production/config/config-fieldon-407951-0-4.txt",
                        const char *confb = "production/config/config-taebong-p2-v8.txt"
                                //const char *confb = "production/config/config-offset-test.txt"
                       )
{
    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.03);

    ifstream fin;
    string line;
    char pfa[80];
    char pfb[80];

    //-- Read config file a
    float beamcentera[2];
    float eastToWesta[3];
    float vtxToCNTa[3];
    char geoFilea[80];

    cout << endl;
    cout << "--> Reading config file A from " << confa << endl;
    fin.open(confa);
    assert(fin);

    while (getline(fin, line))
    {
        if (line.find("beamcenter:") != string::npos)
        {
            sscanf(line.c_str(), "beamcenter: %f %f", &beamcentera[0], &beamcentera[1]);
            cout << "beamcenter=("
                 << beamcentera[0] << ", "
                 << beamcentera[1] << ")"
                 << endl;

        }
        if (line.find("east-to-west:") != string::npos)
        {
            sscanf(line.c_str(), "east-to-west: %f %f %f", &eastToWesta[0], &eastToWesta[1], &eastToWesta[2]);
            cout << "east-to-west=("
                 << eastToWesta[0] << ", "
                 << eastToWesta[1] << ", "
                 << eastToWesta[2] << ")"
                 << endl;
        }
        if (line.find("vtx-to-cnt:") != string::npos)
        {
            sscanf(line.c_str(), "vtx-to-cnt: %f %f %f", &vtxToCNTa[0], &vtxToCNTa[1], &vtxToCNTa[2]);
            cout << "vtx-to-cnt=("
                 << vtxToCNTa[0] << ", "
                 << vtxToCNTa[1] << ", "
                 << vtxToCNTa[2] << ")"
                 << endl;
        }
        if (line.find("geomfile:") != string::npos)
        {
            sscanf(line.c_str(), "geomfile: %s", geoFilea);
            sprintf(pfa, "geom/%s", geoFilea);
            cout << "geomFile= " << pfa << endl;
        }
    }
    fin.close();


    //-- Read config file b
    float beamcenterb[2];
    float eastToWestb[3];
    float vtxToCNTb[3];
    char geoFileb[80];

    cout << endl;
    cout << "--> Reading config file B from " << confb << endl;
    fin.open(confb);
    assert(fin);

    while (getline(fin, line))
    {
        if (line.find("beamcenter:") != string::npos)
        {
            sscanf(line.c_str(), "beamcenter: %f %f", &beamcenterb[0], &beamcenterb[1]);
            cout << "beamcenter=("
                 << beamcenterb[0] << ", "
                 << beamcenterb[1] << ")"
                 << endl;

        }
        if (line.find("east-to-west:") != string::npos)
        {
            sscanf(line.c_str(), "east-to-west: %f %f %f", &eastToWestb[0], &eastToWestb[1], &eastToWestb[2]);
            cout << "east-to-west=("
                 << eastToWestb[0] << ", "
                 << eastToWestb[1] << ", "
                 << eastToWestb[2] << ")"
                 << endl;
        }
        if (line.find("vtx-to-cnt:") != string::npos)
        {
            sscanf(line.c_str(), "vtx-to-cnt: %f %f %f", &vtxToCNTb[0], &vtxToCNTb[1], &vtxToCNTb[2]);
            cout << "vtx-to-cnt=("
                 << vtxToCNTb[0] << ", "
                 << vtxToCNTb[1] << ", "
                 << vtxToCNTb[2] << ")"
                 << endl;
        }
        if (line.find("geomfile:") != string::npos)
        {
            sscanf(line.c_str(), "geomfile: %s", geoFileb);
            cout << "   " << geoFileb << endl;
            sprintf(pfb, "geom/%s", geoFileb);
            cout << "geomFile= " << pfb << endl;
        }
    }
    fin.close();




    // Get ladder positions in "a" and "b" geometry
    vecd xa; vecd ya; vecd za;
    vecd xb; vecd yb; vecd zb;
    //GetXYZ(pfa, xa, ya, za);
    //GetXYZ(pfb, xb, yb, zb);
    SvxTGeo *geoa = GetXYZOffset(pfa, xa, ya, za, eastToWesta, vtxToCNTa);
    SvxTGeo *geob = GetXYZOffset(pfb, xb, yb, zb, eastToWestb, vtxToCNTb);

    //SvxTGeo *geo = VTXModel(pfa);

    TString a = TString(TString(gSystem->BaseName(confa)).ReplaceAll(".txt", "")).ReplaceAll("production/config/", "");
    TString b = TString(TString(gSystem->BaseName(confb)).ReplaceAll(".txt", "")).ReplaceAll("production/config/", "");
    TString ab = Form("%s-vs-%s", a.Data(), b.Data());

    DrawXY(geoa, "s_diff", "#Deltas_{ab}", "L, dead, faint");
    DrawDiffs(xa, ya, za, xb, yb, zb, "s");
    ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/dsxy%s.pdf", ab.Data()));

    DrawXY(geoa, "z_diff", "#Deltaz_{ab}", "L, dead, faint");
    DrawDiffs(xa, ya, za, xb, yb, zb, "z");
    ltx.DrawLatex(0.6, 0.95, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/dzxy%s.pdf", ab.Data()));

    ltx.SetTextSize(0.05);

    DrawDiffsLinear(xa, ya, za, xb, yb, zb, "xdiff", "#Delta x", "x"/*, 0.028*/);
    ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/dx%s.pdf", ab.Data()));

    DrawDiffsLinear(xa, ya, za, xb, yb, zb, "ydiff", "#Delta y", "y"/*, 0.028*/);
    ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/dy%s.pdf", ab.Data()));

    DrawDiffsLinear(xa, ya, za, xb, yb, zb, "zdiff", "#Delta z", "z"/*, 0.028*/);
    ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/dz%s.pdf", ab.Data()));

    DrawDiffsLinear(xa, ya, za, xb, yb, zb, "sdiff", "#Delta s", "s"/*, 0.083*/);
    ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/ds%s.pdf", ab.Data()));

    DrawDiffsLinear(xa, ya, za, xb, yb, zb, "rdiff", "#Delta r", "r"/*, 0.083*/);
    ltx.DrawLatex(0.655, 0.92, Form("#splitline{a: %s}{b: %s}", a.Data(), b.Data()));
    gPad->Print(Form("pdfs/dr%s.pdf", ab.Data()));

    return;
}

void GetXYZ(const char *parfile, vecd &x, vecd &y, vecd &z)
{
    SvxTGeo *geo = VTXModel(parfile);
    GetLadderXYZ(geo, x, y, z);
    delete geo;
    return;
}


SvxTGeo *GetXYZOffset(const char *parfile, vecd &x, vecd &y, vecd &z, float eastToWest[], float vtxToCNT[])
{
    SvxTGeo *geo = VTXModel(parfile);
    //-- shift the ladders by the offsets
    for (int i = 0; i < geo->GetNLayers(); i++)
        for (int j = 0; j < geo->GetNLadders(i); j++)
        {
            //West arm - only gets shifted by vtxToCNT
            if ( (i == 0 && j < 5) ||
                    (i == 1 && j < 10) ||
                    (i == 2 && j < 8) ||
                    (i == 3 && j < 12)
               )
            {
                geo->TranslateLadder(i, j,
                                     vtxToCNT[0],
                                     vtxToCNT[1],
                                     vtxToCNT[2]);
            }
            //East arm - gets shifted by eastToWest+vtxToCNT
            else
            {
                geo->TranslateLadder(i, j,
                                     vtxToCNT[0] + eastToWest[0],
                                     vtxToCNT[1] + eastToWest[1],
                                     vtxToCNT[2] + eastToWest[2]);
            }
        }
    GetLadderXYZ(geo, x, y, z);
    //delete geo;
    return geo;
}
