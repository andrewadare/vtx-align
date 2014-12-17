#include "VtxAlignBase.h"
#include "ParameterDefs.h"

// void ModifyParFile(const char *pfa = "geom/svxPISA-hubert.par",
//                    const char *pfb = "geom/svxPISA-hubert-mod.par")
void ModifyParFile(const char *pfa = "geom/411768-22-11.par",
                   const char *pfb = "geom/411768-22-12.par")
{

  // For 1-3-minus --> 1-3-mm
  // double dr0[10] = {0};
  // dr0[0]=-26.2435;
  // dr0[1]=-63.9824;
  // dr0[2]=-42.9239;
  // dr0[3]=23.0914;
  // dr0[4]=-67.7311;
  // dr0[5]=25.6575;
  // dr0[6]=49.7413;
  // dr0[7]=89.6801;
  // dr0[8]=30.111;
  // dr0[9]=99.5612;

  // double dr1[20] = {0};
  // dr1[0]=-17.2166;
  // dr1[1]=-9.06245;
  // dr1[2]=-82.1011;
  // dr1[3]=-117.766;
  // dr1[4]=-45.7159;
  // dr1[5]=-255.079;
  // dr1[6]=54.2334;
  // dr1[7]=30.3887;
  // dr1[8]=-87.8737;
  // dr1[9]=-71.7099;
  // dr1[10]=-50.3088;
  // dr1[11]=0;
  // dr1[12]=39.1525;
  // dr1[13]=-0.360942;
  // dr1[14]=155.995;
  // dr1[15]=15.3893;
  // dr1[16]=-32.6419;
  // dr1[17]=-7.11585;
  // dr1[18]=52.0381;
  // dr1[19]=111.116;

  // For 1-3 --> 1-3-minus
  // double dr0[10] = {0};
  // dr0[0]=-57.0658;
  // dr0[1]=-121.838;
  // dr0[2]=-68.8414;
  // dr0[3]=43.2215;
  // dr0[4]=-96.8919;
  // dr0[5]=49.3943;
  // dr0[6]=116.691;
  // dr0[7]=168.643;
  // dr0[8]=75.1137;
  // dr0[9]=163.471;

  // double dr1[20] = {0};
  // dr1[0]=-90.2062;
  // dr1[1]=-53.3554;
  // dr1[2]=-244.663;
  // dr1[3]=-277.171;
  // dr1[4]=-147.133;
  // dr1[5]=-411.86;
  // dr1[6]=150.925;
  // dr1[7]=73.8811;
  // dr1[8]=-222.157;
  // dr1[9]=-230.166;
  // dr1[10]=-131.434;
  // dr1[11]=0;
  // dr1[12]=62.4062;
  // dr1[13]=-6.25391;
  // dr1[14]=215.383;
  // dr1[15]=77.837;
  // dr1[16]=-97.7478;
  // dr1[17]=8.99105;
  // dr1[18]=165.61;
  // dr1[19]=182.292;

  SvxTGeo *geo = VTXModel(pfa);

  // geo->MoveLadderRadially(1, 5, +0.015);
  // geo->MoveLadderRadially(1, 14, -0.01);
  // for (int ldr=0; ldr<10; ldr++)
  //   geo->MoveLadderRadially(0, ldr, -1e-4*dr0[ldr]);
  // for (int ldr=0; ldr<20; ldr++)
  //   geo->MoveLadderRadially(1, ldr, -1e-4*dr1[ldr]);
  // geo->MoveLadderRadially(2, 2, -0.2950);
  // geo->MoveLadderRadially(2, 5, -0.2950);
  // geo->MoveLadderRadially(2, 10, -0.2950);
  // geo->MoveLadderRadially(2, 13, -0.2950);


  // for (int i = 0; i < geo->GetNLayers(); i++)
  //   for (int j = 0; j < geo->GetNLadders(i); j++)
  //   {
  //     int l = -1;
  //     // Ladder Phi correction from ds
  //     l = Label(i, j, "s");
  //     if (mpc.find(l) != mpc.end())
  //       geo->RotateLadderRPhi(i, j, mpc[l]);
  //   }

  geo->RotateLadderRPhi(0, 5, -0.0196884   );
  geo->RotateLadderRPhi(0, 6, -0.00965128  );
  geo->RotateLadderRPhi(0, 7, -0.00153134  );
  geo->RotateLadderRPhi(0, 8, +0.00229421 );  
  geo->RotateLadderRPhi(0, 9, +0.00563761 );  

  geo->RotateLadderRPhi(1, 10, -0.0429664  );  
  geo->RotateLadderRPhi(1, 11, -0          ); 
  geo->RotateLadderRPhi(1, 12, -0.0251557  ); 
  geo->RotateLadderRPhi(1, 13, -0.0142688  ); 
  geo->RotateLadderRPhi(1, 14, -0.00781458 );  
  geo->RotateLadderRPhi(1, 15, +0.00302756);
  geo->RotateLadderRPhi(1, 16, +0.00582653);
  geo->RotateLadderRPhi(1, 17, +0.00722349);
  geo->RotateLadderRPhi(1, 18, +0.00816761);

  geo->RotateLadderRPhi(2, 8,  -0.0780687   );   
  geo->RotateLadderRPhi(2, 9,  -0.0619226   );    
  geo->RotateLadderRPhi(2, 10, -0.0427677  );    
  geo->RotateLadderRPhi(2, 11, -0.012706   );   
  geo->RotateLadderRPhi(2, 12, +0.00821891);
  geo->RotateLadderRPhi(2, 13, +0.0188989 );  
  geo->RotateLadderRPhi(2, 14, +0.0178322 );  

  geo->RotateLadderRPhi(3, 13, -0.120821   );  
  geo->RotateLadderRPhi(3, 14, -0.0911046  );   
  geo->RotateLadderRPhi(3, 15, -0.0759552  );   
  geo->RotateLadderRPhi(3, 16, 0          );   
  geo->RotateLadderRPhi(3, 17, +0.0187691  );   
  geo->RotateLadderRPhi(3, 18, +0.00645408);   
  geo->RotateLadderRPhi(3, 19, +0.0260592 );  
  geo->RotateLadderRPhi(3, 20, +0.0252907 );  
  geo->RotateLadderRPhi(3, 21, +0.0401362 );  








  Printf("Writing %s", pfb);
  geo->WriteParFile(pfb);

  return;
}

