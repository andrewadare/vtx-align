#include "Mille.h"
#include "TRandom3.h"

void WriteSteerFile(const char *filename, const char *binfile);

void SimpleExample()
{
  // Align a single detector plane in one direction (y) using residuals dy.
  WriteBinFile("example.bin");
  WriteSteerFile("steer.txt", "example.bin");
  gSystem->Exec("pede steer.txt");
  return;
}

void
WriteBinFile(const char *filename)
{
  // Put one global parameter in the list.
  // The label is an arbitrary positive integer.
  vector<int> labels = {333};

  // One local derivative and one global derivative: d(Delta_y)/dy = 1
  vector<float> ld = {1.0};
  vector<float> gd = {1.0};

  // Width and mean of residual distribution.
  // The mean is what Millepede is supposed to recover.
  float sigma = 0.001;
  float mean = 0.01;

  Mille m(filename);
  TRandom3 ran;

  int nres = 10000;
  while (nres > 0)
  {
    // Residual
    float dy = mean + sigma*ran.Gaus();

    m.mille(0, &ld[0], 1, &gd[0], &labels[0], dy, sigma);
    m.mille(1, &ld[0], 0, &gd[0], &labels[0], 0,  sigma);
    m.end();
    --nres;
  }

  return;
}

void
WriteSteerFile(const char *filename, const char *binfile)
{
  cout << "Writing " << filename << "..." << flush;

  ofstream fs(filename);

  fs << Form("* This is %s", filename) << endl;
  fs << Form("* Pass this file to pede: pede %s", filename) << endl;
  fs << endl;

  fs << "Fortranfiles  ! Fortran/text inputs listed here:" << endl;
  // <None>
  fs << endl;

  fs << "Cfiles  ! c/c++ binary input files listed here:" << endl;
  fs << binfile << " ! binary data file" << endl;
  fs << endl;

  // Uncomment to change from default of 10
  // fs << "entries  1 ! lower limit on number of entries/parameter" << endl;

  fs << "method inversion 5 0.0001  ! Gauss. elim., #iterations, tol." << endl;

  fs << "end" << endl;

  fs.close();

  Printf("done.");

  return;
}
