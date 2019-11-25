#include <iostream>
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"

{

int N_bin = 50;

ifstream in("corr_matrix.txt");

TArrayD data(N_bin*N_bin);
TMatrixD V(N_bin,N_bin);

for (int i=0; i<N_bin*N_bin; i++) {
    in >> data[i];
    }

V.SetMatrixArray(data.GetArray());
//V_inv = V.Invert();

TDecompLU lu(V);
cout << "cond: " << lu.Condition() << endl;

}
