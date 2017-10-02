/*****************  project2 ************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "jacobi_method.h"

using namespace std;
using namespace arma;






int main(int argc, char *argv[])
{
   int n = atoi(argv[argc-1]);
   mat matA=zeros<mat>(n, n);
   mat jacobi_diag_matA=zeros<mat>(n, n);

   double** A;
   A = new double*[n];
   for(int i = 0; i < n;i++){
       A[i]=new double [n];
   }
   //specify matrix A
   for(int i=0;i<n;i++){
       for(int j=0;j<n;j++)
       {
          if(i==j)
          A[i][j]=1;
          else if(fabs(i-j)==1)
            A[i][j]=-2;
          else
              A[i][j]=0;
       }
   }

   for(int i=0;i<n;i++){
       for(int j=0;j<n;j++)
       {
           matA(i,j)=A[i][j];
       }
   }
   cout <<"matA:  \n"<<matA<< endl;
    jacobi_method(A,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
        {
            jacobi_diag_matA(i,j)=A[i][j];
        }
    }
    cout <<"jacobi_diag_matA:  \n"<<jacobi_diag_matA<< endl;
    vec arma_eigval;
    mat arma_eigvec;
    eig_sym(arma_eigval, arma_eigvec, matA);

    arma_eigval.print();
    //arma_eigvec.print();

    delete []A;
    return 0;
}
