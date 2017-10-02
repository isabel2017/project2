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
   int l_obital;
   double ro_m, ro_M ,h;
   mat matA=zeros<mat>(n, n);
   mat jacobi_diag_matA=zeros<mat>(n, n);

   //define variables
   ro_m = 0;
   ro_M = 8.0;
   l_obital = 0;

   //step
   h = ro_M / (double)n;

   double* ro_v = new double[n];
   double* ocilator_potential = new double[n];
   for (int i = 0; i < n; i++) {
	   ro_v[i] = ro_m + i*h;
	   ocilator_potential[i] = ro_v[i] * ro_v[i];
   }

   double* diag = new double[n];
   double* non_diag = new double[n];
   for (int i = 0; i < n; i++) {
	   diag[i] = 2.0 / (h*h) + ocilator_potential[i];
	   non_diag[i] = -1.0 / (h*h);
   }

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
          A[i][j]=diag[i];
          else if(fabs(i-j)==1)
            A[i][j]=non_diag[i];
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
	/*
    vec arma_eigval;
    mat arma_eigvec;
    eig_sym(arma_eigval, arma_eigvec, matA);

    arma_eigval.print();
    //arma_eigvec.print();
	*/
    delete []A;
    return 0;
}
