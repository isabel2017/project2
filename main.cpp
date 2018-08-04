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
	int n = atoi(argv[argc - 1]);
	int l_orbital;
	double ro_m, ro_M, h;
	mat matA = zeros<mat>(n, n);
	mat jacobi_diag_matA = zeros<mat>(n, n);
	mat matA2 = zeros<mat>(n, n);
	mat jacobi_diag_matA2 = zeros<mat>(n, n);
	double* ro_v = new double[n];

	//define variables
	ro_m = 0;
	ro_M = 8.0;
	l_orbital = 0;

	//step
	h = ro_M / (double)n;

	//flags
	//-o : one particle problem
		//-l : lowest three eigen value
		//-a : use armadillo function
	//-t : two particle problem

	// check flag
	for (int i = 1; i < argc; i++) {
		//one particle problem
		if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("o") != string::npos)) {
			
			double* ocilator_potential = new double[n];
			double* diag = new double[n];
			double* non_diag = new double[n];

			for (int i = 0; i < n; i++) {
				ro_v[i] = ro_m + i*h;
				ocilator_potential[i] = ro_v[i] * ro_v[i];
			}

			
			for (int i = 0; i < n; i++) {
				diag[i] = 2.0 / (h*h) + ocilator_potential[i];
				non_diag[i] = -1.0 / (h*h);
			}
			double** A;
			A = new double*[n];
			for (int i = 0; i < n; i++) {
				A[i] = new double[n];
			}
			//specify matrix A
			for (int i = 0; i<n; i++) {
				for (int j = 0; j<n; j++)
				{
					if (i == j)
						A[i][j] = diag[i];
					else if (fabs(i - j) == 1)
						A[i][j] = non_diag[i];
					else
						A[i][j] = 0;
				}
			}

			for (int i = 0; i<n; i++) {
				for (int j = 0; j<n; j++)
				{
					matA(i, j) = A[i][j];
				}
			}
			cout << "matA:  \n" << matA << endl;
			jacobi_method(A, n);

			//copy matrix to arma matrix
			for (int i = 0; i<n; i++) {
				for (int j = 0; j<n; j++)
				{
					jacobi_diag_matA(i, j) = A[i][j];
					if (i == j) {
						//keep diagonal elements in a seperate matrix
						diag[i] = jacobi_diag_matA(i, j);
					}
				}
			}

			cout << "jacobi_diag_matA:  \n" << jacobi_diag_matA << endl;

			delete[]A;

			//check flag -l and -a
			for (int j = 1; j < argc; j++) {
				if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("l") != string::npos)) {
					//in order to find lowest 3 eigen values
					//repeat finding the min
					int i = 0;
					int check = 0;
					int cnt = n;
					double temp;
					double min;
					double lowest[3];

					while (cnt != n - 3) {
						min = diag[0];
						for (i = 0; i < cnt; i++) {
							if (min > diag[i]) {
								min = diag[i];
								check = i;
							}
						}
						lowest[n - cnt] = min;

						temp = diag[check];
						diag[check] = diag[cnt - 1];
						diag[cnt - 1] = temp;

						cnt--;
					}
					cout << endl << "lowest 3 : ";
					for (i = 0; i < 3; i++) {
						cout << lowest[i] << "     ";
					}
					cout << endl << endl;
				}
				else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("a") != string::npos)) {
					//eigen value by armadillo
					/*
					vec arma_eigval;
					mat arma_eigvec;
					eig_sym(arma_eigval, arma_eigvec, matA);

					arma_eigval.print();
					//arma_eigvec.print();
					*/
				}
			}
			
		}

		//two particle problem
		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("t") != string::npos)) {
			//two particle problem
			double w_r;
			l_orbital = 0;

			cout << "choose 0.01, 0.5, 1.0" << endl;
			cin >> w_r;

			double* ocilator_potential_2 = new double[n];
			for (int i = 0; i < n; i++) {
				ro_v[i] = ro_m + i*h;
				if (i == 0) {//since 1.0/ro_v[0] is inf
					ocilator_potential_2[i] = 0;//?bcus ro_m is 0?
				}
				else {
					ocilator_potential_2[i] = w_r * w_r * ro_v[i] * ro_v[i] + (1.0 / ro_v[i]);
				}
			}

			double* diag_2 = new double[n];
			double* non_diag_2 = new double[n];
			for (int i = 0; i < n; i++) {
				diag_2[i] = (2.0 / (h*h)) + ocilator_potential_2[i];
				non_diag_2[i] = -1.0 / (h*h);

			}

			double** A_2;
			A_2 = new double*[n];
			for (int i = 0; i < n; i++) {
				A_2[i] = new double[n];
			}
			//specify matrix A_2
			for (int i = 0; i<n; i++) {
				for (int j = 0; j<n; j++)
				{
					if (i == j)
						A_2[i][j] = diag_2[i];
					else if (fabs(i - j) == 1)
						A_2[i][j] = non_diag_2[i];
					else
						A_2[i][j] = 0;
				}
			}

			for (int i = 0; i<n; i++) {
				for (int j = 0; j<n; j++)
				{
					matA2(i, j) = A_2[i][j];
				}
			}
			cout << "matA2:  \n" << matA2 << endl;
			jacobi_method(A_2, n);

			//copy matrix to arma matrix
			for (int i = 0; i<n; i++) {
				for (int j = 0; j<n; j++)
				{
					jacobi_diag_matA2(i, j) = A_2[i][j];
					if (i == j) {
						//keep diagonal elements in a seperate matrix
						diag_2[i] = jacobi_diag_matA2(i, j);
					}
				}
			}
			cout << "diag_matA2:  \n" << jacobi_diag_matA2 << endl;
			delete[]A_2;
		}
	}
	
    return 0;
}
