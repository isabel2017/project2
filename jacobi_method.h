#ifndef JACOBI_METHOD_H
#define JACOBI_METHOD_H

#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
using namespace std;

double max_nondiag(double** A,int* k, int* l,int n);
void jacobi_rotation(double** A,int k, int l,int n);
void jacobi_method(double** A, int n);







#endif // JACOBI_METHOD_H
