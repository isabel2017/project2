//unit test project

#include "catch.hpp"
#include "../project2/jacobi_method.h"

//unit test
TEST_CASE("max_nondiag") {
	//test max_nondiag
	int size = 3;
	double max;

	double** mat = new double*[size];
	for (int i = 0; i < size; i++) {
		mat[i] = new double[size];
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j <= i; j++) {
			mat[i][j] = 5 * i + j;
			mat[j][i] = 5 * i + j;
		}
	}
	//mat[0] = {0, 5, 10}
	//mat[1] = {5, 6, 11}
	//mat[2] = {10,11,12}

	int p, q;

	max = max_nondiag(mat, &p, &q, size);

	REQUIRE((max == 11));

	for (int i = 0; i < size; i++) {
		delete[] mat[i];
	}
	delete[] mat;

}

TEST_CASE("jacobi_method") {
	//test jacobi_method
	int size = 3;

	double epsilon = 1.0e-8;

	double** mat = new double*[size];
	for (int i = 0; i < size; i++) {
		mat[i] = new double[size];
	}

	mat[0][0] = 6;
	mat[0][1] = -2;
	mat[0][2] = -1;
	mat[1][0] = -2;
	mat[1][1] = 6;
	mat[1][2] = -1;
	mat[2][0] = -1;
	mat[2][1] = -1;
	mat[2][2] = 5;

	jacobi_method(mat, size);

	REQUIRE(mat[0][0] - 3 < epsilon);
	REQUIRE((mat[0][1] < epsilon));
	REQUIRE(mat[0][2] < epsilon);
	REQUIRE(mat[1][0] < epsilon);
	REQUIRE(mat[1][1] - 8 < epsilon);
	REQUIRE(mat[1][2] < epsilon);
	REQUIRE(mat[2][0] < epsilon); 
	REQUIRE(mat[2][1] < epsilon);
	REQUIRE(mat[2][2] - 6 < epsilon);

	for (int i = 0; i < size; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}