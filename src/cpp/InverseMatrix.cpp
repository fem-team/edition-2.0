#include<cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include "InverseMatrix.h"

//#include <wincrypt.h>
using namespace std;

void Inverse(double A[12][12], double B[12][12], int M)
{
	static const int N = 12;
	double s[N][N];
	double m[N][N];
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
			m[i][j] = A[i][j];
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 12; j++)
			s[i][j] = 0.0;
		s[i][i] = 1.0;
	}


	for (unsigned column = 0; column < N; ++column) {
		// Swap row in case our pivot point is not working
		if (m[column][column] == 0) {
			unsigned big = column;
			for (unsigned row = 0; row < N; ++row)
				if (fabs(m[row][column]) > fabs(m[big][column])) big = row;
			// Print this is a singular matrix, return identity ?
			if (big == column) fprintf(stderr, "Singular matrix\n");
			// Swap rows                               
			else for (unsigned j = 0; j < N; ++j) {
				std::swap(m[column][j], m[big][j]);
				std::swap(s[column][j], s[big][j]);
			}
		}
		// Set each row in the column to 0  
		for (unsigned row = 0; row < N; ++row) {
			if (row != column) {
				double coeff = m[row][column] / m[column][column];
				if (coeff != 0) {
					for (unsigned j = 0; j < N; ++j) {
						m[row][j] -= coeff * m[column][j];
						s[row][j] -= coeff * s[column][j];
					}
					// Set the element to 0 for safety
					m[row][column] = 0.0;
				}
			}
		}
	}
	// Set each element of the diagonal to 1
	for (unsigned row = 0; row < N; ++row) {
		for (unsigned column = 0; column < N; ++column) {
			s[row][column] /= m[row][row];
		}
	}
	for (int i = 0; i < 12; i++)
	{
		for (int j = 0; j < 12; j++)
			B[i][j] = s[i][j];
	}
}
