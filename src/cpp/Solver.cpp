/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"
#include "Outputter.h"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#ifdef MKL
#include "mkl.h"
#endif

using namespace std;

CSolver::CSolver(CSkylineMatrix<double>* K) : K(K) {};

// LDLT facterization
void CLDLTSolver::LDLT()
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights
//	ofstream out("E:/GitHub/edition-1.2/data/file.txt");

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
        // Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += (*K)(r,i) * (*K)(r,j);		// C += L_ri * U_rj

			(*K)(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = (*K)(r,j) / (*K)(r,r);	// L_rj = U_rj / D_rr
			(*K)(j,j) -= Lrj * (*K)(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			(*K)(r,j) = Lrj;
		}

        if (fabs((*K)(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << (*K)(j,j) << endl;
            
            exit(4);
        }
		
//		out << "j=" << j << setw(14) << (*K)(j, j) << endl;
    }
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights
	/*
	ofstream out("E:/GitHub/edition-1.2/data/file.txt");
	out << "Force before calculation" << endl;
	for( int i=1; i<N;i++)
		out << "Force" << i - 1 << setw(14) << Force[i - 1] << endl;
	out << "Force after calculation" << endl;
	*/

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= (*K)(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
//		out << "Force" << i - 1 << setw(14) << Force[i-1] << endl;
	}
	

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= (*K)(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= (*K)(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};

#ifdef MKL
void CSRSolver::solve(double* Force, unsigned NLCase)
{
	void* pt[64];
	for (unsigned _ = 0; _ < 64; _++) pt[_] = 0;

	const int mtype = 2;
	int iparm[64] = { 0 };

	pardisoinit(pt, &mtype, iparm);
	iparm[1] = 2; // The parallel (OpenMP) version of the nested dissection algorithm.
	iparm[5] = 1; // write back to Force
	iparm[59] = 1; // use OOC if needed

	const int one = 1;
	const int size = K.size;
	double* values = K.values;
	int* columns = K.columns;
	int* rowIndexs = K.rowIndexs;

	const int rhsCount = NLCase;
	double* rhs = Force;

	int phase = 13;
	double* res = new double[rhsCount*size];
	for (std::size_t _ = 0; _ < rhsCount*size; _++) res[_] = 0;

	int msglvl = 0; // print info
#if defined(_DEBUG_) || defined(_RUN_)
	msglvl = 1;
#endif // _DEBUG_
	int* perm = new int[size];
	int error;

	pardiso(
		pt, // handle to some shit
		&one, // maxfct
		&one, // mnum
		&mtype, // sym pos matrix
		&phase, // go through all
		&size,  // size of matrix
		values,
		rowIndexs,
		columns,
		perm, // idk wtf this is
		&rhsCount,
		iparm, // sort like settings
		&msglvl, // print info or not
		rhs,
		res,
		&error // see if any error
	);
	if (error)
	{
		std::cerr << "ERROR IN PARDISO SOLVER: " << error << std::endl;
		exit(8);
	}

#ifdef _DEBUG_
	for (int _ = 0; _ < size; _++)
		std::cout << "res[" << _ << "] = " << res[_] << std::endl;
	for (int _ = 0; _ < size; _++)
		std::cout << "rhs[" << _ << "] = " << rhs[_] << std::endl;
#endif // _DEBUG_
	delete[] perm;
	delete[] res;
}
#endif

#ifdef _VIB_
//! only to be used before LDLT() factorization
void CLDLTSolver::Multiple(double* acc, double* force, unsigned int numeq, unsigned int vib_m) {
	unsigned int* diag = K.GetDiagonalAddress();
	unsigned int* colh = K.GetColumnHeights();
	unsigned int diag_sc = 0;
	unsigned int row_ck = 0;

	for (unsigned int i = 0; i < numeq*vib_m; ++i) {
		force[i] = 0.0;
	}

	for (unsigned int i = 1; i < diag[numeq]; ++i) {

		if (i == diag[diag_sc]) {
			for (unsigned int j = 0; j < vib_m; ++j) {
				force[j*numeq + diag_sc] += K(diag_sc + 1, diag_sc + 1)*acc[j*numeq + diag_sc];
			}


		}
		else {
			for (unsigned int j = 0; j < vib_m; ++j) {
				force[j*numeq + row_ck] += K(diag_sc + 1, row_ck + 1)*acc[j*numeq + diag_sc];
				force[j*numeq + diag_sc] += K(diag_sc + 1, row_ck + 1)*acc[j*numeq + row_ck];
			}

		}
		if (diag_sc - row_ck < colh[diag_sc]) row_ck--;
		else {
			diag_sc++;
			row_ck = diag_sc;
		}

	}

}

#endif
