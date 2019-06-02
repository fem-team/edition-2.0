/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Jinliang Kang promoted                                                */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Element.h"
#include "Q5.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//	Constructor
CQ5::CQ5()
{
	NEN_ = 5; // Each element has 5 nodes
	nodes_ = new CNode * [NEN_];

	ND_ = 10;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CQ5::~CQ5() {}

//	Read element data from stream Input
bool CQ5::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N; // element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "   Expected element : " << Ele + 1 << endl
			<< "   Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;               // Material property set number
	unsigned int N1, N2, N3, N4, N5; // numbers of 4 Nodes

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> MSet;
	ElementMaterial_ = dynamic_cast<CQ5Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	return true;
}

//	Write element data to stream OutputFile
void CQ5::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber << setw(9)
		<< nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber << setw(9)
		<< nodes_[3]->NodeNumber << setw(9) << nodes_[4]->NodeNumber << setw(12)
		<< ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CQ5::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 2; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 5Q element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int CQ5::SizeOfStiffnessMatrix() { return 55; }

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQ5::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//	Calculate element stiffness matrix

	CQ5Material* material_ = dynamic_cast<CQ5Material*>(ElementMaterial_); // Pointer to material of the element

	double E = material_->E;
	double nu = material_->nu;
	double eta, psi, Je; // ×ø±ê×ª»»£¬eta-x£¬psi-y£¬Jaccobi ¾ØÕó
	double J[2][2];
	double JT[2][2];
	double B[2][5];
	double D[4];
	double x[5], y[5];
	double weight;     //w

	x[0] = nodes_[0]->XYZ[0];
	y[0] = nodes_[0]->XYZ[1];
	x[1] = nodes_[1]->XYZ[0];
	y[1] = nodes_[1]->XYZ[1];
	x[2] = nodes_[2]->XYZ[0];
	y[2] = nodes_[2]->XYZ[1];
	x[3] = nodes_[3]->XYZ[0];
	y[3] = nodes_[3]->XYZ[1];
	x[4] = nodes_[4]->XYZ[0];
	y[4] = nodes_[4]->XYZ[1];

	for (unsigned i = 0; i < 3; i++)
	{
		for (unsigned j = 0; j < 3; j++)
		{
			if ((i == 0) && (j == 0)) {
				psi = -sqrt(3.0 / 5.0);
				eta = -sqrt(3.0 / 5.0);
				weight = 25.0 / 81.0;
			}

			if ((i == 1) && (j == 0)) {
				psi = 0.0;
				eta = -sqrt(3.0 / 5.0);
				weight = 40.0 / 81.0;
			}

			if ((i == 2) && (j == 0)) {
				psi = sqrt(3.0 / 5.0);
				eta = -sqrt(3.0 / 5.0);
				weight = 25.0 / 81.0;
			}

			if ((i == 0) && (j == 1)) {
				psi = -sqrt(3.0 / 5.0);
				eta = 0.0;
				weight = 40.0 / 81.0;
			}

			if ((i == 1) && (j == 1)) {
				psi = 0.0;
				eta = 0.0;
				weight = 64.0 / 81.0;
			}

			if ((i == 2) && (j == 1)) {
				psi = sqrt(3.0 / 5.0);
				eta = 0.0;
				weight = 40.0 / 81.0;
			}

			if ((i == 0) && (j == 2)) {
				psi = -sqrt(3.0 / 5.0);
				eta = sqrt(3.0 / 5.0);
				weight = 25.0 / 81.0;
			}

			if ((i == 1) && (j == 2)) {
				psi = 0.0;
				eta = sqrt(3.0 / 5.0);
				weight = 40.0 / 81.0;
			}

			if ((i == 2) && (j == 2)) {
				psi = sqrt(3.0 / 5.0);
				eta = sqrt(3.0 / 5.0);
				weight = 25.0 / 81.0;
			}

			J[0][0] = x[0] * (eta - 1) * eta / 4 + x[1] * (1 + eta) * eta / 4 +
				x[2] * (-1 - eta) / 4 + x[3] * (1 + eta) / 4 + x[4] * (1 - eta * eta) / 2;
			J[0][1] = y[0] * (eta - 1) * eta / 4 + y[1] * (1 + eta) * eta / 4 +
				y[2] * (-1 - eta) / 4 + y[3] * (1 + eta) / 4 + y[4] * (1 - eta * eta) / 2;
			J[1][0] = x[0] * (1 + psi) * (2 * eta - 1) / 4 + x[1] * (1 + psi) * (2 * eta + 1) / 4 +
				x[2] * (1 - psi) / 4 + x[3] * (psi - 1) / 4 + x[4] * (1 + psi) * (-eta);
			J[1][1] = y[0] * (1 + psi) * (2 * eta - 1) / 4 + y[1] * (1 + psi) * (2 * eta + 1) / 4 +
				y[2] * (1 - psi) / 4 + y[3] * (psi - 1) / 4 + y[4] * (1 + psi) * (-eta);
			Je = fabs(J[0][0] * J[1][1] - J[0][1] * J[1][0]);
			JT[0][0] = J[1][1] / Je;
			JT[0][1] = -J[0][1] / Je;
			JT[1][0] = -J[1][0] / Je;
			JT[1][1] = J[0][0] / Je;

			B[0][0] = JT[0][0] * eta * (eta - 1) / 4 + JT[0][1] * (1 + psi) * (2 * eta - 1) / 4;
			B[0][1] = JT[0][0] * eta * (eta + 1) / 4 + JT[0][1] * (1 + psi) * (2 * eta + 1) / 4;
			B[0][2] = JT[0][0] * (-1 - eta) / 4 + JT[0][1] * (1 - psi) / 4;
			B[0][3] = JT[0][0] * (eta - 1) / 4 + JT[0][1] * (psi - 1) / 4;
			B[0][4] = JT[0][0] * (1 - eta * eta) / 2 + JT[0][1] * (1 + psi) * (-eta);
			B[1][0] = JT[1][0] * eta * (eta - 1) / 4 + JT[1][1] * (1 + psi) * (2 * eta - 1) / 4;
			B[1][1] = JT[1][0] * eta * (eta + 1) / 4 + JT[1][1] * (1 + psi) * (2 * eta + 1) / 4;
			B[1][2] = JT[1][0] * (-1 - eta) / 4 + JT[1][1] * (1 - psi) / 4;
			B[1][3] = JT[1][0] * (eta - 1) / 4 + JT[1][1] * (psi - 1) / 4;
			B[1][4] = JT[1][0] * (1 - eta * eta) / 2 + JT[1][1] * (1 + psi) * (-eta);

			D[0] = E / (1.0 - nu * nu);
			D[1] = nu * E / (1.0 - nu * nu);
			D[2] = E / (1.0 - nu * nu);
			D[3] = (1.0 - nu) * E / (2.0 - 2.0 * nu * nu);

			Matrix[0] += weight * Je * (B[0][0] * B[0][0] * D[0] + B[1][0] * B[1][0] * D[3]);
			Matrix[1] += weight * Je * (B[1][0] * B[1][0] * D[2] + B[0][0] * B[0][0] * D[3]);
			Matrix[2] += weight * Je * (B[0][0] * B[1][0] * (D[1] + D[3]));
			Matrix[3] += weight * Je * (B[0][1] * B[0][1] * D[0] + B[1][1] * B[1][1] * D[3]);
			Matrix[4] += weight * Je * (B[1][0] * B[0][1] * D[1] + B[0][0] * B[1][1] * D[3]);
			Matrix[5] += weight * Je * (B[0][0] * B[0][1] * D[0] + B[1][0] * B[1][1] * D[3]);
			Matrix[6] += weight * Je * (B[1][1] * B[1][1] * D[2] + B[0][1] * B[0][1] * D[3]);
			Matrix[7] += weight * Je * (B[0][1] * B[1][1] * D[1] + B[1][1] * B[0][1] * D[3]);
			Matrix[8] += weight * Je * (B[1][0] * B[1][1] * D[2] + B[0][0] * B[0][1] * D[3]);
			Matrix[9] += weight * Je * (B[0][0] * B[1][1] * D[1] + B[1][0] * B[0][1] * D[3]);
			Matrix[10] += weight * Je * (B[0][2] * B[0][2] * D[0] + B[1][2] * B[1][2] * D[3]);
			Matrix[11] += weight * Je * (B[1][1] * B[0][2] * D[1] + B[0][1] * B[1][2] * D[3]);
			Matrix[12] += weight * Je * (B[0][1] * B[0][2] * D[0] + B[1][1] * B[1][2] * D[3]);
			Matrix[13] += weight * Je * (B[1][0] * B[0][2] * D[1] + B[0][0] * B[1][2] * D[3]);
			Matrix[14] += weight * Je * (B[0][0] * B[0][2] * D[0] + B[1][0] * B[1][2] * D[3]);
			Matrix[15] += weight * Je * (B[1][2] * B[1][2] * D[2] + B[0][2] * B[0][2] * D[3]);
			Matrix[16] += weight * Je * (B[0][2] * B[1][2] * D[1] + B[1][2] * B[0][2] * D[3]);
			Matrix[17] += weight * Je * (B[1][1] * B[1][2] * D[2] + B[0][1] * B[0][2] * D[3]);
			Matrix[18] += weight * Je * (B[0][1] * B[1][2] * D[1] + B[1][1] * B[0][2] * D[3]);
			Matrix[19] += weight * Je * (B[1][0] * B[1][2] * D[2] + B[0][0] * B[0][2] * D[3]);
			Matrix[20] += weight * Je * (B[0][0] * B[1][2] * D[1] + B[1][0] * B[0][2] * D[3]);
			Matrix[21] += weight * Je * (B[0][3] * B[0][3] * D[0] + B[1][3] * B[1][3] * D[3]);
			Matrix[22] += weight * Je * (B[1][2] * B[0][3] * D[1] + B[0][2] * B[1][3] * D[3]);
			Matrix[23] += weight * Je * (B[0][2] * B[0][3] * D[0] + B[1][2] * B[1][3] * D[3]);
			Matrix[24] += weight * Je * (B[1][1] * B[0][3] * D[1] + B[0][1] * B[1][3] * D[3]);
			Matrix[25] += weight * Je * (B[0][1] * B[0][3] * D[0] + B[1][1] * B[1][3] * D[3]);
			Matrix[26] += weight * Je * (B[1][0] * B[0][3] * D[1] + B[0][0] * B[1][3] * D[3]);
			Matrix[27] += weight * Je * (B[0][0] * B[0][3] * D[0] + B[1][0] * B[1][3] * D[3]);
			Matrix[28] += weight * Je * (B[1][3] * B[1][3] * D[2] + B[0][3] * B[0][3] * D[3]);
			Matrix[29] += weight * Je * (B[0][3] * B[1][3] * D[1] + B[1][3] * B[0][3] * D[3]);
			Matrix[30] += weight * Je * (B[1][2] * B[1][3] * D[2] + B[0][2] * B[0][3] * D[3]);
			Matrix[31] += weight * Je * (B[0][2] * B[1][3] * D[1] + B[1][2] * B[0][3] * D[3]);
			Matrix[32] += weight * Je * (B[1][1] * B[1][3] * D[2] + B[0][1] * B[0][3] * D[3]);
			Matrix[33] += weight * Je * (B[0][1] * B[1][3] * D[1] + B[1][1] * B[0][3] * D[3]);
			Matrix[34] += weight * Je * (B[1][0] * B[1][3] * D[2] + B[0][0] * B[0][3] * D[3]);
			Matrix[35] += weight * Je * (B[0][0] * B[1][3] * D[1] + B[1][0] * B[0][3] * D[3]);
			Matrix[36] += weight * Je * (B[0][4] * B[0][4] * D[0] + B[1][4] * B[1][4] * D[3]);
			Matrix[37] += weight * Je * (B[0][4] * B[1][3] * D[1] + B[0][3] * B[1][4] * D[3]);
			Matrix[38] += weight * Je * (B[0][3] * B[0][4] * D[0] + B[1][3] * B[1][4] * D[3]);
			Matrix[39] += weight * Je * (B[0][4] * B[1][2] * D[1] + B[0][2] * B[1][4] * D[3]);
			Matrix[40] += weight * Je * (B[0][2] * B[0][4] * D[0] + B[1][2] * B[1][4] * D[3]);
			Matrix[41] += weight * Je * (B[0][4] * B[1][1] * D[1] + B[0][1] * B[1][4] * D[3]);
			Matrix[42] += weight * Je * (B[0][1] * B[0][4] * D[0] + B[1][1] * B[1][4] * D[3]);
			Matrix[43] += weight * Je * (B[0][4] * B[1][0] * D[1] + B[0][0] * B[1][4] * D[3]);
			Matrix[44] += weight * Je * (B[0][0] * B[0][4] * D[0] + B[1][0] * B[1][4] * D[3]);
			Matrix[45] += weight * Je * (B[0][4] * B[0][4] * D[3] + B[1][4] * B[1][4] * D[2]);
			Matrix[46] += weight * Je * (B[0][4] * B[1][4] * D[1] + B[0][4] * B[1][4] * D[3]);
			Matrix[47] += weight * Je * (B[0][3] * B[0][4] * D[3] + B[1][3] * B[1][4] * D[2]);
			Matrix[48] += weight * Je * (B[0][3] * B[1][4] * D[1] + B[0][4] * B[1][3] * D[3]);
			Matrix[49] += weight * Je * (B[0][2] * B[0][4] * D[3] + B[1][2] * B[1][4] * D[2]);
			Matrix[50] += weight * Je * (B[0][2] * B[1][4] * D[1] + B[0][4] * B[1][2] * D[3]);
			Matrix[51] += weight * Je * (B[0][1] * B[0][4] * D[3] + B[1][1] * B[1][4] * D[2]);
			Matrix[52] += weight * Je * (B[0][1] * B[1][4] * D[1] + B[0][4] * B[1][1] * D[3]);
			Matrix[53] += weight * Je * (B[0][0] * B[0][4] * D[3] + B[1][0] * B[1][4] * D[2]);
			Matrix[54] += weight * Je * (B[0][0] * B[1][4] * D[1] + B[0][4] * B[1][0] * D[3]);
		}
	}
}

//	Calculate element stress
void CQ5::ElementStress(double* Q4stress, double* Displacement) {}
void CQ5::ElementGauss(double* Coordinate) {}
// Gravity
double CQ5::Gravity()
{
	return 0;
}