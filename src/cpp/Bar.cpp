/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBar::CBar()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CBar::~CBar()
{
}

//	Read element data from stream Input
bool CBar::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
    ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBar::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CBar::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CBar::SizeOfStiffnessMatrix() { return 21; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

//	Calculate element stiffness matrix

	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double k = material_->E * material_->Area / L / L2;

	Matrix[0] = k*DX2[0];
	Matrix[1] = k*DX2[1];
	Matrix[2] = k*DX2[3];
	Matrix[3] = k*DX2[2];
	Matrix[4] = k*DX2[4];
	Matrix[5] = k*DX2[5];
	Matrix[6] = k*DX2[0];
	Matrix[7] = -k*DX2[5];
	Matrix[8] = -k*DX2[3];
	Matrix[9] = -k*DX2[0];
	Matrix[10] = k*DX2[1];
	Matrix[11] = k*DX2[3];
	Matrix[12] = -k*DX2[4];
	Matrix[13] = -k*DX2[1];
	Matrix[14] = -k*DX2[3];
	Matrix[15] = k*DX2[2];
	Matrix[16] = k*DX2[4];
	Matrix[17] = k*DX2[5];
	Matrix[18] = -k*DX2[2];
	Matrix[19] = -k*DX2[4];
	Matrix[20] = -k*DX2[5];
}

//	Calculate element stress 
void CBar::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0.0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i]-1];
	}


}

	// 新加 计算某个单元受到的重力

double CBar::Gravity()
	{
		CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);
		double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
			for (unsigned int i = 0; i < 3; i++)
			{
				DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
			}
		  L = sqrt(DX[0]*DX[0] + DX[1]*DX[1] + DX[2]*DX[2]);

		double G = 9.8*material_->Area * material_->rho_3D * L;
		return G;
}

void CBar::ElementPostInfo(double* stress, double* Displacement, double* PrePositions, double* PostPositions)
{
	double DX[3];
	for (int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
	}

	double Thetay[3],Thetaz[3];
	if (abs(DX[1]) > DBL_EPSILON)
	{
		Thetay[0] = 0.0;
		Thetay[1] = DX[2];
		Thetay[2] = -DX[1];
		Thetaz[0] = -DX[2] * DX[2] - DX[1] * DX[1];
		Thetaz[1] = DX[0] * DX[1];
		Thetaz[2] = DX[0] * DX[2];
	}
	else
	{
		Thetay[0] = DX[2];
		Thetay[1] = 0.0;
		Thetay[2] = -DX[0];
		Thetaz[0] = -DX[0] * DX[1];
		Thetaz[1] = DX[0] * DX[0] + DX[2] * DX[2];
		Thetaz[2] = -DX[1] * DX[2];
	}

	double leny = sqrt(Thetay[0] * Thetay[0] + Thetay[1] * Thetay[1] + Thetay[2] * Thetay[2]);
	Thetay[0] /= leny;
	Thetay[1] /= leny;
	Thetay[2] /= leny;
	double lenz = sqrt(Thetaz[0] * Thetaz[0] + Thetaz[1] * Thetaz[1] + Thetaz[2] * Thetaz[2]);
	Thetaz[0] /= lenz;
	Thetaz[1] /= lenz;
	Thetaz[2] /= lenz;

	// Determine the quasi-prepositions for POSTPROCESS
	// PrePositions = [X[0]  Y[0]  Z[0]  X[1]  Y[1] ... Z[7]]
	// Define the scale of co-dimension
	double magCodim = 1E-1;
	for (unsigned i = 0; i < 2; i++)
	{
		PrePositions[12 * i + 0] = nodes_[i]->XYZ[0] + magCodim * ( Thetay[0] + Thetaz[0]);
		PrePositions[12 * i + 1] = nodes_[i]->XYZ[1] + magCodim * ( Thetay[1] + Thetaz[1]);
		PrePositions[12 * i + 2] = nodes_[i]->XYZ[2] + magCodim * ( Thetay[2] + Thetaz[2]);

		PrePositions[12 * i + 3] = nodes_[i]->XYZ[0] + magCodim * ( - Thetay[0] + Thetaz[0]);
		PrePositions[12 * i + 4] = nodes_[i]->XYZ[1] + magCodim * ( - Thetay[1] + Thetaz[1]);
		PrePositions[12 * i + 5] = nodes_[i]->XYZ[2] + magCodim * ( - Thetay[2] + Thetaz[2]);

		PrePositions[12 * i + 6] = nodes_[i]->XYZ[0] + magCodim * ( - Thetay[0] - Thetaz[0]);
		PrePositions[12 * i + 7] = nodes_[i]->XYZ[1] + magCodim * ( - Thetay[1] - Thetaz[1]);
		PrePositions[12 * i + 8] = nodes_[i]->XYZ[2] + magCodim * ( - Thetay[2] - Thetaz[2]);

		PrePositions[12 * i + 9] = nodes_[i]->XYZ[0] + magCodim * ( Thetay[0] - Thetaz[0]);
		PrePositions[12 * i + 10] = nodes_[i]->XYZ[1] + magCodim * ( Thetay[1] - Thetaz[1]);
		PrePositions[12 * i + 11] = nodes_[i]->XYZ[2] + magCodim * ( Thetay[2] - Thetaz[2]);
	}
	
	// Calculate the postpositions for POSTPROCESS
	// PostPositions = [X[0]  Y[0]  Z[0]  X[1]  Y[1] ... Z[7]]
	double Disp[6];
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
		{
			Disp[i] = Displacement[LocationMatrix_[i] - 1];
		}
		else
		{
			Disp[i] = 0.0;
		}
	}
	for (unsigned i = 0; i < 2; i++)
	{
		for (unsigned j = 0; j < 4; j++)
		{
			PostPositions[12 * i + 3 * j] = PrePositions[12 * i + 3 * j] + Disp[3 * i];
			PostPositions[12 * i + 3 * j + 1] = PrePositions[12 * i + 3 * j + 1] + Disp[3 * i + 1];
			PostPositions[12 * i + 3 * j + 2] = PrePositions[12 * i + 3 * j + 2] + Disp[3 * i + 2];
		}
	}

	// Calculate element stress for POSTPROCESS
	CBarMaterial* material = dynamic_cast<CBarMaterial*>(ElementMaterial_);
	double L2 = DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2];
    double S[6];
	double ElementStress;
    for (unsigned int i = 0; i < 3; i++)
    {
        S[i] = -DX[i] * material->E / L2;
        S[i + 3] = -S[i];
    }
    ElementStress = 0.0;
    for (unsigned int i = 0; i < 6; i++)
    {
        if (LocationMatrix_[i])
			ElementStress += S[i] * Displacement[LocationMatrix_[i] - 1];
    }
	// Allocate the stress to each quasi-node
	// stress = [XX[0]  YY[0]  ZZ[0]  YZ[0]  ZX[0]  XY[0] ... ]
	clear(stress, 48);
	for (unsigned int i = 0; i < 8; i++)
	{
		stress[6 * i] = ElementStress;
	}
}
