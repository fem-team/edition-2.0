/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "3T.h"
#include "Bar.h"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
C3T::C3T()
{
	NEN_ = 3;	// Each element has 3 nodes
	nodes_ = new CNode*[NEN_];
    
    ND_ = 6;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
C3T::~C3T()
{
}

//	Read element data from stream Input
bool C3T::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2,N3;	// Left node number and right node number

	Input >> N1 >> N2 >>N3>> MSet;
    ElementMaterial_ = dynamic_cast<C3TMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];

	return true;
}

//	Write element data to stream
void C3T::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(9)<<nodes_[2]->NodeNumber << setw(9)
		   << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C3T::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 3 node 3T element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int C3T::SizeOfStiffnessMatrix() { return 21; }

//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C3T::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());
	C3TMaterial* material_ = dynamic_cast<C3TMaterial*>(ElementMaterial_);	// Pointer to material of the element
	

	Eigen::MatrixXd B(3,6);
	Eigen::MatrixXd KK(6,6);
	Eigen::MatrixXd DD(3,3);
	Eigen::MatrixXd BB(6,3);
	
	double E1;
	double poi;
	double D1;

	E1 = material_->E;
	poi = material_->Nu;
	D1 = E1 / (1 - poi * poi);

	DD(0, 0) = D1 * 1;
	DD(0, 1) = D1 * poi;
	DD(0, 2) = 0;
	DD(1, 0) = D1 * poi;
	DD(1, 1) = D1;
	DD(1, 2) = D1 * 0;
	DD(2, 0) = D1 * 0;
	DD(2, 1) = D1 * 0;
	DD(2, 2) = D1 * (1 - poi) / 2;

	double area ;
	
	double t = material_->thick;

	double DX[3];
	double DY[3];
	double DZ[3];

			DX[0] = (nodes_[2]->XYZ[0]) - (nodes_[1]->XYZ[0]);
			DX[1] = (nodes_[0]->XYZ[0]) - (nodes_[2]->XYZ[0]);
			DX[2] = (nodes_[1]->XYZ[0]) - (nodes_[0]->XYZ[0]);

			DY[0] = (nodes_[1]->XYZ[1]) - (nodes_[2]->XYZ[1]);
			DY[1] = (nodes_[2]->XYZ[1]) - (nodes_[0]->XYZ[1]);
			DY[2] = (nodes_[0]->XYZ[1]) - (nodes_[1]->XYZ[1]);
	
			DZ[0] = (nodes_[1]->XYZ[2]) - (nodes_[2]->XYZ[2]);
			DZ[1] = (nodes_[2]->XYZ[2]) - (nodes_[0]->XYZ[2]);
			DZ[2] = (nodes_[0]->XYZ[2]) - (nodes_[1]->XYZ[2]);
	
	double L12;
	double L22;
	double L32;
	double L1;
	double L2;
	double L3;
	double p;

	L12 = DX[0] * DX[0] + DY[0] * DY[0] + DZ[0] * DZ[0];
	L22 = DX[1] * DX[1] + DY[1] * DY[1] + DZ[1] * DZ[1];
	L32 = DX[2] * DX[2] + DY[2] * DY[2] + DZ[2] * DZ[2];
	L1 = sqrt(L12);
	L2 = sqrt(L22);
	L3 = sqrt(L32);
	p = (L1 + L2 + L3) / 2;
	area = sqrt( p * ( p - L1 )*( p - L2 )*( p - L3 ));
	
		B(0, 0) = DY[0];
		B(0, 1) = 0;
		B(0, 2) = DY[1];
		B(0, 3) = 0;
		B(0, 4) = DY[2];
		B(0, 5) = 0;
		B(1, 0) = 0;
		B(1, 1) = DX[0];
		B(1, 2) = 0;
		B(1, 3) = DX[1];
		B(1, 4) = 0;
		B(1, 5) = DX[2];
		B(2, 0) = DX[0];
		B(2, 1) = DY[0];
	    B(2, 2) = DX[1];
		B(2, 3) = DY[1];
		B(2, 4) = DX[2];
		B(2, 5) = DY[2];
		
		B =B / (2*area);
		BB = B.transpose();
			
		KK = t * area * BB * DD * B ;
		Matrix[0] = KK(0, 0);
		Matrix[1] = KK(1, 1);
		Matrix[2] = KK(1, 0);
		Matrix[3] = KK(2, 2);
		Matrix[4] = KK(2, 1);
		Matrix[5] = KK(2, 0);
		Matrix[6] = KK(3, 3);
		Matrix[7] = KK(3, 2);
		Matrix[8] = KK(3, 1);
		Matrix[9] = KK(3, 0);
		Matrix[10] = KK(4, 4);
		Matrix[11] = KK(4, 3);
		Matrix[12] = KK(4, 2);
		Matrix[13] = KK(4, 1);
		Matrix[14] = KK(4, 0);
		Matrix[15] = KK(5, 5);
		Matrix[16] = KK(5, 4);
		Matrix[17] = KK(5, 3);
		Matrix[18] = KK(5, 2);
		Matrix[19] = KK(5, 1);
		Matrix[20] = KK(5, 0);
	}
	
	
	//	应变 Calculate element stress 
	void C3T::ElementStress(double* stress, double* Displacement)
	{
		C3TMaterial* material_ = dynamic_cast<C3TMaterial*>(ElementMaterial_);	// Pointer to material of the element
		//	Calculate D1
		double D1 = material_->E / (1 - material_->Nu*material_->Nu);
		
	
		double poi = material_->Nu;
		

			Eigen::MatrixXd  stress1(3,1);
			Eigen::MatrixXd  DD(3,3);
			Eigen::MatrixXd  strain(3,1);     //计算应变strain
			Eigen::MatrixXd  d(1,6);
			Eigen::MatrixXd  B(3,6);
			Eigen::MatrixXd  dz(6,1);
			Eigen::MatrixXd  stress2(3,1);
			DD(0, 0) = D1 * 1 ;
			DD(0, 1) = D1 * poi ;
			DD(0, 2) = 0 ;
			DD(1, 0) = D1 * poi ;
			DD(1, 1) = D1  ;
			DD(1, 2) = D1 * 0 ;
			DD(2, 0) = D1 * 0 ;
			DD(2, 1) = D1 * 0 ;
			DD(2, 2) = D1 * (1-poi)/2;

			for (unsigned int i = 0; i < 6; i++)
			{

				if ( LocationMatrix_[i])
					d(0,i) = Displacement[ LocationMatrix_[i] - 1];
				else if ( LocationMatrix_[i] == 0)
					d(0,i) = 0;

			}

			double DX[3];
			double DY[3];
			double DZ[3];

			DX[0] = nodes_[2]->XYZ[0] - nodes_[1]->XYZ[0];
			DX[1] = nodes_[0]->XYZ[0] - nodes_[2]->XYZ[0];
			DX[2] = nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0];

			DY[0] = nodes_[1]->XYZ[1] - nodes_[2]->XYZ[1];
			DY[1] = nodes_[2]->XYZ[1] - nodes_[0]->XYZ[1];
			DY[2] = nodes_[0]->XYZ[1] - nodes_[1]->XYZ[1];

			DZ[0] = nodes_[1]->XYZ[2] - nodes_[2]->XYZ[2];
			DZ[1] = nodes_[2]->XYZ[2] - nodes_[0]->XYZ[2];
			DZ[2] = nodes_[0]->XYZ[2] - nodes_[1]->XYZ[2];

			double L12;
			double L22;
			double L32;
			double L1;
			double L2;
			double L3;
			double p;

			L12 = DX[0] * DX[0] + DY[0] * DY[0] + DZ[0] * DZ[0];
			L22 = DX[1] * DX[1] + DY[1] * DY[1] + DZ[1] * DZ[1];
			L32 = DX[2] * DX[2] + DY[2] * DY[2] + DZ[2] * DZ[2];
			L1 = sqrt(L12);
			L2 = sqrt(L22);
			L3 = sqrt(L32);
			p = (L1 + L2 + L3) / 2;

			double area = sqrt(p*(p - L1)*(p - L2)*(p - L3));

			B(0, 0) = DY[0];
			B(0, 1) = 0;
			B(0, 2) = DY[1];
			B(0, 3) = 0;
			B(0, 4) = DY[2];
			B(0, 5) = 0;
			B(1, 0) = 0;
			B(1, 1) = DX[0];
			B(1, 2) = 0;
			B(1, 3) = DX[1];
			B(1, 4) = 0;
			B(1, 5) = DX[2];
			B(2, 0) = DX[0];
			B(2, 1) = DY[0];
			B(2, 2) = DX[1];
			B(2, 3) = DY[1];
			B(2, 4) = DX[2];
			B(2, 5) = DY[2];

			dz = d.transpose();
			strain = (B * dz)/(2*area);
			stress1 = DD * strain;
			double X0 = nodes_[0]->XYZ[0];
			double X1 = nodes_[1]->XYZ[0];
			double X2 = nodes_[2]->XYZ[0];
			double Y0 = nodes_[0]->XYZ[1];
			double Y1 = nodes_[1]->XYZ[1];
			double Y2 = nodes_[2]->XYZ[1];
			double Z0 = nodes_[0]->XYZ[0];
			double Z1 = nodes_[1]->XYZ[1];
			double Z2 = nodes_[2]->XYZ[2];

			stress[0] = X0;
			stress[1] = Y0;
			stress[2] = Z0;
			stress[3] = stress1(0, 0);
			stress[4] = stress1(1, 0);
			stress[5] = stress1(2,0);
			stress[6] = X1;
			stress[7] = Y1;
			stress[8] = Z1;
			stress[9] = stress1(0, 0);
			stress[10] = stress1(1, 0);
			stress[11] = stress1(2, 0);
			stress[12] = X2;
			stress[13] = Y2;
			stress[14] = Z2;
			stress[15] = stress1(0, 0);
			stress[16] = stress1(1, 0);
			stress[17] = stress1(2, 0);
			
		}

	

void C3T::ElementPostInfo(double* stress, double* Displacement, double* PrePositions,
                                double* PostPositions)
{
    ElementStress(stress, Displacement);
    for (unsigned index = 0; index < 9; ++index)
    {
        if (LocationMatrix_[index])
        {
            PrePositions[index] = nodes_[index / 3]->XYZ[index % 3];
            PostPositions[index] = PrePositions[index] + Displacement[LocationMatrix_[index] - 1];
        }
        else
        {
            PrePositions[index] = PostPositions[index] = nodes_[index / 3]->XYZ[index % 3];
        }
    }
}

//Gravity
double C3T::Gravity()
{
	return 0;
}
