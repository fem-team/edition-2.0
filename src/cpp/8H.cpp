/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include "8H.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
C8H::C8H()
{
	NEN_ = 8;
	nodes_ = new CNode*[NEN_];
    
    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
C8H::~C8H()
{
	delete [] nodes_;
    delete [] LocationMatrix_;
}
//	Read element data from stream Input
bool C8H::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int MSet;
	unsigned int N1, N2, N3 ,N4 ,N5 ,N6 ,N7 ,N8;
	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial_ = dynamic_cast<C8HMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];
	return true;
}
//	Write element data to stream
void C8H::Write(COutputter& output, unsigned int Ele)
{
  output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber 
		 << setw(9) << nodes_[1]->NodeNumber 
		 << setw(9) << nodes_[2]->NodeNumber 
		 << setw(9) << nodes_[3]->NodeNumber
		 << setw(9) << nodes_[4]->NodeNumber
		 << setw(9) << nodes_[5]->NodeNumber
		 << setw(9) << nodes_[6]->NodeNumber
		 << setw(9) << nodes_[7]->NodeNumber
		 << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C8H::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For Plate element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int C8H::SizeOfStiffnessMatrix() { return 300; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element

void C8H::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	// Construct constitutive matrix
	C8HMaterial* material = dynamic_cast<C8HMaterial*>(ElementMaterial_);	// Pointer to material of the element
//	double v = material->nu;
//	double k = material->E * (1 - v) / (1 + v) / (1 - 2 * v);
//	double D[3];
//	D[0] = k;
//	D[1] = k * v / (1 - v);
//	D[2] = k * (1 - 2 * v) / 2.0 / (1 - v);
	double E = material->E;
	double nv = material->Nu;
	double G = 0.5*E/(1+nv);//D[2]
	double lambda = nv*E/((1.0 +nv)*(1.0 -2.0 *nv));//D[1]
	double mu=lambda+2.0 *G;//mu
	// Construct coordinate matrix
	double COORXYZ[24];
	for (unsigned int i = 0; i < 8; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			COORXYZ[3 * i + j] = nodes_[i]->XYZ[j];
		}
	}

	// construct Jacobi matrix
	const double xi8[8]   = { 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 };
	const double eta8[8]  = {-0.577350269189626 , 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,0.577350269189626 , 0.577350269189626 ,-0.577350269189626 };
	const double zeta8[8] = {-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 , 0.577350269189626 , 0.577350269189626 };

	for (unsigned p = 0; p < 8; p++)
	{
		double xi = xi8[p];
		double eta = eta8[p];
		double zeta = zeta8[p];

		double GN[12];
		GN[0] = (1.0-eta)*(1.0-zeta) / 8.0;
		GN[1] = (1.0+eta)*(1.0-zeta) / 8.0;
		GN[2] = (1.0-eta)*(1.0+zeta) / 8.0;
		GN[3] = (1.0+eta)*(1.0+zeta) / 8.0;
		GN[4] = (1.0+xi)*(1.0-zeta) / 8.0;
		GN[5] = (1.0-xi)*(1.0-zeta) / 8.0;
		GN[6] = (1.0+xi)*(1.0+zeta) / 8.0;
		GN[7] = (1.0-xi)*(1.0+zeta) / 8.0;
		GN[8] = (1.0+xi)*(1.0-eta) / 8.0;
		GN[9] = (1.0+xi)*(1.0+eta) / 8.0;
		GN[10] = (1.0-xi)*(1.0+eta) / 8.0;
		GN[11] = (1.0-xi)*(1.0-eta) / 8.0;

		double J[9];
		J[0] = COORXYZ[0]*GN[0]+COORXYZ[3]*GN[1]-COORXYZ[6]*GN[1]-COORXYZ[9]*GN[0]+COORXYZ[12]*GN[2]+COORXYZ[15]*GN[3]-COORXYZ[18]*GN[3]-COORXYZ[21]*GN[2];
		J[1] = -COORXYZ[0]*GN[4]+COORXYZ[3]*GN[4]+COORXYZ[6]*GN[5]-COORXYZ[9]*GN[5]-COORXYZ[12]*GN[6]+COORXYZ[15]*GN[6]+COORXYZ[18]*GN[7]-COORXYZ[21]*GN[7];
		J[2] = -COORXYZ[0]*GN[8]-COORXYZ[3]*GN[9]-COORXYZ[6]*GN[10]-COORXYZ[9]*GN[11]+COORXYZ[12]*GN[8]+COORXYZ[15]*GN[9]+COORXYZ[18]*GN[10]+COORXYZ[21]*GN[11];
		J[3] = COORXYZ[1]*GN[0]+COORXYZ[4]*GN[1]-COORXYZ[7]*GN[1]-COORXYZ[10]*GN[0]+COORXYZ[13]*GN[2]+COORXYZ[16]*GN[3]-COORXYZ[19]*GN[3]-COORXYZ[22]*GN[2];
		J[4] = -COORXYZ[1]*GN[4]+COORXYZ[4]*GN[4]+COORXYZ[7]*GN[5]-COORXYZ[10]*GN[5]-COORXYZ[13]*GN[6]+COORXYZ[16]*GN[6]+COORXYZ[19]*GN[7]-COORXYZ[22]*GN[7];
		J[5] = -COORXYZ[1]*GN[8]-COORXYZ[4]*GN[9]-COORXYZ[7]*GN[10]-COORXYZ[10]*GN[11]+COORXYZ[13]*GN[8]+COORXYZ[16]*GN[9]+COORXYZ[19]*GN[10]+COORXYZ[22]*GN[11];
		J[6] = COORXYZ[2]*GN[0]+COORXYZ[5]*GN[1]-COORXYZ[8]*GN[1]-COORXYZ[11]*GN[0]+COORXYZ[14]*GN[2]+COORXYZ[17]*GN[3]-COORXYZ[20]*GN[3]-COORXYZ[23]*GN[2];
		J[7] = -COORXYZ[2]*GN[4]+COORXYZ[5]*GN[4]+COORXYZ[8]*GN[5]-COORXYZ[11]*GN[5]-COORXYZ[14]*GN[6]+COORXYZ[17]*GN[6]+COORXYZ[20]*GN[7]-COORXYZ[23]*GN[7];
		J[8] = -COORXYZ[2]*GN[8]-COORXYZ[5]*GN[9]-COORXYZ[8]*GN[10]-COORXYZ[11]*GN[11]+COORXYZ[14]*GN[8]+COORXYZ[17]*GN[9]+COORXYZ[20]*GN[10]+COORXYZ[23]*GN[11];

		double detJ = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];

		double InvJ[9];
		InvJ[0] = (J[4]*J[8]-J[5]*J[7])/detJ;
		InvJ[1] = -(J[1]*J[8]-J[2]*J[7])/detJ;
		InvJ[2] = (J[1]*J[5]-J[2]*J[4])/detJ;
		InvJ[3] = -(J[3]*J[8]-J[5]*J[6])/detJ;
		InvJ[4] = (J[0]*J[8]-J[2]*J[6])/detJ;
		InvJ[5] = -(J[0]*J[5]-J[2]*J[3])/detJ;
		InvJ[6] = (J[3]*J[7]-J[4]*J[6])/detJ;
		InvJ[7] = -(J[0]*J[7]-J[1]*J[6])/detJ;
		InvJ[8] = (J[0]*J[4]-J[1]*J[3])/detJ;

		double B1[24];
		B1[0] = GN[0] * InvJ[0] - GN[4] * InvJ[3] - GN[8] * InvJ[6];
		B1[1] = GN[0] * InvJ[1] - GN[4] * InvJ[4] - GN[8] * InvJ[7];
		B1[2] = GN[0] * InvJ[2] - GN[4] * InvJ[5] - GN[8] * InvJ[8];
		B1[3] = GN[1] * InvJ[0] + GN[4] * InvJ[3] - GN[9] * InvJ[6];
		B1[4] = GN[1] * InvJ[1] + GN[4] * InvJ[4] - GN[9] * InvJ[7];
		B1[5] = GN[1] * InvJ[2] + GN[4] * InvJ[5] - GN[9] * InvJ[8];
		B1[6] = -GN[1] * InvJ[0] + GN[5] * InvJ[3] - GN[10] * InvJ[6];
		B1[7] = -GN[1] * InvJ[1] + GN[5] * InvJ[4] - GN[10] * InvJ[7];
		B1[8] = -GN[1] * InvJ[2] + GN[5] * InvJ[5] - GN[10] * InvJ[8];
		B1[9] = -GN[0] * InvJ[0] - GN[5] * InvJ[3] - GN[11] * InvJ[6];
		B1[10] = -GN[0] * InvJ[1] - GN[5] * InvJ[4] - GN[11] * InvJ[7];
		B1[11] = -GN[0] * InvJ[2] - GN[5] * InvJ[5] - GN[11] * InvJ[8];
		B1[12] = GN[2] * InvJ[0] - GN[6] * InvJ[3] + GN[8] * InvJ[6];
		B1[13] = GN[2] * InvJ[1] - GN[6] * InvJ[4] + GN[8] * InvJ[7];
		B1[14] = GN[2] * InvJ[2] - GN[6] * InvJ[5] + GN[8] * InvJ[8];
		B1[15] = GN[3] * InvJ[0] + GN[6] * InvJ[3] + GN[9] * InvJ[6];
		B1[16] = GN[3] * InvJ[1] + GN[6] * InvJ[4] + GN[9] * InvJ[7];
		B1[17] = GN[3] * InvJ[2] + GN[6] * InvJ[5] + GN[9] * InvJ[8];
		B1[18] = -GN[3] * InvJ[0] + GN[7] * InvJ[3] + GN[10] * InvJ[6];
		B1[19] = -GN[3] * InvJ[1] + GN[7] * InvJ[4] + GN[10] * InvJ[7];
		B1[20] = -GN[3] * InvJ[2] + GN[7] * InvJ[5] + GN[10] * InvJ[8];
		B1[21] = -GN[2] * InvJ[0] - GN[7] * InvJ[3] + GN[11] * InvJ[6];
		B1[22] = -GN[2] * InvJ[1] - GN[7] * InvJ[4] + GN[11] * InvJ[7];
		B1[23] = -GN[2] * InvJ[2] - GN[7] * InvJ[5] + GN[11] * InvJ[8];

		// construct element stiffness matrix
		Matrix[0] += detJ*(mu * (B1[0] * B1[0]) + G * (B1[1] * B1[1]) + G * (B1[2] * B1[2]));
		Matrix[1] += detJ*(mu * (B1[1] * B1[1]) + G * (B1[0] * B1[0]) + G * (B1[2] * B1[2]));
		Matrix[2] += detJ*(lambda * B1[0] * B1[1] + G * B1[0] * B1[1]);
		Matrix[3] += detJ*(mu * (B1[2] * B1[2]) + G * (B1[0] * B1[0]) + G * (B1[1] * B1[1]));
		Matrix[4] += detJ*(lambda * B1[1] * B1[2] + G * B1[1] * B1[2]);
		Matrix[5] += detJ*(lambda * B1[0] * B1[2] + G * B1[0] * B1[2]);
		Matrix[6] += detJ*(mu * (B1[3] * B1[3]) + G * (B1[4] * B1[4]) + G * (B1[5] * B1[5]));
		Matrix[7] += detJ*(lambda * B1[2] * B1[3] + G * B1[0] * B1[5]);
		Matrix[8] += detJ*(lambda * B1[1] * B1[3] + G * B1[0] * B1[4]);
		Matrix[9] += detJ*(mu * B1[0] * B1[3] + G * B1[1] * B1[4] + G * B1[2] * B1[5]);
		Matrix[10] += detJ*(mu * (B1[4] * B1[4]) + G * (B1[3] * B1[3]) + G * (B1[5] * B1[5]));
		Matrix[11] += detJ*(lambda * B1[3] * B1[4] + G * B1[3] * B1[4]);
		Matrix[12] += detJ*(lambda * B1[2] * B1[4] + G * B1[1] * B1[5]);
		Matrix[13] += detJ*(mu * B1[1] * B1[4] + G * B1[0] * B1[3] + G * B1[2] * B1[5]);
		Matrix[14] += detJ*(lambda * B1[0] * B1[4] + G * B1[1] * B1[3]);
		Matrix[15] += detJ*(mu * (B1[5] * B1[5]) + G * (B1[3] * B1[3]) + G * (B1[4] * B1[4]));
		Matrix[16] += detJ*(lambda * B1[4] * B1[5] + G * B1[4] * B1[5]);
		Matrix[17] += detJ*(lambda * B1[3] * B1[5] + G * B1[3] * B1[5]);
		Matrix[18] += detJ*(G * B1[0] * B1[3] + mu * B1[2] * B1[5] + G * B1[1] * B1[4]);
		Matrix[19] += detJ*(lambda * B1[1] * B1[5] + G * B1[2] * B1[4]);
		Matrix[20] += detJ*(lambda * B1[0] * B1[5] + G * B1[2] * B1[3]);
		Matrix[21] += detJ*(mu * (B1[6] * B1[6]) + G * (B1[7] * B1[7]) + G * (B1[8] * B1[8]));
		Matrix[22] += detJ*(lambda * B1[5] * B1[6] + G * B1[3] * B1[8]);
		Matrix[23] += detJ*(lambda * B1[4] * B1[6] + G * B1[3] * B1[7]);
		Matrix[24] += detJ*(mu * B1[3] * B1[6] + G * B1[4] * B1[7] + G * B1[5] * B1[8]);
		Matrix[25] += detJ*(lambda * B1[2] * B1[6] + G * B1[0] * B1[8]);
		Matrix[26] += detJ*(lambda * B1[1] * B1[6] + G * B1[0] * B1[7]);
		Matrix[27] += detJ*(mu * B1[0] * B1[6] + G * B1[1] * B1[7] + G * B1[2] * B1[8]);
		Matrix[28] += detJ*(mu * (B1[7] * B1[7]) + G * (B1[6] * B1[6]) + G * (B1[8] * B1[8]));
		Matrix[29] += detJ*(lambda * B1[6] * B1[7] + G * B1[6] * B1[7]);
		Matrix[30] += detJ*(lambda * B1[5] * B1[7] + G * B1[4] * B1[8]);
		Matrix[31] += detJ*(mu * B1[4] * B1[7] + G * B1[3] * B1[6] + G * B1[5] * B1[8]);
		Matrix[32] += detJ*(lambda * B1[3] * B1[7] + G * B1[4] * B1[6]);
		Matrix[33] += detJ*(lambda * B1[2] * B1[7] + G * B1[1] * B1[8]);
		Matrix[34] += detJ*(mu * B1[1] * B1[7] + G * B1[0] * B1[6] + G * B1[2] * B1[8]);
		Matrix[35] += detJ*(lambda * B1[0] * B1[7] + G * B1[1] * B1[6]);
		Matrix[36] += detJ*(mu * (B1[8] * B1[8]) + G * (B1[6] * B1[6]) + G * (B1[7] * B1[7]));
		Matrix[37] += detJ*(lambda * B1[7] * B1[8] + G * B1[7] * B1[8]);
		Matrix[38] += detJ*(lambda * B1[6] * B1[8] + G * B1[6] * B1[8]);
		Matrix[39] += detJ*(G * B1[3] * B1[6] + mu * B1[5] * B1[8] + G * B1[4] * B1[7]);
		Matrix[40] += detJ*(lambda * B1[4] * B1[8] + G * B1[5] * B1[7]);
		Matrix[41] += detJ*(lambda * B1[3] * B1[8] + G * B1[5] * B1[6]);
		Matrix[42] += detJ*(G * B1[0] * B1[6] + mu * B1[2] * B1[8] + G * B1[1] * B1[7]);
		Matrix[43] += detJ*(lambda * B1[1] * B1[8] + G * B1[2] * B1[7]);
		Matrix[44] += detJ*(lambda * B1[0] * B1[8] + G * B1[2] * B1[6]);
		Matrix[45] += detJ*(mu * (B1[9] * B1[9]) + G * (B1[10] * B1[10]) + G * (B1[11] * B1[11]));
		Matrix[46] += detJ*(lambda * B1[8] * B1[9] + G * B1[6] * B1[11]);
		Matrix[47] += detJ*(lambda * B1[7] * B1[9] + G * B1[6] * B1[10]);
		Matrix[48] += detJ*(mu * B1[6] * B1[9] + G * B1[7] * B1[10] + G * B1[8] * B1[11]);
		Matrix[49] += detJ*(lambda * B1[5] * B1[9] + G * B1[3] * B1[11]);
		Matrix[50] += detJ*(lambda * B1[4] * B1[9] + G * B1[3] * B1[10]);
		Matrix[51] += detJ*(mu * B1[3] * B1[9] + G * B1[4] * B1[10] + G * B1[5] * B1[11]);
		Matrix[52] += detJ*(lambda * B1[2] * B1[9] + G * B1[0] * B1[11]);
		Matrix[53] += detJ*(lambda * B1[1] * B1[9] + G * B1[0] * B1[10]);
		Matrix[54] += detJ*(mu * B1[0] * B1[9] + G * B1[1] * B1[10] + G * B1[2] * B1[11]);
		Matrix[55] += detJ*(mu * (B1[10] * B1[10]) + G * (B1[9] * B1[9]) + G * (B1[11] * B1[11]));
		Matrix[56] += detJ*(lambda * B1[9] * B1[10] + G * B1[9] * B1[10]);
		Matrix[57] += detJ*(lambda * B1[8] * B1[10] + G * B1[7] * B1[11]);
		Matrix[58] += detJ*(mu * B1[7] * B1[10] + G * B1[6] * B1[9] + G * B1[8] * B1[11]);
		Matrix[59] += detJ*(lambda * B1[6] * B1[10] + G * B1[7] * B1[9]);
		Matrix[60] += detJ*(lambda * B1[5] * B1[10] + G * B1[4] * B1[11]);
		Matrix[61] += detJ*(mu * B1[4] * B1[10] + G * B1[3] * B1[9] + G * B1[5] * B1[11]);
		Matrix[62] += detJ*(lambda * B1[3] * B1[10] + G * B1[4] * B1[9]);
		Matrix[63] += detJ*(lambda * B1[2] * B1[10] + G * B1[1] * B1[11]);
		Matrix[64] += detJ*(mu * B1[1] * B1[10] + G * B1[0] * B1[9] + G * B1[2] * B1[11]);
		Matrix[65] += detJ*(lambda * B1[0] * B1[10] + G * B1[1] * B1[9]);
		Matrix[66] += detJ*(mu * (B1[11] * B1[11]) + G * (B1[9] * B1[9]) + G * (B1[10] * B1[10]));
		Matrix[67] += detJ*(lambda * B1[10] * B1[11] + G * B1[10] * B1[11]);
		Matrix[68] += detJ*(lambda * B1[9] * B1[11] + G * B1[9] * B1[11]);
		Matrix[69] += detJ*(G * B1[6] * B1[9] + mu * B1[8] * B1[11] + G * B1[7] * B1[10]);
		Matrix[70] += detJ*(lambda * B1[7] * B1[11] + G * B1[8] * B1[10]);
		Matrix[71] += detJ*(lambda * B1[6] * B1[11] + G * B1[8] * B1[9]);
		Matrix[72] += detJ*(G * B1[3] * B1[9] + mu * B1[5] * B1[11] + G * B1[4] * B1[10]);
		Matrix[73] += detJ*(lambda * B1[4] * B1[11] + G * B1[5] * B1[10]);
		Matrix[74] += detJ*(lambda * B1[3] * B1[11] + G * B1[5] * B1[9]);
		Matrix[75] += detJ*(G * B1[0] * B1[9] + mu * B1[2] * B1[11] + G * B1[1] * B1[10]);
		Matrix[76] += detJ*(lambda * B1[1] * B1[11] + G * B1[2] * B1[10]);
		Matrix[77] += detJ*(lambda * B1[0] * B1[11] + G * B1[2] * B1[9]);
		Matrix[78] += detJ*(mu * (B1[12] * B1[12]) + G * (B1[13] * B1[13]) + G * (B1[14] * B1[14]));
		Matrix[79] += detJ*(lambda * B1[11] * B1[12] + G * B1[9] * B1[14]);
		Matrix[80] += detJ*(lambda * B1[10] * B1[12] + G * B1[9] * B1[13]);
		Matrix[81] += detJ*(mu * B1[9] * B1[12] + G * B1[10] * B1[13] + G * B1[11] * B1[14]);
		Matrix[82] += detJ*(lambda * B1[8] * B1[12] + G * B1[6] * B1[14]);
		Matrix[83] += detJ*(lambda * B1[7] * B1[12] + G * B1[6] * B1[13]);
		Matrix[84] += detJ*(mu * B1[6] * B1[12] + G * B1[7] * B1[13] + G * B1[8] * B1[14]);
		Matrix[85] += detJ*(lambda * B1[5] * B1[12] + G * B1[3] * B1[14]);
		Matrix[86] += detJ*(lambda * B1[4] * B1[12] + G * B1[3] * B1[13]);
		Matrix[87] += detJ*(mu * B1[3] * B1[12] + G * B1[4] * B1[13] + G * B1[5] * B1[14]);
		Matrix[88] += detJ*(lambda * B1[2] * B1[12] + G * B1[0] * B1[14]);
		Matrix[89] += detJ*(lambda * B1[1] * B1[12] + G * B1[0] * B1[13]);
		Matrix[90] += detJ*(mu * B1[0] * B1[12] + G * B1[1] * B1[13] + G * B1[2] * B1[14]);
		Matrix[91] += detJ*(mu * (B1[13] * B1[13]) + G * (B1[12] * B1[12]) + G * (B1[14] * B1[14]));
		Matrix[92] += detJ*(lambda * B1[12] * B1[13] + G * B1[12] * B1[13]);
		Matrix[93] += detJ*(lambda * B1[11] * B1[13] + G * B1[10] * B1[14]);
		Matrix[94] += detJ*(mu * B1[10] * B1[13] + G * B1[9] * B1[12] + G * B1[11] * B1[14]);
		Matrix[95] += detJ*(lambda * B1[9] * B1[13] + G * B1[10] * B1[12]);
		Matrix[96] += detJ*(lambda * B1[8] * B1[13] + G * B1[7] * B1[14]);
		Matrix[97] += detJ*(mu * B1[7] * B1[13] + G * B1[6] * B1[12] + G * B1[8] * B1[14]);
		Matrix[98] += detJ*(lambda * B1[6] * B1[13] + G * B1[7] * B1[12]);
		Matrix[99] += detJ*(lambda * B1[5] * B1[13] + G * B1[4] * B1[14]);
		Matrix[100] += detJ*(mu * B1[4] * B1[13] + G * B1[3] * B1[12] + G * B1[5] * B1[14]);
		Matrix[101] += detJ*(lambda * B1[3] * B1[13] + G * B1[4] * B1[12]);
		Matrix[102] += detJ*(lambda * B1[2] * B1[13] + G * B1[1] * B1[14]);
		Matrix[103] += detJ*(mu * B1[1] * B1[13] + G * B1[0] * B1[12] + G * B1[2] * B1[14]);
		Matrix[104] += detJ*(lambda * B1[0] * B1[13] + G * B1[1] * B1[12]);
		Matrix[105] += detJ*(mu * (B1[14] * B1[14]) + G * (B1[12] * B1[12]) + G * (B1[13] * B1[13]));
		Matrix[106] += detJ*(lambda * B1[13] * B1[14] + G * B1[13] * B1[14]);
		Matrix[107] += detJ*(lambda * B1[12] * B1[14] + G * B1[12] * B1[14]);
		Matrix[108] += detJ*(G * B1[9] * B1[12] + mu * B1[11] * B1[14] + G * B1[10] * B1[13]);
		Matrix[109] += detJ*(lambda * B1[10] * B1[14] + G * B1[11] * B1[13]);
		Matrix[110] += detJ*(lambda * B1[9] * B1[14] + G * B1[11] * B1[12]);
		Matrix[111] += detJ*(G * B1[6] * B1[12] + mu * B1[8] * B1[14] + G * B1[7] * B1[13]);
		Matrix[112] += detJ*(lambda * B1[7] * B1[14] + G * B1[8] * B1[13]);
		Matrix[113] += detJ*(lambda * B1[6] * B1[14] + G * B1[8] * B1[12]);
		Matrix[114] += detJ*(G * B1[3] * B1[12] + mu * B1[5] * B1[14] + G * B1[4] * B1[13]);
		Matrix[115] += detJ*(lambda * B1[4] * B1[14] + G * B1[5] * B1[13]);
		Matrix[116] += detJ*(lambda * B1[3] * B1[14] + G * B1[5] * B1[12]);
		Matrix[117] += detJ*(G * B1[0] * B1[12] + mu * B1[2] * B1[14] + G * B1[1] * B1[13]);
		Matrix[118] += detJ*(lambda * B1[1] * B1[14] + G * B1[2] * B1[13]);
		Matrix[119] += detJ*(lambda * B1[0] * B1[14] + G * B1[2] * B1[12]);
		Matrix[120] += detJ*(mu * (B1[15] * B1[15]) + G * (B1[16] * B1[16]) + G * (B1[17] * B1[17]));
		Matrix[121] += detJ*(lambda * B1[14] * B1[15] + G * B1[12] * B1[17]);
		Matrix[122] += detJ*(lambda * B1[13] * B1[15] + G * B1[12] * B1[16]);
		Matrix[123] += detJ*(mu * B1[12] * B1[15] + G * B1[13] * B1[16] + G * B1[14] * B1[17]);
		Matrix[124] += detJ*(lambda * B1[11] * B1[15] + G * B1[9] * B1[17]);
		Matrix[125] += detJ*(lambda * B1[10] * B1[15] + G * B1[9] * B1[16]);
		Matrix[126] += detJ*(mu * B1[9] * B1[15] + G * B1[10] * B1[16] + G * B1[11] * B1[17]);
		Matrix[127] += detJ*(lambda * B1[8] * B1[15] + G * B1[6] * B1[17]);
		Matrix[128] += detJ*(lambda * B1[7] * B1[15] + G * B1[6] * B1[16]);
		Matrix[129] += detJ*(mu * B1[6] * B1[15] + G * B1[7] * B1[16] + G * B1[8] * B1[17]);
		Matrix[130] += detJ*(lambda * B1[5] * B1[15] + G * B1[3] * B1[17]);
		Matrix[131] += detJ*(lambda * B1[4] * B1[15] + G * B1[3] * B1[16]);
		Matrix[132] += detJ*(mu * B1[3] * B1[15] + G * B1[4] * B1[16] + G * B1[5] * B1[17]);
		Matrix[133] += detJ*(lambda * B1[2] * B1[15] + G * B1[0] * B1[17]);
		Matrix[134] += detJ*(lambda * B1[1] * B1[15] + G * B1[0] * B1[16]);
		Matrix[135] += detJ*(mu * B1[0] * B1[15] + G * B1[1] * B1[16] + G * B1[2] * B1[17]);
		Matrix[136] += detJ*(mu * (B1[16] * B1[16]) + G * (B1[15] * B1[15]) + G * (B1[17] * B1[17]));
		Matrix[137] += detJ*(lambda * B1[15] * B1[16] + G * B1[15] * B1[16]);
		Matrix[138] += detJ*(lambda * B1[14] * B1[16] + G * B1[13] * B1[17]);
		Matrix[139] += detJ*(mu * B1[13] * B1[16] + G * B1[12] * B1[15] + G * B1[14] * B1[17]);
		Matrix[140] += detJ*(lambda * B1[12] * B1[16] + G * B1[13] * B1[15]);
		Matrix[141] += detJ*(lambda * B1[11] * B1[16] + G * B1[10] * B1[17]);
		Matrix[142] += detJ*(mu * B1[10] * B1[16] + G * B1[9] * B1[15] + G * B1[11] * B1[17]);
		Matrix[143] += detJ*(lambda * B1[9] * B1[16] + G * B1[10] * B1[15]);
		Matrix[144] += detJ*(lambda * B1[8] * B1[16] + G * B1[7] * B1[17]);
		Matrix[145] += detJ*(mu * B1[7] * B1[16] + G * B1[6] * B1[15] + G * B1[8] * B1[17]);
		Matrix[146] += detJ*(lambda * B1[6] * B1[16] + G * B1[7] * B1[15]);
		Matrix[147] += detJ*(lambda * B1[5] * B1[16] + G * B1[4] * B1[17]);
		Matrix[148] += detJ*(mu * B1[4] * B1[16] + G * B1[3] * B1[15] + G * B1[5] * B1[17]);
		Matrix[149] += detJ*(lambda * B1[3] * B1[16] + G * B1[4] * B1[15]);
		Matrix[150] += detJ*(lambda * B1[2] * B1[16] + G * B1[1] * B1[17]);
		Matrix[151] += detJ*(mu * B1[1] * B1[16] + G * B1[0] * B1[15] + G * B1[2] * B1[17]);
		Matrix[152] += detJ*(lambda * B1[0] * B1[16] + G * B1[1] * B1[15]);
		Matrix[153] += detJ*(mu * (B1[17] * B1[17]) + G * (B1[15] * B1[15]) + G * (B1[16] * B1[16]));
		Matrix[154] += detJ*(lambda * B1[16] * B1[17] + G * B1[16] * B1[17]);
		Matrix[155] += detJ*(lambda * B1[15] * B1[17] + G * B1[15] * B1[17]);
		Matrix[156] += detJ*(G * B1[12] * B1[15] + mu * B1[14] * B1[17] + G * B1[13] * B1[16]);
		Matrix[157] += detJ*(lambda * B1[13] * B1[17] + G * B1[14] * B1[16]);
		Matrix[158] += detJ*(lambda * B1[12] * B1[17] + G * B1[14] * B1[15]);
		Matrix[159] += detJ*(G * B1[9] * B1[15] + mu * B1[11] * B1[17] + G * B1[10] * B1[16]);
		Matrix[160] += detJ*(lambda * B1[10] * B1[17] + G * B1[11] * B1[16]);
		Matrix[161] += detJ*(lambda * B1[9] * B1[17] + G * B1[11] * B1[15]);
		Matrix[162] += detJ*(G * B1[6] * B1[15] + mu * B1[8] * B1[17] + G * B1[7] * B1[16]);
		Matrix[163] += detJ*(lambda * B1[7] * B1[17] + G * B1[8] * B1[16]);
		Matrix[164] += detJ*(lambda * B1[6] * B1[17] + G * B1[8] * B1[15]);
		Matrix[165] += detJ*(G * B1[3] * B1[15] + mu * B1[5] * B1[17] + G * B1[4] * B1[16]);
		Matrix[166] += detJ*(lambda * B1[4] * B1[17] + G * B1[5] * B1[16]);
		Matrix[167] += detJ*(lambda * B1[3] * B1[17] + G * B1[5] * B1[15]);
		Matrix[168] += detJ*(G * B1[0] * B1[15] + mu * B1[2] * B1[17] + G * B1[1] * B1[16]);
		Matrix[169] += detJ*(lambda * B1[1] * B1[17] + G * B1[2] * B1[16]);
		Matrix[170] += detJ*(lambda * B1[0] * B1[17] + G * B1[2] * B1[15]);
		Matrix[171] += detJ*(mu * (B1[18] * B1[18]) + G * (B1[19] * B1[19]) + G * (B1[20] * B1[20]));
		Matrix[172] += detJ*(lambda * B1[17] * B1[18] + G * B1[15] * B1[20]);
		Matrix[173] += detJ*(lambda * B1[16] * B1[18] + G * B1[15] * B1[19]);
		Matrix[174] += detJ*(mu * B1[15] * B1[18] + G * B1[16] * B1[19] + G * B1[17] * B1[20]);
		Matrix[175] += detJ*(lambda * B1[14] * B1[18] + G * B1[12] * B1[20]);
		Matrix[176] += detJ*(lambda * B1[13] * B1[18] + G * B1[12] * B1[19]);
		Matrix[177] += detJ*(mu * B1[12] * B1[18] + G * B1[13] * B1[19] + G * B1[14] * B1[20]);
		Matrix[178] += detJ*(lambda * B1[11] * B1[18] + G * B1[9] * B1[20]);
		Matrix[179] += detJ*(lambda * B1[10] * B1[18] + G * B1[9] * B1[19]);
		Matrix[180] += detJ*(mu * B1[9] * B1[18] + G * B1[10] * B1[19] + G * B1[11] * B1[20]);
		Matrix[181] += detJ*(lambda * B1[8] * B1[18] + G * B1[6] * B1[20]);
		Matrix[182] += detJ*(lambda * B1[7] * B1[18] + G * B1[6] * B1[19]);
		Matrix[183] += detJ*(mu * B1[6] * B1[18] + G * B1[7] * B1[19] + G * B1[8] * B1[20]);
		Matrix[184] += detJ*(lambda * B1[5] * B1[18] + G * B1[3] * B1[20]);
		Matrix[185] += detJ*(lambda * B1[4] * B1[18] + G * B1[3] * B1[19]);
		Matrix[186] += detJ*(mu * B1[3] * B1[18] + G * B1[4] * B1[19] + G * B1[5] * B1[20]);
		Matrix[187] += detJ*(lambda * B1[2] * B1[18] + G * B1[0] * B1[20]);
		Matrix[188] += detJ*(lambda * B1[1] * B1[18] + G * B1[0] * B1[19]);
		Matrix[189] += detJ*(mu * B1[0] * B1[18] + G * B1[1] * B1[19] + G * B1[2] * B1[20]);
		Matrix[190] += detJ*(mu * (B1[19] * B1[19]) + G * (B1[18] * B1[18]) + G * (B1[20] * B1[20]));
		Matrix[191] += detJ*(lambda * B1[18] * B1[19] + G * B1[18] * B1[19]);
		Matrix[192] += detJ*(lambda * B1[17] * B1[19] + G * B1[16] * B1[20]);
		Matrix[193] += detJ*(mu * B1[16] * B1[19] + G * B1[15] * B1[18] + G * B1[17] * B1[20]);
		Matrix[194] += detJ*(lambda * B1[15] * B1[19] + G * B1[16] * B1[18]);
		Matrix[195] += detJ*(lambda * B1[14] * B1[19] + G * B1[13] * B1[20]);
		Matrix[196] += detJ*(mu * B1[13] * B1[19] + G * B1[12] * B1[18] + G * B1[14] * B1[20]);
		Matrix[197] += detJ*(lambda * B1[12] * B1[19] + G * B1[13] * B1[18]);
		Matrix[198] += detJ*(lambda * B1[11] * B1[19] + G * B1[10] * B1[20]);
		Matrix[199] += detJ*(mu * B1[10] * B1[19] + G * B1[9] * B1[18] + G * B1[11] * B1[20]);
		Matrix[200] += detJ*(lambda * B1[9] * B1[19] + G * B1[10] * B1[18]);
		Matrix[201] += detJ*(lambda * B1[8] * B1[19] + G * B1[7] * B1[20]);
		Matrix[202] += detJ*(mu * B1[7] * B1[19] + G * B1[6] * B1[18] + G * B1[8] * B1[20]);
		Matrix[203] += detJ*(lambda * B1[6] * B1[19] + G * B1[7] * B1[18]);
		Matrix[204] += detJ*(lambda * B1[5] * B1[19] + G * B1[4] * B1[20]);
		Matrix[205] += detJ*(mu * B1[4] * B1[19] + G * B1[3] * B1[18] + G * B1[5] * B1[20]);
		Matrix[206] += detJ*(lambda * B1[3] * B1[19] + G * B1[4] * B1[18]);
		Matrix[207] += detJ*(lambda * B1[2] * B1[19] + G * B1[1] * B1[20]);
		Matrix[208] += detJ*(mu * B1[1] * B1[19] + G * B1[0] * B1[18] + G * B1[2] * B1[20]);
		Matrix[209] += detJ*(lambda * B1[0] * B1[19] + G * B1[1] * B1[18]);
		Matrix[210] += detJ*(mu * (B1[20] * B1[20]) + G * (B1[18] * B1[18]) + G * (B1[19] * B1[19]));
		Matrix[211] += detJ*(lambda * B1[19] * B1[20] + G * B1[19] * B1[20]);
		Matrix[212] += detJ*(lambda * B1[18] * B1[20] + G * B1[18] * B1[20]);
		Matrix[213] += detJ*(G * B1[15] * B1[18] + mu * B1[17] * B1[20] + G * B1[16] * B1[19]);
		Matrix[214] += detJ*(lambda * B1[16] * B1[20] + G * B1[17] * B1[19]);
		Matrix[215] += detJ*(lambda * B1[15] * B1[20] + G * B1[17] * B1[18]);
		Matrix[216] += detJ*(G * B1[12] * B1[18] + mu * B1[14] * B1[20] + G * B1[13] * B1[19]);
		Matrix[217] += detJ*(lambda * B1[13] * B1[20] + G * B1[14] * B1[19]);
		Matrix[218] += detJ*(lambda * B1[12] * B1[20] + G * B1[14] * B1[18]);
		Matrix[219] += detJ*(G * B1[9] * B1[18] + mu * B1[11] * B1[20] + G * B1[10] * B1[19]);
		Matrix[220] += detJ*(lambda * B1[10] * B1[20] + G * B1[11] * B1[19]);
		Matrix[221] += detJ*(lambda * B1[9] * B1[20] + G * B1[11] * B1[18]);
		Matrix[222] += detJ*(G * B1[6] * B1[18] + mu * B1[8] * B1[20] + G * B1[7] * B1[19]);
		Matrix[223] += detJ*(lambda * B1[7] * B1[20] + G * B1[8] * B1[19]);
		Matrix[224] += detJ*(lambda * B1[6] * B1[20] + G * B1[8] * B1[18]);
		Matrix[225] += detJ*(G * B1[3] * B1[18] + mu * B1[5] * B1[20] + G * B1[4] * B1[19]);
		Matrix[226] += detJ*(lambda * B1[4] * B1[20] + G * B1[5] * B1[19]);
		Matrix[227] += detJ*(lambda * B1[3] * B1[20] + G * B1[5] * B1[18]);
		Matrix[228] += detJ*(G * B1[0] * B1[18] + mu * B1[2] * B1[20] + G * B1[1] * B1[19]);
		Matrix[229] += detJ*(lambda * B1[1] * B1[20] + G * B1[2] * B1[19]);
		Matrix[230] += detJ*(lambda * B1[0] * B1[20] + G * B1[2] * B1[18]);
		Matrix[231] += detJ*(mu * (B1[21] * B1[21]) + G * (B1[22] * B1[22]) + G * (B1[23] * B1[23]));
		Matrix[232] += detJ*(lambda * B1[20] * B1[21] + G * B1[18] * B1[23]);
		Matrix[233] += detJ*(lambda * B1[19] * B1[21] + G * B1[18] * B1[22]);
		Matrix[234] += detJ*(mu * B1[18] * B1[21] + G * B1[19] * B1[22] + G * B1[20] * B1[23]);
		Matrix[235] += detJ*(lambda * B1[17] * B1[21] + G * B1[15] * B1[23]);
		Matrix[236] += detJ*(lambda * B1[16] * B1[21] + G * B1[15] * B1[22]);
		Matrix[237] += detJ*(mu * B1[15] * B1[21] + G * B1[16] * B1[22] + G * B1[17] * B1[23]);
		Matrix[238] += detJ*(lambda * B1[14] * B1[21] + G * B1[12] * B1[23]);
		Matrix[239] += detJ*(lambda * B1[13] * B1[21] + G * B1[12] * B1[22]);
		Matrix[240] += detJ*(mu * B1[12] * B1[21] + G * B1[13] * B1[22] + G * B1[14] * B1[23]);
		Matrix[241] += detJ*(lambda * B1[11] * B1[21] + G * B1[9] * B1[23]);
		Matrix[242] += detJ*(lambda * B1[10] * B1[21] + G * B1[9] * B1[22]);
		Matrix[243] += detJ*(mu * B1[9] * B1[21] + G * B1[10] * B1[22] + G * B1[11] * B1[23]);
		Matrix[244] += detJ*(lambda * B1[8] * B1[21] + G * B1[6] * B1[23]);
		Matrix[245] += detJ*(lambda * B1[7] * B1[21] + G * B1[6] * B1[22]);
		Matrix[246] += detJ*(mu * B1[6] * B1[21] + G * B1[7] * B1[22] + G * B1[8] * B1[23]);
		Matrix[247] += detJ*(lambda * B1[5] * B1[21] + G * B1[3] * B1[23]);
		Matrix[248] += detJ*(lambda * B1[4] * B1[21] + G * B1[3] * B1[22]);
		Matrix[249] += detJ*(mu * B1[3] * B1[21] + G * B1[4] * B1[22] + G * B1[5] * B1[23]);
		Matrix[250] += detJ*(lambda * B1[2] * B1[21] + G * B1[0] * B1[23]);
		Matrix[251] += detJ*(lambda * B1[1] * B1[21] + G * B1[0] * B1[22]);
		Matrix[252] += detJ*(mu * B1[0] * B1[21] + G * B1[1] * B1[22] + G * B1[2] * B1[23]);
		Matrix[253] += detJ*(mu * (B1[22] * B1[22]) + G * (B1[21] * B1[21]) + G * (B1[23] * B1[23]));
		Matrix[254] += detJ*(lambda * B1[21] * B1[22] + G * B1[21] * B1[22]);
		Matrix[255] += detJ*(lambda * B1[20] * B1[22] + G * B1[19] * B1[23]);
		Matrix[256] += detJ*(mu * B1[19] * B1[22] + G * B1[18] * B1[21] + G * B1[20] * B1[23]);
		Matrix[257] += detJ*(lambda * B1[18] * B1[22] + G * B1[19] * B1[21]);
		Matrix[258] += detJ*(lambda * B1[17] * B1[22] + G * B1[16] * B1[23]);
		Matrix[259] += detJ*(mu * B1[16] * B1[22] + G * B1[15] * B1[21] + G * B1[17] * B1[23]);
		Matrix[260] += detJ*(lambda * B1[15] * B1[22] + G * B1[16] * B1[21]);
		Matrix[261] += detJ*(lambda * B1[14] * B1[22] + G * B1[13] * B1[23]);
		Matrix[262] += detJ*(mu * B1[13] * B1[22] + G * B1[12] * B1[21] + G * B1[14] * B1[23]);
		Matrix[263] += detJ*(lambda * B1[12] * B1[22] + G * B1[13] * B1[21]);
		Matrix[264] += detJ*(lambda * B1[11] * B1[22] + G * B1[10] * B1[23]);
		Matrix[265] += detJ*(mu * B1[10] * B1[22] + G * B1[9] * B1[21] + G * B1[11] * B1[23]);
		Matrix[266] += detJ*(lambda * B1[9] * B1[22] + G * B1[10] * B1[21]);
		Matrix[267] += detJ*(lambda * B1[8] * B1[22] + G * B1[7] * B1[23]);
		Matrix[268] += detJ*(mu * B1[7] * B1[22] + G * B1[6] * B1[21] + G * B1[8] * B1[23]);
		Matrix[269] += detJ*(lambda * B1[6] * B1[22] + G * B1[7] * B1[21]);
		Matrix[270] += detJ*(lambda * B1[5] * B1[22] + G * B1[4] * B1[23]);
		Matrix[271] += detJ*(mu * B1[4] * B1[22] + G * B1[3] * B1[21] + G * B1[5] * B1[23]);
		Matrix[272] += detJ*(lambda * B1[3] * B1[22] + G * B1[4] * B1[21]);
		Matrix[273] += detJ*(lambda * B1[2] * B1[22] + G * B1[1] * B1[23]);
		Matrix[274] += detJ*(mu * B1[1] * B1[22] + G * B1[0] * B1[21] + G * B1[2] * B1[23]);
		Matrix[275] += detJ*(lambda * B1[0] * B1[22] + G * B1[1] * B1[21]);
		Matrix[276] += detJ*(mu * (B1[23] * B1[23]) + G * (B1[21] * B1[21]) + G * (B1[22] * B1[22]));
		Matrix[277] += detJ*(lambda * B1[22] * B1[23] + G * B1[22] * B1[23]);
		Matrix[278] += detJ*(lambda * B1[21] * B1[23] + G * B1[21] * B1[23]);
		Matrix[279] += detJ*(G * B1[18] * B1[21] + mu * B1[20] * B1[23] + G * B1[19] * B1[22]);
		Matrix[280] += detJ*(lambda * B1[19] * B1[23] + G * B1[20] * B1[22]);
		Matrix[281] += detJ*(lambda * B1[18] * B1[23] + G * B1[20] * B1[21]);
		Matrix[282] += detJ*(G * B1[15] * B1[21] + mu * B1[17] * B1[23] + G * B1[16] * B1[22]);
		Matrix[283] += detJ*(lambda * B1[16] * B1[23] + G * B1[17] * B1[22]);
		Matrix[284] += detJ*(lambda * B1[15] * B1[23] + G * B1[17] * B1[21]);
		Matrix[285] += detJ*(G * B1[12] * B1[21] + mu * B1[14] * B1[23] + G * B1[13] * B1[22]);
		Matrix[286] += detJ*(lambda * B1[13] * B1[23] + G * B1[14] * B1[22]);
		Matrix[287] += detJ*(lambda * B1[12] * B1[23] + G * B1[14] * B1[21]);
		Matrix[288] += detJ*(G * B1[9] * B1[21] + mu * B1[11] * B1[23] + G * B1[10] * B1[22]);
		Matrix[289] += detJ*(lambda * B1[10] * B1[23] + G * B1[11] * B1[22]);
		Matrix[290] += detJ*(lambda * B1[9] * B1[23] + G * B1[11] * B1[21]);
		Matrix[291] += detJ*(G * B1[6] * B1[21] + mu * B1[8] * B1[23] + G * B1[7] * B1[22]);
		Matrix[292] += detJ*(lambda * B1[7] * B1[23] + G * B1[8] * B1[22]);
		Matrix[293] += detJ*(lambda * B1[6] * B1[23] + G * B1[8] * B1[21]);
		Matrix[294] += detJ*(G * B1[3] * B1[21] + mu * B1[5] * B1[23] + G * B1[4] * B1[22]);
		Matrix[295] += detJ*(lambda * B1[4] * B1[23] + G * B1[5] * B1[22]);
		Matrix[296] += detJ*(lambda * B1[3] * B1[23] + G * B1[5] * B1[21]);
		Matrix[297] += detJ*(G * B1[0] * B1[21] + mu * B1[2] * B1[23] + G * B1[1] * B1[22]);
		Matrix[298] += detJ*(lambda * B1[1] * B1[23] + G * B1[2] * B1[22]);
		Matrix[299] += detJ*(lambda * B1[0] * B1[23] + G * B1[2] * B1[21]);
	}
}


//	Calculate element stress 
void C8H::ElementStress(double* stress8H, double* Displacement)
{
	// Get nodal displacements [LM can be used here]
	double Disp[24];
	for (unsigned int i = 0; i < 24; i++)
	{
		if (LocationMatrix_[i])
			Disp[i] = Displacement[LocationMatrix_[i]-1];
		else
			Disp[i] = 0.0;
	}

	// Construct constitutive matrix
	C8HMaterial* material = static_cast<C8HMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double v = material->Nu;
	double k = material->E * (1.0 -v)/(1.0 +v)/(1.0 -2.0 *v);
	double D[3];
	D[0] = k;
	D[1] = k * v / (1.0 - v);
	D[2] = k * (1.0 - 2.0 * v) / 2.0 / (1.0 - v);

	// Construct coordinate matrix
	double COORXYZ[24];
	for (unsigned int i=0;i<8;i++)
	{
		for (unsigned int j=0;j<3;j++)
		{
			COORXYZ[3 * i + j]=nodes_[i]->XYZ[j];
		}
	}

	// Construct Jacobi matrix
	const double xi8[8] = { 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 };
	const double eta8[8] = { -0.577350269189626 , 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,0.577350269189626 , 0.577350269189626 ,-0.577350269189626 };
	const double zeta8[8] = { -0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 , 0.577350269189626 , 0.577350269189626 };

	double stressXYZ[6][8];	// 8 gauss points, 6 stress components
	for (unsigned p = 0; p < 8; p++)
	{
		double xi   = xi8[p];
		double eta  = eta8[p];
		double zeta = zeta8[p];

		double GN[12];
		GN[0] = (1.0-eta)*(1.0-zeta) / 8.0;
		GN[1] = (1.0+eta)*(1.0-zeta) / 8.0;
		GN[2] = (1.0-eta)*(1.0+zeta) / 8.0;
		GN[3] = (1.0+eta)*(1.0+zeta) / 8.0;
		GN[4] = (1.0+xi)*(1.0-zeta) / 8.0;
		GN[5] = (1.0-xi)*(1.0-zeta) / 8.0;
		GN[6] = (1.0+xi)*(1.0+zeta) / 8.0;
		GN[7] = (1.0-xi)*(1.0+zeta) / 8.0;
		GN[8] = (1.0+xi)*(1.0-eta) / 8.0;
		GN[9] = (1.0+xi)*(1.0+eta) / 8.0;
		GN[10] = (1.0-xi)*(1.0+eta) / 8.0;
		GN[11] = (1.0-xi)*(1.0-eta) / 8.0;

		double J[9];
		J[0] = COORXYZ[0] * GN[0] + COORXYZ[3] * GN[1] - COORXYZ[6] * GN[1] - COORXYZ[9] * GN[0] + COORXYZ[12] * GN[2] + COORXYZ[15] * GN[3] - COORXYZ[18] * GN[3] - COORXYZ[21] * GN[2];
		J[1] = -COORXYZ[0] * GN[4] + COORXYZ[3] * GN[4] + COORXYZ[6] * GN[5] - COORXYZ[9] * GN[5] - COORXYZ[12] * GN[6] + COORXYZ[15] * GN[6] + COORXYZ[18] * GN[7] - COORXYZ[21] * GN[7];
		J[2] = -COORXYZ[0] * GN[8] - COORXYZ[3] * GN[9] - COORXYZ[6] * GN[10] - COORXYZ[9] * GN[11] + COORXYZ[12] * GN[8] + COORXYZ[15] * GN[9] + COORXYZ[18] * GN[10] + COORXYZ[21] * GN[11];
		J[3] = COORXYZ[1] * GN[0] + COORXYZ[4] * GN[1] - COORXYZ[7] * GN[1] - COORXYZ[10] * GN[0] + COORXYZ[13] * GN[2] + COORXYZ[16] * GN[3] - COORXYZ[19] * GN[3] - COORXYZ[22] * GN[2];
		J[4] = -COORXYZ[1] * GN[4] + COORXYZ[4] * GN[4] + COORXYZ[7] * GN[5] - COORXYZ[10] * GN[5] - COORXYZ[13] * GN[6] + COORXYZ[16] * GN[6] + COORXYZ[19] * GN[7] - COORXYZ[22] * GN[7];
		J[5] = -COORXYZ[1] * GN[8] - COORXYZ[4] * GN[9] - COORXYZ[7] * GN[10] - COORXYZ[10] * GN[11] + COORXYZ[13] * GN[8] + COORXYZ[16] * GN[9] + COORXYZ[19] * GN[10] + COORXYZ[22] * GN[11];
		J[6] = COORXYZ[2] * GN[0] + COORXYZ[5] * GN[1] - COORXYZ[8] * GN[1] - COORXYZ[11] * GN[0] + COORXYZ[14] * GN[2] + COORXYZ[17] * GN[3] - COORXYZ[20] * GN[3] - COORXYZ[23] * GN[2];
		J[7] = -COORXYZ[2] * GN[4] + COORXYZ[5] * GN[4] + COORXYZ[8] * GN[5] - COORXYZ[11] * GN[5] - COORXYZ[14] * GN[6] + COORXYZ[17] * GN[6] + COORXYZ[20] * GN[7] - COORXYZ[23] * GN[7];
		J[8] = -COORXYZ[2] * GN[8] - COORXYZ[5] * GN[9] - COORXYZ[8] * GN[10] - COORXYZ[11] * GN[11] + COORXYZ[14] * GN[8] + COORXYZ[17] * GN[9] + COORXYZ[20] * GN[10] + COORXYZ[23] * GN[11];

		double detJ = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];

		double InvJ[9];
		InvJ[0] = (J[4]*J[8]-J[5]*J[7])/detJ;
		InvJ[1] = -(J[1]*J[8]-J[2]*J[7])/detJ;
		InvJ[2] = (J[1]*J[5]-J[2]*J[4])/detJ;
		InvJ[3] = -(J[3]*J[8]-J[5]*J[6])/detJ;
		InvJ[4] = (J[0]*J[8]-J[2]*J[6])/detJ;
		InvJ[5] = -(J[0]*J[5]-J[2]*J[3])/detJ;
		InvJ[6] = (J[3]*J[7]-J[4]*J[6])/detJ;
		InvJ[7] = -(J[0]*J[7]-J[1]*J[6])/detJ;
		InvJ[8] = (J[0]*J[4]-J[1]*J[3])/detJ;

		double kerB[24];
		kerB[0] = GN[0] * InvJ[0] - GN[4] * InvJ[3] - GN[8] * InvJ[6];
		kerB[1] = GN[0] * InvJ[1] - GN[4] * InvJ[4] - GN[8] * InvJ[7];
		kerB[2] = GN[0] * InvJ[2] - GN[4] * InvJ[5] - GN[8] * InvJ[8];
		kerB[3] = GN[1] * InvJ[0] + GN[4] * InvJ[3] - GN[9] * InvJ[6];
		kerB[4] = GN[1] * InvJ[1] + GN[4] * InvJ[4] - GN[9] * InvJ[7];
		kerB[5] = GN[1] * InvJ[2] + GN[4] * InvJ[5] - GN[9] * InvJ[8];
		kerB[6] = -GN[1] * InvJ[0] + GN[5] * InvJ[3] - GN[10] * InvJ[6];
		kerB[7] = -GN[1] * InvJ[1] + GN[5] * InvJ[4] - GN[10] * InvJ[7];
		kerB[8] = -GN[1] * InvJ[2] + GN[5] * InvJ[5] - GN[10] * InvJ[8];
		kerB[9] = -GN[0] * InvJ[0] - GN[5] * InvJ[3] - GN[11] * InvJ[6];
		kerB[10] = -GN[0] * InvJ[1] - GN[5] * InvJ[4] - GN[11] * InvJ[7];
		kerB[11] = -GN[0] * InvJ[2] - GN[5] * InvJ[5] - GN[11] * InvJ[8];
		kerB[12] = GN[2] * InvJ[0] - GN[6] * InvJ[3] + GN[8] * InvJ[6];
		kerB[13] = GN[2] * InvJ[1] - GN[6] * InvJ[4] + GN[8] * InvJ[7];
		kerB[14] = GN[2] * InvJ[2] - GN[6] * InvJ[5] + GN[8] * InvJ[8];
		kerB[15] = GN[3] * InvJ[0] + GN[6] * InvJ[3] + GN[9] * InvJ[6];
		kerB[16] = GN[3] * InvJ[1] + GN[6] * InvJ[4] + GN[9] * InvJ[7];
		kerB[17] = GN[3] * InvJ[2] + GN[6] * InvJ[5] + GN[9] * InvJ[8];
		kerB[18] = -GN[3] * InvJ[0] + GN[7] * InvJ[3] + GN[10] * InvJ[6];
		kerB[19] = -GN[3] * InvJ[1] + GN[7] * InvJ[4] + GN[10] * InvJ[7];
		kerB[20] = -GN[3] * InvJ[2] + GN[7] * InvJ[5] + GN[10] * InvJ[8];
		kerB[21] = -GN[2] * InvJ[0] - GN[7] * InvJ[3] + GN[11] * InvJ[6];
		kerB[22] = -GN[2] * InvJ[1] - GN[7] * InvJ[4] + GN[11] * InvJ[7];
		kerB[23] = -GN[2] * InvJ[2] - GN[7] * InvJ[5] + GN[11] * InvJ[8];
		
		double EPS[6];
		clear(EPS,6);
		EPS[0]=Disp[0] * kerB[0]+ Disp[3] * kerB[3]+Disp[6] * kerB[6]+Disp[9] * kerB[9]+Disp[12] * kerB[12]+ Disp[15] * kerB[15]+Disp[18] * kerB[18]+Disp[21] * kerB[21];
		EPS[1]=Disp[1] * kerB[1]+ Disp[4] * kerB[4]+Disp[7] * kerB[7]+Disp[10] * kerB[10]+Disp[13] * kerB[13]+ Disp[16] * kerB[16]+Disp[19] * kerB[19]+Disp[22] * kerB[22];
		EPS[2]=Disp[2] * kerB[2]+ Disp[5] * kerB[5]+Disp[8] * kerB[8]+Disp[11] * kerB[11]+Disp[14] * kerB[14]+ Disp[17] * kerB[17]+Disp[20] * kerB[20]+Disp[23] * kerB[23];
		stressXYZ[0][p]=EPS[0]*D[0]+(EPS[1]+EPS[2])*D[1];
		stressXYZ[1][p]=EPS[1]*D[0]+(EPS[0]+EPS[2])*D[1];
		stressXYZ[2][p]=EPS[2]*D[0]+(EPS[0]+EPS[1])*D[1];
//		stressXYZ[0][p] = D[0] * Disp[0] * kerB[0] + D[1] * Disp[1] * kerB[1] + D[1] * Disp[2] * kerB[2] + D[0] * Disp[3] * kerB[3] + D[1] * Disp[4] * kerB[4] + D[1] * Disp[5] * kerB[5] + D[0] * Disp[6] * kerB[6] + D[1] * Disp[7] * kerB[7] + D[1] * Disp[8] * kerB[8] + D[0] * Disp[9] * kerB[9] + D[1] * Disp[10] * kerB[10] + D[1] * Disp[11] * kerB[11] + D[0] * Disp[12] * kerB[12] + D[1] * Disp[13] * kerB[13] + D[1] * Disp[14] * kerB[14] + D[0] * Disp[15] * kerB[15] + D[1] * Disp[16] * kerB[16] + D[1] * Disp[17] * kerB[17] + D[0] * Disp[18] * kerB[18] + D[1] * Disp[19] * kerB[19] + D[1] * Disp[20] * kerB[20] + D[0] * Disp[21] * kerB[21] + D[1] * Disp[22] * kerB[22] + D[1] * Disp[23] * kerB[23];
//		stressXYZ[1][p] = D[1] * Disp[0] * kerB[0] + D[0] * Disp[1] * kerB[1] + D[1] * Disp[2] * kerB[2] + D[1] * Disp[3] * kerB[3] + D[0] * Disp[4] * kerB[4] + D[1] * Disp[5] * kerB[5] + D[1] * Disp[6] * kerB[6] + D[0] * Disp[7] * kerB[7] + D[1] * Disp[8] * kerB[8] + D[1] * Disp[9] * kerB[9] + D[0] * Disp[10] * kerB[10] + D[1] * Disp[11] * kerB[11] + D[1] * Disp[12] * kerB[12] + D[0] * Disp[13] * kerB[13] + D[1] * Disp[14] * kerB[14] + D[1] * Disp[15] * kerB[15] + D[0] * Disp[16] * kerB[16] + D[1] * Disp[17] * kerB[17] + D[1] * Disp[18] * kerB[18] + D[0] * Disp[19] * kerB[19] + D[1] * Disp[20] * kerB[20] + D[1] * Disp[21] * kerB[21] + D[0] * Disp[22] * kerB[22] + D[1] * Disp[23] * kerB[23];
//		stressXYZ[2][p] = D[1] * Disp[0] * kerB[0] + D[1] * Disp[1] * kerB[1] + D[0] * Disp[2] * kerB[2] + D[1] * Disp[3] * kerB[3] + D[1] * Disp[4] * kerB[4] + D[0] * Disp[5] * kerB[5] + D[1] * Disp[6] * kerB[6] + D[1] * Disp[7] * kerB[7] + D[0] * Disp[8] * kerB[8] + D[1] * Disp[9] * kerB[9] + D[1] * Disp[10] * kerB[10] + D[0] * Disp[11] * kerB[11] + D[1] * Disp[12] * kerB[12] + D[1] * Disp[13] * kerB[13] + D[0] * Disp[14] * kerB[14] + D[1] * Disp[15] * kerB[15] + D[1] * Disp[16] * kerB[16] + D[0] * Disp[17] * kerB[17] + D[1] * Disp[18] * kerB[18] + D[1] * Disp[19] * kerB[19] + D[0] * Disp[20] * kerB[20] + D[1] * Disp[21] * kerB[21] + D[1] * Disp[22] * kerB[22] + D[0] * Disp[23] * kerB[23];
		stressXYZ[3][p] = D[2] * Disp[0] * kerB[1] + D[2] * Disp[1] * kerB[0] + D[2] * Disp[3] * kerB[4] + D[2] * Disp[4] * kerB[3] + D[2] * Disp[6] * kerB[7] + D[2] * Disp[7] * kerB[6] + D[2] * Disp[9] * kerB[10] + D[2] * Disp[10] * kerB[9] + D[2] * Disp[12] * kerB[13] + D[2] * Disp[13] * kerB[12] + D[2] * Disp[15] * kerB[16] + D[2] * Disp[16] * kerB[15] + D[2] * Disp[18] * kerB[19] + D[2] * Disp[19] * kerB[18] + D[2] * Disp[21] * kerB[22] + D[2] * Disp[22] * kerB[21];
		stressXYZ[4][p] = D[2] * Disp[1] * kerB[2] + D[2] * Disp[2] * kerB[1] + D[2] * Disp[4] * kerB[5] + D[2] * Disp[5] * kerB[4] + D[2] * Disp[7] * kerB[8] + D[2] * Disp[8] * kerB[7] + D[2] * Disp[10] * kerB[11] + D[2] * Disp[11] * kerB[10] + D[2] * Disp[13] * kerB[14] + D[2] * Disp[14] * kerB[13] + D[2] * Disp[16] * kerB[17] + D[2] * Disp[17] * kerB[16] + D[2] * Disp[19] * kerB[20] + D[2] * Disp[20] * kerB[19] + D[2] * Disp[22] * kerB[23] + D[2] * Disp[23] * kerB[22];
		stressXYZ[5][p] = D[2] * Disp[0] * kerB[2] + D[2] * Disp[2] * kerB[0] + D[2] * Disp[3] * kerB[5] + D[2] * Disp[5] * kerB[3] + D[2] * Disp[6] * kerB[8] + D[2] * Disp[8] * kerB[6] + D[2] * Disp[9] * kerB[11] + D[2] * Disp[11] * kerB[9] + D[2] * Disp[12] * kerB[14] + D[2] * Disp[14] * kerB[12] + D[2] * Disp[15] * kerB[17] + D[2] * Disp[17] * kerB[15] + D[2] * Disp[18] * kerB[20] + D[2] * Disp[20] * kerB[18] + D[2] * Disp[21] * kerB[23] + D[2] * Disp[23] * kerB[21];
	}

	// stress recovery for stress on nodes
	double interpo[4] = {2.549038105676658, -0.683012701892219, 0.183012701892219, -0.049038105676658};
	double recovery[8];
	for (unsigned i = 0; i < 6; i++)
	{
		recovery[0] = interpo[0]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][3] + interpo[2]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][7] + interpo[3]*stressXYZ[i][6];
		recovery[1] = interpo[0]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][6] + interpo[3]*stressXYZ[i][7];
		recovery[2] = interpo[0]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][5] + interpo[3]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][7];
		recovery[3] = interpo[1]*stressXYZ[i][0] + interpo[0]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][4] + interpo[1]*stressXYZ[i][7] + interpo[2]*stressXYZ[i][6] + interpo[3]*stressXYZ[i][5];
		recovery[4] = interpo[1]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][1] + interpo[0]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][3] + interpo[3]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][5] + interpo[1]*stressXYZ[i][7] + interpo[2]*stressXYZ[i][6];
		recovery[5] = interpo[1]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][2] + interpo[0]*stressXYZ[i][5] + interpo[1]*stressXYZ[i][4] + interpo[3]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][7];
		recovery[6] = interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][1] + interpo[3]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][3] + interpo[0]*stressXYZ[i][6] + interpo[1]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][4] + interpo[1]*stressXYZ[i][7];
		recovery[7] = interpo[2]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][3] + interpo[2]*stressXYZ[i][2] + interpo[3]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][4] + interpo[0]*stressXYZ[i][7] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][5];
		for (unsigned j = 0; j < 8; j++)
		{
			stress8H[6 * j + i] = recovery[j];
		}
	}
}





void  C8H::ElementPostInfo(double* stress, double* Displacement , double* PrePositions, double* PostPositions)
{
	// get original position: preposition
	for (unsigned int i =0 ; i<3; i++)
	{
		for (unsigned int j=0;j < 8; j++)
		{
			PrePositions[i+3*j] = nodes_[j]->XYZ[i];	 			
		}
	}
    double Disp[24];
	// Get nodal displacements [LM can be used here]
	for (unsigned int i = 0; i < 24; i++)
	{
	
		if (LocationMatrix_[i])
			//locatiion matrix start from 1 not 0

		{Disp[i] = Displacement[LocationMatrix_[i]-1];}
		else
		{Disp[i] = 0.0;}

		PostPositions[i] = PrePositions[i] + Disp[i];

	}

	// Construct constitutive matrix
	C8HMaterial* material = static_cast<C8HMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double v = material->Nu;
	double k = material->E * (1.0 -v)/(1.0 +v)/(1.0 -2.0 *v);
	double D[3];
	D[0] = k;
	D[1] = k * v / (1.0 - v);
	D[2] = k * (1.0 - 2.0 * v) / 2.0 / (1.0 - v);


	// Construct Jacobi matrix
	const double xi8[8] = { 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 };
	const double eta8[8] = { -0.577350269189626 , 0.577350269189626 , 0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,0.577350269189626 , 0.577350269189626 ,-0.577350269189626 };
	const double zeta8[8] = { -0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 ,-0.577350269189626 , 0.577350269189626 ,0.577350269189626 , 0.577350269189626 , 0.577350269189626 };

	double stressXYZ[6][8];	// 8 gauss points, 6 stress components
	for (unsigned p = 0; p < 8; p++)
	{
		double xi   = xi8[p];
		double eta  = eta8[p];
		double zeta = zeta8[p];

		double GN[12];
		GN[0] = (1.0-eta)*(1.0-zeta) / 8.0;
		GN[1] = (1.0+eta)*(1.0-zeta) / 8.0;
		GN[2] = (1.0-eta)*(1.0+zeta) / 8.0;
		GN[3] = (1.0+eta)*(1.0+zeta) / 8.0;
		GN[4] = (1.0+xi)*(1.0-zeta) / 8.0;
		GN[5] = (1.0-xi)*(1.0-zeta) / 8.0;
		GN[6] = (1.0+xi)*(1.0+zeta) / 8.0;
		GN[7] = (1.0-xi)*(1.0+zeta) / 8.0;
		GN[8] = (1.0+xi)*(1.0-eta) / 8.0;
		GN[9] = (1.0+xi)*(1.0+eta) / 8.0;
		GN[10] = (1.0-xi)*(1.0+eta) / 8.0;
		GN[11] = (1.0-xi)*(1.0-eta) / 8.0;

		double J[9];
		J[0] = PrePositions[0] * GN[0] + PrePositions[3] * GN[1] - PrePositions[6] * GN[1] - PrePositions[9] * GN[0] + PrePositions[12] * GN[2] + PrePositions[15] * GN[3] - PrePositions[18] * GN[3] - PrePositions[21] * GN[2];
		J[1] = -PrePositions[0] * GN[4] + PrePositions[3] * GN[4] + PrePositions[6] * GN[5] - PrePositions[9] * GN[5] - PrePositions[12] * GN[6] + PrePositions[15] * GN[6] + PrePositions[18] * GN[7] - PrePositions[21] * GN[7];
		J[2] = -PrePositions[0] * GN[8] - PrePositions[3] * GN[9] - PrePositions[6] * GN[10] - PrePositions[9] * GN[11] + PrePositions[12] * GN[8] + PrePositions[15] * GN[9] + PrePositions[18] * GN[10] + PrePositions[21] * GN[11];
		J[3] = PrePositions[1] * GN[0] + PrePositions[4] * GN[1] - PrePositions[7] * GN[1] - PrePositions[10] * GN[0] + PrePositions[13] * GN[2] + PrePositions[16] * GN[3] - PrePositions[19] * GN[3] - PrePositions[22] * GN[2];
		J[4] = -PrePositions[1] * GN[4] + PrePositions[4] * GN[4] + PrePositions[7] * GN[5] - PrePositions[10] * GN[5] - PrePositions[13] * GN[6] + PrePositions[16] * GN[6] + PrePositions[19] * GN[7] - PrePositions[22] * GN[7];
		J[5] = -PrePositions[1] * GN[8] - PrePositions[4] * GN[9] - PrePositions[7] * GN[10] - PrePositions[10] * GN[11] + PrePositions[13] * GN[8] + PrePositions[16] * GN[9] + PrePositions[19] * GN[10] + PrePositions[22] * GN[11];
		J[6] = PrePositions[2] * GN[0] + PrePositions[5] * GN[1] - PrePositions[8] * GN[1] - PrePositions[11] * GN[0] + PrePositions[14] * GN[2] + PrePositions[17] * GN[3] - PrePositions[20] * GN[3] - PrePositions[23] * GN[2];
		J[7] = -PrePositions[2] * GN[4] + PrePositions[5] * GN[4] + PrePositions[8] * GN[5] - PrePositions[11] * GN[5] - PrePositions[14] * GN[6] + PrePositions[17] * GN[6] + PrePositions[20] * GN[7] - PrePositions[23] * GN[7];
		J[8] = -PrePositions[2] * GN[8] - PrePositions[5] * GN[9] - PrePositions[8] * GN[10] - PrePositions[11] * GN[11] + PrePositions[14] * GN[8] + PrePositions[17] * GN[9] + PrePositions[20] * GN[10] + PrePositions[23] * GN[11];

		double detJ = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];

		double InvJ[9];
		InvJ[0] = (J[4]*J[8]-J[5]*J[7])/detJ;
		InvJ[1] = -(J[1]*J[8]-J[2]*J[7])/detJ;
		InvJ[2] = (J[1]*J[5]-J[2]*J[4])/detJ;
		InvJ[3] = -(J[3]*J[8]-J[5]*J[6])/detJ;
		InvJ[4] = (J[0]*J[8]-J[2]*J[6])/detJ;
		InvJ[5] = -(J[0]*J[5]-J[2]*J[3])/detJ;
		InvJ[6] = (J[3]*J[7]-J[4]*J[6])/detJ;
		InvJ[7] = -(J[0]*J[7]-J[1]*J[6])/detJ;
		InvJ[8] = (J[0]*J[4]-J[1]*J[3])/detJ;

		double kerB[24];
		kerB[0] = GN[0] * InvJ[0] - GN[4] * InvJ[3] - GN[8] * InvJ[6];
		kerB[1] = GN[0] * InvJ[1] - GN[4] * InvJ[4] - GN[8] * InvJ[7];
		kerB[2] = GN[0] * InvJ[2] - GN[4] * InvJ[5] - GN[8] * InvJ[8];
		kerB[3] = GN[1] * InvJ[0] + GN[4] * InvJ[3] - GN[9] * InvJ[6];
		kerB[4] = GN[1] * InvJ[1] + GN[4] * InvJ[4] - GN[9] * InvJ[7];
		kerB[5] = GN[1] * InvJ[2] + GN[4] * InvJ[5] - GN[9] * InvJ[8];
		kerB[6] = -GN[1] * InvJ[0] + GN[5] * InvJ[3] - GN[10] * InvJ[6];
		kerB[7] = -GN[1] * InvJ[1] + GN[5] * InvJ[4] - GN[10] * InvJ[7];
		kerB[8] = -GN[1] * InvJ[2] + GN[5] * InvJ[5] - GN[10] * InvJ[8];
		kerB[9] = -GN[0] * InvJ[0] - GN[5] * InvJ[3] - GN[11] * InvJ[6];
		kerB[10] = -GN[0] * InvJ[1] - GN[5] * InvJ[4] - GN[11] * InvJ[7];
		kerB[11] = -GN[0] * InvJ[2] - GN[5] * InvJ[5] - GN[11] * InvJ[8];
		kerB[12] = GN[2] * InvJ[0] - GN[6] * InvJ[3] + GN[8] * InvJ[6];
		kerB[13] = GN[2] * InvJ[1] - GN[6] * InvJ[4] + GN[8] * InvJ[7];
		kerB[14] = GN[2] * InvJ[2] - GN[6] * InvJ[5] + GN[8] * InvJ[8];
		kerB[15] = GN[3] * InvJ[0] + GN[6] * InvJ[3] + GN[9] * InvJ[6];
		kerB[16] = GN[3] * InvJ[1] + GN[6] * InvJ[4] + GN[9] * InvJ[7];
		kerB[17] = GN[3] * InvJ[2] + GN[6] * InvJ[5] + GN[9] * InvJ[8];
		kerB[18] = -GN[3] * InvJ[0] + GN[7] * InvJ[3] + GN[10] * InvJ[6];
		kerB[19] = -GN[3] * InvJ[1] + GN[7] * InvJ[4] + GN[10] * InvJ[7];
		kerB[20] = -GN[3] * InvJ[2] + GN[7] * InvJ[5] + GN[10] * InvJ[8];
		kerB[21] = -GN[2] * InvJ[0] - GN[7] * InvJ[3] + GN[11] * InvJ[6];
		kerB[22] = -GN[2] * InvJ[1] - GN[7] * InvJ[4] + GN[11] * InvJ[7];
		kerB[23] = -GN[2] * InvJ[2] - GN[7] * InvJ[5] + GN[11] * InvJ[8];

		stressXYZ[0][p] = D[0] * Disp[0] * kerB[0] + D[1] * Disp[1] * kerB[1] + D[1] * Disp[2] * kerB[2] + D[0] * Disp[3] * kerB[3] + D[1] * Disp[4] * kerB[4] + D[1] * Disp[5] * kerB[5] + D[0] * Disp[6] * kerB[6] + D[1] * Disp[7] * kerB[7] + D[1] * Disp[8] * kerB[8] + D[0] * Disp[9] * kerB[9] + D[1] * Disp[10] * kerB[10] + D[1] * Disp[11] * kerB[11] + D[0] * Disp[12] * kerB[12] + D[1] * Disp[13] * kerB[13] + D[1] * Disp[14] * kerB[14] + D[0] * Disp[15] * kerB[15] + D[1] * Disp[16] * kerB[16] + D[1] * Disp[17] * kerB[17] + D[0] * Disp[18] * kerB[18] + D[1] * Disp[19] * kerB[19] + D[1] * Disp[20] * kerB[20] + D[0] * Disp[21] * kerB[21] + D[1] * Disp[22] * kerB[22] + D[1] * Disp[23] * kerB[23];
		stressXYZ[1][p] = D[1] * Disp[0] * kerB[0] + D[0] * Disp[1] * kerB[1] + D[1] * Disp[2] * kerB[2] + D[1] * Disp[3] * kerB[3] + D[0] * Disp[4] * kerB[4] + D[1] * Disp[5] * kerB[5] + D[1] * Disp[6] * kerB[6] + D[0] * Disp[7] * kerB[7] + D[1] * Disp[8] * kerB[8] + D[1] * Disp[9] * kerB[9] + D[0] * Disp[10] * kerB[10] + D[1] * Disp[11] * kerB[11] + D[1] * Disp[12] * kerB[12] + D[0] * Disp[13] * kerB[13] + D[1] * Disp[14] * kerB[14] + D[1] * Disp[15] * kerB[15] + D[0] * Disp[16] * kerB[16] + D[1] * Disp[17] * kerB[17] + D[1] * Disp[18] * kerB[18] + D[0] * Disp[19] * kerB[19] + D[1] * Disp[20] * kerB[20] + D[1] * Disp[21] * kerB[21] + D[0] * Disp[22] * kerB[22] + D[1] * Disp[23] * kerB[23];
		stressXYZ[2][p] = D[1] * Disp[0] * kerB[0] + D[1] * Disp[1] * kerB[1] + D[0] * Disp[2] * kerB[2] + D[1] * Disp[3] * kerB[3] + D[1] * Disp[4] * kerB[4] + D[0] * Disp[5] * kerB[5] + D[1] * Disp[6] * kerB[6] + D[1] * Disp[7] * kerB[7] + D[0] * Disp[8] * kerB[8] + D[1] * Disp[9] * kerB[9] + D[1] * Disp[10] * kerB[10] + D[0] * Disp[11] * kerB[11] + D[1] * Disp[12] * kerB[12] + D[1] * Disp[13] * kerB[13] + D[0] * Disp[14] * kerB[14] + D[1] * Disp[15] * kerB[15] + D[1] * Disp[16] * kerB[16] + D[0] * Disp[17] * kerB[17] + D[1] * Disp[18] * kerB[18] + D[1] * Disp[19] * kerB[19] + D[0] * Disp[20] * kerB[20] + D[1] * Disp[21] * kerB[21] + D[1] * Disp[22] * kerB[22] + D[0] * Disp[23] * kerB[23];
		stressXYZ[3][p] = D[2] * Disp[0] * kerB[1] + D[2] * Disp[1] * kerB[0] + D[2] * Disp[3] * kerB[4] + D[2] * Disp[4] * kerB[3] + D[2] * Disp[6] * kerB[7] + D[2] * Disp[7] * kerB[6] + D[2] * Disp[9] * kerB[10] + D[2] * Disp[10] * kerB[9] + D[2] * Disp[12] * kerB[13] + D[2] * Disp[13] * kerB[12] + D[2] * Disp[15] * kerB[16] + D[2] * Disp[16] * kerB[15] + D[2] * Disp[18] * kerB[19] + D[2] * Disp[19] * kerB[18] + D[2] * Disp[21] * kerB[22] + D[2] * Disp[22] * kerB[21];
		stressXYZ[4][p] = D[2] * Disp[1] * kerB[2] + D[2] * Disp[2] * kerB[1] + D[2] * Disp[4] * kerB[5] + D[2] * Disp[5] * kerB[4] + D[2] * Disp[7] * kerB[8] + D[2] * Disp[8] * kerB[7] + D[2] * Disp[10] * kerB[11] + D[2] * Disp[11] * kerB[10] + D[2] * Disp[13] * kerB[14] + D[2] * Disp[14] * kerB[13] + D[2] * Disp[16] * kerB[17] + D[2] * Disp[17] * kerB[16] + D[2] * Disp[19] * kerB[20] + D[2] * Disp[20] * kerB[19] + D[2] * Disp[22] * kerB[23] + D[2] * Disp[23] * kerB[22];
		stressXYZ[5][p] = D[2] * Disp[0] * kerB[2] + D[2] * Disp[2] * kerB[0] + D[2] * Disp[3] * kerB[5] + D[2] * Disp[5] * kerB[3] + D[2] * Disp[6] * kerB[8] + D[2] * Disp[8] * kerB[6] + D[2] * Disp[9] * kerB[11] + D[2] * Disp[11] * kerB[9] + D[2] * Disp[12] * kerB[14] + D[2] * Disp[14] * kerB[12] + D[2] * Disp[15] * kerB[17] + D[2] * Disp[17] * kerB[15] + D[2] * Disp[18] * kerB[20] + D[2] * Disp[20] * kerB[18] + D[2] * Disp[21] * kerB[23] + D[2] * Disp[23] * kerB[21];
	}

	// stress recovery for stress on nodes
	double interpo[4] = {2.549038105676658, -0.683012701892219, 0.183012701892219, -0.049038105676658};
	double recovery[8];
	for (unsigned i = 0; i < 6; i++)
	{
		recovery[0] = interpo[0]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][3] + interpo[2]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][7] + interpo[3]*stressXYZ[i][6];
		recovery[1] = interpo[0]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][6] + interpo[3]*stressXYZ[i][7];
		recovery[2] = interpo[0]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][5] + interpo[3]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][7];
		recovery[3] = interpo[1]*stressXYZ[i][0] + interpo[0]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][4] + interpo[1]*stressXYZ[i][7] + interpo[2]*stressXYZ[i][6] + interpo[3]*stressXYZ[i][5];
		recovery[4] = interpo[1]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][1] + interpo[0]*stressXYZ[i][4] + interpo[2]*stressXYZ[i][3] + interpo[3]*stressXYZ[i][2] + interpo[1]*stressXYZ[i][5] + interpo[1]*stressXYZ[i][7] + interpo[2]*stressXYZ[i][6];
		recovery[5] = interpo[1]*stressXYZ[i][1] + interpo[2]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][2] + interpo[0]*stressXYZ[i][5] + interpo[1]*stressXYZ[i][4] + interpo[3]*stressXYZ[i][3] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][7];
		recovery[6] = interpo[1]*stressXYZ[i][2] + interpo[2]*stressXYZ[i][1] + interpo[3]*stressXYZ[i][0] + interpo[2]*stressXYZ[i][3] + interpo[0]*stressXYZ[i][6] + interpo[1]*stressXYZ[i][5] + interpo[2]*stressXYZ[i][4] + interpo[1]*stressXYZ[i][7];
		recovery[7] = interpo[2]*stressXYZ[i][0] + interpo[1]*stressXYZ[i][3] + interpo[2]*stressXYZ[i][2] + interpo[3]*stressXYZ[i][1] + interpo[1]*stressXYZ[i][4] + interpo[0]*stressXYZ[i][7] + interpo[1]*stressXYZ[i][6] + interpo[2]*stressXYZ[i][5];
		for (unsigned j = 0; j < 8; j++)
		{
			stress[6 * j + i] = recovery[j];
		}
	}

}
//Gravity
double C8H::Gravity()
{
	return 0;
}




