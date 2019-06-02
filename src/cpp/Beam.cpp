#include "Beam.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

CBeam::CBeam()
{
	NEN_ = 2;
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

CBeam::~CBeam()
{}


bool CBeam::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
    ElementMaterial_ = dynamic_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

void CBeam::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

void CBeam::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int	 D = 0; D < 6; D++)
		{
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
		}
		

}

unsigned int CBeam::SizeOfStiffnessMatrix() { return 78; }

void CBeam::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

    
    double DX[3]; //	
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
	double a=material_->a;
	double b=material_->b;
	double x1=material_->x1;
	double y1=material_->y1;
	double x2=material_->x2;
	double y2=material_->y2;
	double Iz = a*a*a*b/12.0 -(a-x1-x2)*(a-x1-x2)*(a-x1-x2)*(b-y1-y2)/12.0;
	double Iy = b*b*b*a/12.0 -(b-y1-y2)*(b-y1-y2)*(b-y1-y2)*(a-x1-x2)/12.0;
	
	

	double k[8];
	 k[0] = material_->E * material_->a * material_->b / L;
     k[1] = 12.0 * material_->E * Iz / (L * L * L);
     k[2] = 12.0 * material_->E * Iy / (L * L * L);
	 k[3] = material_->E * (Iz+Iy) / ((2.0 + 2.0 * material_->Nu) * L);
     k[4] = 4.0 * material_->E * Iy / L;
     k[5] = 4.0 * material_->E * Iz / L;
	 k[6] = 6.0 * material_->E * Iy / (L * L);
     k[7] = 6.0 * material_->E * Iz / (L * L);

		 double n[3][3];
		 n[0][0] = DX[0] / L;
		 n[0][1] = DX[1] / L;
		 n[0][2] = DX[2] / L;
		 n[1][0] = material_->n1;
		 n[1][1] = material_->n2;
		 n[1][2] = material_->n3; 
		 n[2][0] = n[0][1] * n[1][2] - n[0][2] * n[1][1];
		 n[2][1] = n[0][2] * n[1][0] - n[0][0] * n[1][2];
		 n[2][2] = n[0][0] * n[1][1] - n[0][1] * n[1][0]; 

	

	 double N[27];
	 N[0] = n[0][0] * n[0][0];
     N[1] = n[0][1] * n[0][1];
     N[2] = n[0][2] * n[0][2];
     N[3] = n[1][0] * n[1][0];
     N[4] = n[1][1] * n[1][1];
     N[5] = n[1][2] * n[1][2];
     N[6] = n[2][0] * n[2][0];
     N[7] = n[2][1] * n[2][1];
     N[8] = n[2][2] * n[2][2];
     N[9] = n[0][0] * n[0][1];
     N[10] = n[1][0] * n[1][1];
     N[11] = n[2][0] * n[2][1];
     N[12] = n[0][1] * n[0][2];
     N[13] = n[1][1] * n[1][2];
     N[14] = n[2][1] * n[2][2];
     N[15] = n[0][0] * n[0][2];
     N[16] = n[1][0] * n[1][2];
     N[17] = n[2][0] * n[2][2];
     N[18] = n[1][2] * n[2][0];
     N[19] = n[1][0] * n[2][2];
     N[20] = n[1][1] * n[2][0];
     N[21] = n[1][0] * n[2][1];
     N[22] = n[1][2] * n[2][1];
     N[23] = n[1][1] * n[2][2];
     N[24] = n[1][0] * n[2][0];
     N[25] = n[1][1] * n[2][1];
     N[26] = n[1][2] * n[2][2];

    Matrix[0] = k[0] * N[0] + k[1] * N[3] + k[2] * N[6];
    Matrix[1] = k[0] * N[1] + k[1] * N[4] + k[2] * N[7];
    Matrix[2] = k[0] * N[9] + k[1] * N[10] + k[2] * N[11];
    Matrix[3] = k[0] * N[2] + k[1] * N[5] + k[2] * N[8];
    Matrix[4] = k[0] * N[12] + k[1] * N[13] + k[2] * N[14];
    Matrix[5] = k[0] * N[15] + k[1] * N[16] + k[2] * N[17];
    Matrix[6] = k[3] * N[0] + k[4] * N[3] + k[5] * N[6];
    Matrix[7] = k[7] * N[18] - k[6] * N[19];
    Matrix[8] = k[7] * N[20] - k[6] * N[21];
    Matrix[9] = k[7] * N[24] - k[6] * N[24];
    Matrix[10] = k[3] * N[1] + k[4] * N[4] + k[5] * N[7];
    Matrix[11] = k[3] * N[9] + k[4] * N[10] + k[5] * N[11];
    Matrix[12] = k[7] * N[22] - k[6] * N[23];
    Matrix[13] = k[7] * N[25] - k[6] * N[25];
    Matrix[14] = k[7] * N[21] - k[6] * N[20];
    Matrix[15] = k[3] * N[2] + k[4] * N[5] + k[5] * N[8];
    Matrix[16] = k[3] * N[12] + k[4] * N[13] + k[5] * N[14];
    Matrix[17] = k[3] * N[15] + k[4] * N[16] + k[5] * N[17];
    Matrix[18] = k[7] * N[26] - k[6] * N[26];
    Matrix[19] = k[7] * N[23] - k[6] * N[22];
    Matrix[20] = k[7] * N[19] - k[6] * N[18];
    Matrix[21] = k[0] * N[0] + k[1] * N[3] + k[2] * N[6];
    Matrix[22] = k[6] * N[18] - k[7] * N[19];
    Matrix[23] = k[6] * N[20] - k[7] * N[21];
    Matrix[24] = k[6] * N[24] - k[7] * N[24];
    Matrix[25] = -k[0] * N[15] - k[1] * N[16] - k[2] * N[17];
    Matrix[26] = -k[0] * N[9] - k[1] * N[10] - k[2] * N[11];
    Matrix[27] = -k[0] * N[0] - k[1] * N[3] - k[2] * N[6];
    Matrix[28] = k[0] * N[1] + k[1] * N[4] + k[2] * N[7];
    Matrix[29] = k[0] * N[9] + k[1] * N[10] + k[2] * N[11];
    Matrix[30] = k[6] * N[22] - k[7] * N[23];
    Matrix[31] = k[6] * N[25] - k[7] * N[25];
    Matrix[32] = k[6] * N[21] - k[7] * N[20];
    Matrix[33] = -k[0] * N[12] - k[1] * N[13] - k[2] * N[14];
    Matrix[34] = -k[0] * N[1] - k[1] * N[4] - k[2] * N[7];
    Matrix[35] = -k[0] * N[9] - k[1] * N[10] - k[2] * N[11];
    Matrix[36] = k[0] * N[2] + k[1] * N[5] + k[2] * N[8];
    Matrix[37] = k[0] * N[12] + k[1] * N[13] + k[2] * N[14];
    Matrix[38] = k[0] * N[15] + k[1] * N[16] + k[2] * N[17];
    Matrix[39] = k[6] * N[26] - k[7] * N[26];
    Matrix[40] = k[6] * N[23] - k[7] * N[22];
    Matrix[41] = k[6] * N[19] - k[7] * N[18];
    Matrix[42] = -k[0] * N[2] - k[1] * N[5] - k[2] * N[8];
    Matrix[43] = -k[0] * N[12] - k[1] * N[13] - k[2] * N[14];
    Matrix[44] = -k[0] * N[15] - k[1] * N[16] - k[2] * N[17];
    Matrix[45] = k[3] * N[0] + k[4] * N[3] + k[5] * N[6];
    Matrix[46] = k[6] * N[19] - k[7] * N[18];
    Matrix[47] = k[6] * N[21] - k[7] * N[20];
    Matrix[48] = k[6] * N[24] - k[7] * N[24];
    Matrix[49] = (k[4] * N[16]) / 2.0 - k[3] * N[15] + (k[5] * N[17]) / 2.0;
    Matrix[50] = (k[4] * N[10]) / 2.0 - k[3] * N[9] + (k[5] * N[11]) / 2.0;
    Matrix[51] = -k[3] * N[0] + (k[4] * N[3]) / 2.0 + (k[5] * N[6]) / 2.0;
    Matrix[52] = k[7] * N[18] - k[6] * N[19];
    Matrix[53] = k[7] * N[20] - k[6] * N[21];
    Matrix[54] = k[7] * N[24] - k[6] * N[24];
    Matrix[55] = k[3] * N[1] + k[4] * N[4] + k[5] * N[7];
    Matrix[56] = k[3] * N[9] + k[4] * N[10] + k[5] * N[11];
    Matrix[57] = k[6] * N[23] - k[7] * N[22];
    Matrix[58] = k[6] * N[25] - k[7] * N[25];
    Matrix[59] = k[6] * N[20] - k[7] * N[21];
    Matrix[60] = (k[4] * N[13]) / 2.0 - k[3] * N[12] + (k[5] * N[14]) / 2.0;
    Matrix[61] = -k[3] * N[1] + (k[4] * N[4]) / 2.0 + (k[5] * N[7]) / 2.0;
    Matrix[62] = (k[4] * N[10]) / 2.0 - k[3] * N[9] + (k[5] * N[11]) / 2.0;
    Matrix[63] = k[7] * N[22] - k[6] * N[23];
    Matrix[64] = k[7] * N[25] - k[6] * N[25];
    Matrix[65] = k[7] * N[21] - k[6] * N[20];
    Matrix[66] = k[3] * N[2] + k[4] * N[5] + k[5] * N[8];
    Matrix[67] = k[3] * N[12] + k[4] * N[13] + k[5] * N[14];
    Matrix[68] = k[3] * N[15] + k[4] * N[16] + k[5] * N[17];
    Matrix[69] = k[6] * N[26] - k[7] * N[26];
    Matrix[70] = k[6] * N[22] - k[7] * N[23];
    Matrix[71] = k[6] * N[18] - k[7] * N[19];
    Matrix[72] = -k[3] * N[2] + (k[4] * N[5]) / 2.0 + (k[5] * N[8]) / 2.0;
    Matrix[73] = (k[4] * N[13]) / 2.0 - k[3] * N[12] + (k[5] * N[14]) / 2.0;
    Matrix[74] = (k[4] * N[16]) / 2.0 - k[3] * N[15] + (k[5] * N[17]) / 2.0;
    Matrix[75] = k[7] * N[26] - k[6] * N[26];
    Matrix[76] = k[7] * N[23] - k[6] * N[22];
    Matrix[77] = k[7] * N[19] - k[6] * N[18];

	
/*
	for (int i = 0; i < 78; i++)
	{
		cout << "Matrix[["<< i <<"]=" <<Matrix[i]<<endl;
	}
*/

}

void CBeam::ElementStress(double* stress, double* Displacement)
{
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
	clear(stress,3);
    double DX[3]; //	
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
    {
        S[i] = -DX[i] * DX[i] * material_->E / (L * L * L);
        S[i + 3] = -S[i];
    }


    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (LocationMatrix_[i * 6 + j])
            {
                stress[j] += S[i * 3 + j] * Displacement[LocationMatrix_[i * 6 + j] - 1];
            }
        }
    }
}

void CBeam::ElementPostInfo(double* beamstress, double* Displacement, double* prePositionBeam, double* postPositionBeam)
{
	CBeamMaterial* material_ = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
    double DX[3]; //	
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

	 double n[3][3];
	 n[0][0] = DX[0] / L;
	 n[0][1] = DX[1] / L;
	 n[0][2] = DX[2] / L;
	 n[1][0] = material_->n1;
     n[1][1] = material_->n2;
     n[1][2] = material_->n3; 
     n[2][0] = n[0][1] * n[1][2] - n[0][2] * n[1][1];
     n[2][1] = n[0][2] * n[1][0] - n[0][0] * n[1][2];
     n[2][2] = n[0][0] * n[1][1] - n[0][1] * n[1][0]; 
	
    double Loc[2][3]; // preposition
    double d[2][3]; //displacement in the main coordinate
    double D[2][3]; //displacement in the local coordinate
    double theta[2][3]; //rotation in the main coordinate
    double phi[2][3]; //rotation in the local coordinate
    double r[4][3];//vector from center of rectangle to four angles in the local coordinate
    double R[4][3];//vector from center of rectangle to four angles in the main coordinate
	double a=material_->a;
	double b=material_->b;
	double x1=material_->x1;
	double y1=material_->y1;
	double x2=material_->x2;
	double y2=material_->y2;
	double Iz = a*a*a*b/12-(a-x1-x2)*(a-x1-x2)*(a-x1-x2)*(b-y1-y2)/12;
	double Iy = b*b*b*a/12-(b-y1-y2)*(b-y1-y2)*(b-y1-y2)*(a-x1-x2)/12;


	// Define the scale of co-dimension
	double magCodim = 0.1;
    r[0][0] = 0.0;
    r[1][0] = 0.0;
    r[2][0] = 0.0;
    r[3][0] = 0.0;
    r[0][1] = - magCodim * a;
    r[1][1] = magCodim * a;
    r[2][1] = magCodim * a;
    r[3][1] = -magCodim * a;
    r[0][2] = magCodim * b;
    r[1][2] = magCodim * b;
    r[2][2] = - magCodim * b;
    r[3][2] = - magCodim * b;

    for (unsigned int i = 0; i < 4; i++){
        R[i][0] = n[0][0] * r[i][0] + n[1][0] * r[i][1] + n[2][0] * r[i][2];
        R[i][1] = n[0][1] * r[i][0] + n[1][1] * r[i][1] + n[2][1] * r[i][2];
        R[i][2] = n[0][2] * r[i][0] + n[1][2] * r[i][1] + n[2][2] * r[i][2];
    }    

	for (unsigned int i = 0; i < 3; i++){

		if (LocationMatrix_[i]){
		    d[0][i] = Displacement[LocationMatrix_[i]-1];
            Loc[0][i] = nodes_[0]->XYZ[i];
		 }
		else{
            d[0][i] = 0.0;
		    Loc[0][i] = nodes_[0]->XYZ[i];	 
		}

		if (LocationMatrix_[i+6]){
	        d[1][i] = Displacement[LocationMatrix_[i+6]-1];
            Loc[1][i] = nodes_[1]->XYZ[i];
		}
		else{
            d[1][i] = 0.0;
		    Loc[1][i] = nodes_[1]->XYZ[i];
		}
	}
	
	for (unsigned int i = 3; i < 6; i++){

		if (LocationMatrix_[i]){
		  theta[0][i-3] = Displacement[LocationMatrix_[i]-1];
		 }
		else{
		  theta[0][i-3] = nodes_[0]->XYZ[i-3];	 
		}

		if (LocationMatrix_[i+6]){
		  theta[1][i-3] = Displacement[LocationMatrix_[i+6]-1];
		}
		else{		 
		  theta[1][i-3] = nodes_[1]->XYZ[i-3];
		}
	}

    for (unsigned int i = 0; i < 2; i++){
        for (unsigned int j = 0; j < 4; j++){
            postPositionBeam[(i * 4 + j) * 3] = Loc[i][0] + d[i][0] + R[j][0] + R[j][2] * theta[i][1] - R[j][1] * theta[i][2];
            postPositionBeam[(i * 4 + j) * 3 + 1] = Loc[i][1] + d[i][1] + R[j][1] + R[j][0] * theta[i][2] - R[j][2] * theta[i][0];
            postPositionBeam[(i * 4 + j) * 3 + 2] = Loc[i][2] + d[i][2] + R[j][2] + R[j][1] * theta[i][0] - R[j][0] * theta[i][1];
            prePositionBeam[(i * 4 + j) * 3] = Loc[i][0] + R[j][0];
            prePositionBeam[(i * 4 + j) * 3 + 1] = Loc[i][1] + R[j][1];
            prePositionBeam[(i * 4 + j) * 3 + 2] = Loc[i][2] + R[j][2];
        }
    }
    
    //coordinate conversion
    for (unsigned int i = 0; i < 2; i++){
        D[i][0] = n[0][0] * d[i][0] + n[1][0] * d[i][1] + n[2][0] * d[i][2];
        D[i][1] = n[0][1] * d[i][0] + n[1][1] * d[i][1] + n[2][1] * d[i][2];
        D[i][2] = n[0][2] * d[i][0] + n[1][2] * d[i][1] + n[2][2] * d[i][2];
        phi[i][0] = n[0][0] * theta[i][0] + n[1][0] * theta[i][1] + n[2][0] * theta[i][2];
        phi[i][1] = n[0][1] * theta[i][0] + n[1][1] * theta[i][1] + n[2][1] * theta[i][2];
        phi[i][2] = n[0][2] * theta[i][0] + n[1][2] * theta[i][1] + n[2][2] * theta[i][2];
    }

    double dtheta[2][3];//rate of the change of the corner
    dtheta[0][0] = (phi[1][2] - phi[0][2]) / L;
    dtheta[1][0] = dtheta[0][0];
	dtheta[0][1] = (D[0][2] - D[1][2]) * 6.0 / (L*L) - phi[0][1] * 4.0 / L + phi[1][1] * 2.0 / L;
    dtheta[1][1] = (D[1][2] - D[0][2]) * 6.0 / (L*L) + phi[0][1] * 2.0 / L - phi[1][1] * 4.0 / L;
    dtheta[0][2] = (D[1][1] - D[0][1]) * 6.0 / (L*L) - phi[0][2] * 4.0 / L + phi[1][2] * 2.0 / L;
    dtheta[1][2] = (D[0][1] - D[1][1]) * 6.0 / (L*L) + phi[0][2] * 2.0 / L - phi[1][2] * 4.0 / L;

    double sigma1; //Normal stress caused by strech
    double sigma2[2][2]; //Normal stress caused by bending(z-bending and y-bending)
    double tau_xy;
    double tau_xz;

    sigma1 = material_->E * (D[1][0] - D[0][0]) / L;
    tau_xy = material_->E * b * dtheta[0][0] / (4.0 + 4.0 * material_->Nu);
    tau_xz = material_->E * a * dtheta[0][0] / (4.0 + 4.0 * material_->Nu);
    for (unsigned int i = 0; i < 2; i++){
        sigma2[i][0] = material_->E * a * dtheta[i][2] / 2.0;
        sigma2[i][1] = material_->E * b * dtheta[i][1] / 2.0;
        
        beamstress[i * 24] = sigma1 - sigma2[i][0] + sigma2[i][1];
        beamstress[i * 24 + 1] = 0.0;
        beamstress[i * 24 + 2] = 0.0;
        beamstress[i * 24 + 3] = -tau_xy;
        beamstress[i * 24 + 4] = 0.0;
        beamstress[i * 24 + 5] = tau_xz;
        beamstress[i * 24 + 6] = sigma1 - sigma2[i][0] - sigma2[i][1];
        beamstress[i * 24 + 7] = 0.0;
        beamstress[i * 24 + 8] = 0.0;
        beamstress[i * 24 + 9] = tau_xy;
        beamstress[i * 24 + 10] = 0.0;
        beamstress[i * 24 + 11] = tau_xz;
        beamstress[i * 24 + 12] = sigma1 + sigma2[i][0] - sigma2[i][1];
        beamstress[i * 24 + 13] = 0.0;
        beamstress[i * 24 + 14] = 0.0;
        beamstress[i * 24 + 15] = tau_xy;
        beamstress[i * 24 + 16] = 0.0;
        beamstress[i * 24 + 17] = -tau_xz;
        beamstress[i * 24 + 18] = sigma1 + sigma2[i][0] + sigma2[i][1];
        beamstress[i * 24 + 19] = 0.0;
        beamstress[i * 24 + 20] = 0.0;
        beamstress[i * 24 + 21] = -tau_xy;
        beamstress[i * 24 + 22] = 0.0;
        beamstress[i * 24 + 23] = -tau_xz;   

    }
}

double CBeam::Gravity()
	{return 0;}