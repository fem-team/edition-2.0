/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
	bcode[3] = 1;   //将三个rotation自由度默认锁住
	bcode[4] = 1;
	bcode[5] = 1;


};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input, unsigned int np)
{
	unsigned int N;

	Input >> N;	// node number
	if (N != np + 1) 
	{
		cerr << "*** Error *** Nodes must be inputted in order !" << endl 
			 << "   Expected node number : " << np + 1 << endl
			 << "   Provided node number : " << N << endl;

		return false;
	}

	NodeNumber = N;

	string NodeInfo;
	getline(Input,NodeInfo);
	stringstream NodeInfoforInput(NodeInfo);

	if (getArgsNumber(NodeInfo) == 9)
	{
		NodeInfoforInput >> bcode[0] >> bcode[1] >> bcode[2] >> bcode[3] >> bcode[4] >> bcode[5]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];
	}
	else if (getArgsNumber(NodeInfo) == 6)
	{
		NodeInfoforInput >> bcode[0] >> bcode[1] >> bcode[2] 
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];
	}
	else
	{
		cerr << "Please Check the input of NodeInfo! " <<endl
			<< "in this line you have" << getArgsNumber(NodeInfo) << "Args" <<endl;
	}



	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output, unsigned int np)
{
	output << setw(9) << np + 1 << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2] 
		<< setw(5) << bcode[3] << setw(5) << bcode[4] << setw(5) << bcode[5]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output, unsigned int np)
{
	output << setw(9) << np+1 << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, unsigned int np, double* Displacement)
{
	output << setw(5) << np + 1 << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}

int getArgsNumber(std::string buff)
{
    int count = 0;
    bool inNumFlag = false;
    for (unsigned i=0; i<buff.length(); ++i)
    {
        if (buff[i] != ' ' && buff[i] != '\t')
        {
            if (!inNumFlag)
            {
                inNumFlag = true;
                count++;
            }
        }
        else
        {
            inNumFlag = false;
        }
    }
    return count;
}
