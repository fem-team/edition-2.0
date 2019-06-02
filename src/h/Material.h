/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

	double density;

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output, unsigned int mset) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

	//新加 2019.3.24
	double rho_3D;  //体密度
	//    //线密度
	//

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for 3T element
class C3TMaterial : public CMaterial
{
public:

	double Nu;	//!Poisson ratio of a 4Q element

	double thick; // the thickness of element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


class CQ4Material : public CMaterial
{
public:

	double Nu;	//!< Sectional area of a bar element



public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

class CBeamMaterial : public CMaterial
{
public:

	double a;
	double b;
	double Nu;
	double x1;
	double y1;
	double x2;
	double y2;
	double n1;// x component of y' axis
	double n2;// y component of y' axis
	double n3;// z component of y' axis

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);

};
	
class C8HMaterial : public CMaterial
{
public:

	double Nu;	//!< Sectional area of a 8H element

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream OutputFile
	virtual void Write(COutputter& output, unsigned int mset);
};

class CQ5Material : public CMaterial
{
public:

	double nu;	//!< Sectional area of a bar element



public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

class CPlateMaterial : public CMaterial
{
public:

	double nu;   //!< Poisson's ratio

	double h;	//!< Thickness of a plate element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for IEM element
class CIEMMaterial : public CMaterial
{
public:

	double nu;	//!< Sectional area of a bar element



public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//! Material class for Shell element
class CShellMaterial : public CMaterial
{
public:

	double nu; // Poisson's ratio

	double h; // shell thickness


public:
	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

