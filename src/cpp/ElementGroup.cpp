/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Domain.h"

CNode* CElementGroup::NodeList_ = nullptr;

//! Constructor
CElementGroup::CElementGroup()
{
    if (!NodeList_)
    {
        CDomain* FEMData = CDomain::Instance();
        NodeList_ = FEMData->GetNodeList();
    }
    
    ElementType_ = ElementTypes::UNDEFINED;
    
    NUME_ = 0;
    ElementList_ = nullptr;
    
    NUMMAT_ = 0;
    MaterialList_ = nullptr;
}

//! Deconstructor
CElementGroup::~CElementGroup()
{
    if (ElementList_)
        delete [] ElementList_;
    
    if (MaterialList_)
        delete [] MaterialList_;
}

//! operator []
//! For the sake of efficiency, the index bounds are not checked
CElement& CElementGroup::operator[](unsigned int i)
{
    return *(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
}

//! Return index-th material in this element group
CMaterial& CElementGroup::GetMaterial(unsigned int index)
{
    return *(CMaterial*)((std::size_t)(MaterialList_) + index*MaterialSize_);
}

//! Calculate the size of the derived element and material class
void CElementGroup::CalculateMemberSize()
{
    switch (ElementType_)
    {
        case ElementTypes::UNDEFINED:
            std::cerr << "Setting element type to UNDEFINED." << std::endl;
            exit(5);
        case ElementTypes::Bar:
            ElementSize_ = sizeof(CBar);
            MaterialSize_ = sizeof(CBarMaterial);
            break;
		case ElementTypes::T3:
			ElementSize_ = sizeof(C3T);
			MaterialSize_ = sizeof(C3TMaterial);
			break;
		case ElementTypes::H8:
            ElementSize_ = sizeof(C8H);
            MaterialSize_ = sizeof(C8HMaterial);
            break;
		case ElementTypes::Q4:
			ElementSize_ = sizeof(CQ4);
			MaterialSize_ = sizeof(CQ4Material);
			break;
		case ElementTypes::Beam:
			ElementSize_ = sizeof(CBeam);
			MaterialSize_ = sizeof(CBeamMaterial);
			break;
		case ElementTypes::Plate:
			ElementSize_ = sizeof(CPlate);
			MaterialSize_ = sizeof(CPlateMaterial);
			break;
		case ElementTypes::IEM:
			ElementSize_ = sizeof(CIEM);
			MaterialSize_ = sizeof(CIEMMaterial);
			break;
		case ElementTypes::Q5:
			ElementSize_ = sizeof(CQ5);
			MaterialSize_ = sizeof(CQ5Material);
			break;
		case ElementTypes::Shell:
			ElementSize_ = sizeof(CShell);
			MaterialSize_ = sizeof(CShellMaterial);
			break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::CalculateMemberSize." << std::endl;
            exit(5);
            break;
    }
}

//! Allocate array of derived elements
void CElementGroup::AllocateElements(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            ElementList_ = new CBar[size];
            break;
		case ElementTypes::T3:
			ElementList_ = new C3T[size];
			break;
		case ElementTypes::H8:
            ElementList_ = new C8H[size];
            break;
		case ElementTypes::Q4:
			ElementList_ = new CQ4[size];
			break;
		case ElementTypes::Beam:
			ElementList_ = new CBeam[size];
			break;
		case ElementTypes::Plate:
			ElementList_ = new CPlate[size];
			break;
		case ElementTypes::IEM:
			ElementList_ = new CIEM[size];
			break;
		case ElementTypes::Q5:
			ElementList_ = new CQ5[size];
			break;
		case ElementTypes::Shell:
			ElementList_ = new CShell[size];
			break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateElement." << std::endl;
            exit(5);
    }
}

//! Allocate array of derived materials
void CElementGroup::AllocateMaterials(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            MaterialList_ = new CBarMaterial[size];
            break;
		case ElementTypes::T3:
			MaterialList_ = new C3TMaterial[size];
			break;
		case ElementTypes::Q4:
			MaterialList_ = new CQ4Material[size];
			break;
		case ElementTypes::Beam:
			MaterialList_ = new CBeamMaterial[size];
			break;
		case ElementTypes::H8:
			MaterialList_ = new C8HMaterial[size];
			break;	
		case ElementTypes::Plate:
			MaterialList_ = new CPlateMaterial[size];
			break;
		case ElementTypes::IEM:
			MaterialList_ = new CIEMMaterial[size];
			break;
		case ElementTypes::Q5:
			MaterialList_ = new CQ4Material[size];
			break;
		case ElementTypes::Shell:
			MaterialList_ = new CShellMaterial[size];
			break;
        default:
            std::cerr << "Type " << ElementType_ << " not available. See CElementGroup::AllocateMaterial." << std::endl;
            exit(5);
    }
}

//! Read element group data from stream Input
bool CElementGroup::Read(ifstream& Input)
{
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;  // 单元类型、单元数、材料类型
    
    CalculateMemberSize();

    if (!ReadElementData(Input))
        return false;

    return true;
}

//  Read bar element data from the input data file
bool CElementGroup::ReadElementData(ifstream& Input)
{
//  Read material/section property lines
    AllocateMaterials(NUMMAT_);
    
//  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        if (!GetMaterial(mset).Read(Input, mset))
            return false;
    
//  Read element data lines
    AllocateElements(NUME_);
    
//  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
        if (!(*this)[Ele].Read(Input, Ele, MaterialList_, NodeList_))
            return false;
    
    return true;
}
