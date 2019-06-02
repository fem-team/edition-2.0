#include "PostOutputter.h"
#include "Domain.h"

#define Datalength 14
#define coeff 100

CPostOutputter* CPostOutputter::_instance = nullptr;

CPostOutputter::CPostOutputter(string FileName)
{
    OutputFile.open(FileName);

    if (!OutputFile)
    {
        cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
        exit(1);
    }
}

CPostOutputter* CPostOutputter::Instance(string FileName)
{
    if (!_instance)
        _instance = new CPostOutputter(FileName);
    return _instance;
}

void CPostOutputter::OutputElementStress()
{
	CDomain* FEMData = CDomain::Instance();

    double* Displacement = FEMData->GetDisplacement();

    const unsigned int NUMEG = FEMData->GetNUMEG(); // Number of element groups

    *this << "TITLE = \" STAPpp FEM \" " << endl
          << "VARIABLES = \"X_POST\", \"Y_POST\", \"Z_POST\", "
             "\"STRESS_I\", \"STRESS_II\", \"STRESS_III\", \"STRESS_VONMISES\", "
             "\"STRESS_XX\", \"STRESS_YY\", \"STRESS_ZZ\", \"STRESS_XY\", \"STRESS_YZ\", \"STRESS_ZX\""
          << endl;


	    for (unsigned int EleGrpIndex = 0; EleGrpIndex < NUMEG; EleGrpIndex++)
		{

        // Get the ElementGroup and related infos
        CElementGroup& EleGrp = FEMData->GetEleGrpList()[EleGrpIndex];

        ElementTypes ElementType = EleGrp.GetElementType();	// ElementType
        unsigned int NUME = EleGrp.GetNUME();				// Number of elements
		unsigned int NUMNP = FEMData->GetNUMNP();

		switch (ElementType)
		{

		case ElementTypes::Bar:

			 *this << "ZONE T = \"Bridge\", N = " << NUME * 8 << ", E = " << NUME
                  << ", F = FEPOINT , ET = BRICK, C = RED" << endl;

			double PrePositionBar[24];
            double PostPositionBar[24];
            double stressBar[48];
            double cmptStressBar[4];  // cmptStressBar = {stressI, stressII, stressIII, stress_vonMises};

			 for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
				CElement& Element = EleGrp[Ele];
                Element.ElementPostInfo(stressBar, Displacement, PrePositionBar, PostPositionBar);

                for (unsigned nodeIndex = 0; nodeIndex < 8; nodeIndex++)
                {
                    for (unsigned DegOF = 0; DegOF < 3; DegOF++)
                    {
						*this << setw(Datalength)
							  << (1 - coeff) * PrePositionBar[3 * nodeIndex + DegOF] +
                                     coeff * PostPositionBar[3 * nodeIndex + DegOF];
                    }

                    cmptStressBar[0] = stressBar[6 * nodeIndex] + stressBar[6 * nodeIndex + 1] + stressBar[6 * nodeIndex + 2];
                    cmptStressBar[1] = stressBar[6 * nodeIndex]*stressBar[6 * nodeIndex + 1] - stressBar[6 * nodeIndex + 3]*stressBar[6 * nodeIndex + 3]
                                     + stressBar[6 * nodeIndex]*stressBar[6 * nodeIndex + 2] - stressBar[6 * nodeIndex + 5]*stressBar[6 * nodeIndex + 5]
                                     + stressBar[6 * nodeIndex + 1]*stressBar[6 * nodeIndex + 2] - stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex + 4];
                    cmptStressBar[2] = stressBar[6 * nodeIndex]*stressBar[6 * nodeIndex + 1]*stressBar[6 * nodeIndex + 2]
                                     + stressBar[6 * nodeIndex + 3]*stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex + 5]*2
                                     - stressBar[6 * nodeIndex + 1]*stressBar[6 * nodeIndex + 5]*stressBar[6 * nodeIndex + 5]
                                     - stressBar[6 * nodeIndex + 2]*stressBar[6 * nodeIndex + 3]*stressBar[6 * nodeIndex + 3]
                                     - stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex + 4]*stressBar[6 * nodeIndex];
                    cmptStressBar[3] = sqrt(cmptStressBar[0]*cmptStressBar[0] - cmptStressBar[1]);
					*this << setw(Datalength) << cmptStressBar[0]
						  << setw(Datalength) << cmptStressBar[1]
						  << setw(Datalength) << cmptStressBar[2]
						  << setw(Datalength) << cmptStressBar[3];

					for (unsigned DegOF = 0; DegOF < 6; DegOF++)
                    {
						*this << setw(Datalength) << stressBar[6 * nodeIndex + DegOF];
                    }
                    *this << endl;
                }
            }
            // Node numbers corresponding to each element
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned NumEleNode = 0; NumEleNode < 8; NumEleNode++)
                {
					*this << setw(Datalength) << Ele * 8 + NumEleNode + 1;
                }
                *this << endl;
            }
            break;

		case ElementTypes::Q4:

			*this << "ZONE T = \"Bridge\", N = " << NUME * 4 << ",E = " << NUME
                  << " ,F = FEPOINT , ET = QUADRILATERAL, C = RED" << endl;

            double stress4Q[24];
            double PrePosition4Q[12];
            double PostPosition4Q[12];
            double cmptStress4Q[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
				CElement& Element = EleGrp[Ele];
                Element.ElementPostInfo(stress4Q, Displacement, PrePosition4Q, PostPosition4Q);

                for (unsigned ni = 0; ni < 4; ++ni)
                {
                    for (unsigned dof = 0; dof < 3; ++dof)
						*this << setw(Datalength) << PostPosition4Q[ni * 3 + dof];

                    cmptStress4Q[0] = stress4Q[6 * ni] + stress4Q[6 * ni + 1] + stress4Q[6 * ni + 2];
                    cmptStress4Q[1] = stress4Q[6 * ni]*stress4Q[6 * ni + 1] - stress4Q[6 * ni + 3]*stress4Q[6 * ni + 3]
                                     + stress4Q[6 * ni]*stress4Q[6 * ni + 2] - stress4Q[6 * ni + 5]*stress4Q[6 * ni + 5]
                                     + stress4Q[6 * ni + 1]*stress4Q[6 * ni + 2] - stress4Q[6 * ni + 4]*stress4Q[6 * ni + 4];
                    cmptStress4Q[2] = stress4Q[6 * ni]*stress4Q[6 * ni + 1]*stress4Q[6 * ni + 2]
                                     + stress4Q[6 * ni + 3]*stress4Q[6 * ni + 4]*stress4Q[6 * ni + 5]*2
                                     - stress4Q[6 * ni + 1]*stress4Q[6 * ni + 5]*stress4Q[6 * ni + 5]
                                     - stress4Q[6 * ni + 2]*stress4Q[6 * ni + 3]*stress4Q[6 * ni + 3]
                                     - stress4Q[6 * ni + 4]*stress4Q[6 * ni + 4]*stress4Q[6 * ni];
                    cmptStress4Q[3] = sqrt(cmptStress4Q[0]*cmptStress4Q[0] - cmptStress4Q[1]);
					*this << setw(Datalength) << cmptStress4Q[0]
					<< setw(Datalength) << cmptStress4Q[1]
					<< setw(Datalength) << cmptStress4Q[2]
					<< setw(Datalength) << cmptStress4Q[3];

                    for (unsigned dof = 0; dof < 6; ++dof)
						*this << setw(Datalength) << stress4Q[ni * 6 + dof];
                    *this << std::endl;
                }
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned NumEleNode = 0; NumEleNode < 4; NumEleNode++)
					*this << setw(Datalength) << Ele * 4 + NumEleNode + 1;
                
                *this << endl;
            }
			
			*this << endl;
			break;
			 
		case ElementTypes::Beam:

				*this << "ZONE T = \"Bridge\", N = " << NUME * 8 << ",E = " << NUME
                  << " ,F = FEPOINT , ET = BRICK, C = RED" << endl;

            double beamstress[48];
            double prePositionBeam[24];
            double postPositionBeam[24];
            double cmptStressBeam[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
				CElement& Element = EleGrp[Ele];
                Element.ElementPostInfo(beamstress, Displacement, prePositionBeam,postPositionBeam);

                for (unsigned i = 0; i < 8; i++)
                {

                    for (unsigned DegOF = 0; DegOF < 3; DegOF++)
						*this << setw(Datalength)
                              << (1 - coeff) * prePositionBeam[3 * i + DegOF] + coeff * postPositionBeam[3 * i + DegOF];
                    

                    cmptStressBeam[0] = beamstress[6 * i] + beamstress[6 * i + 1] + beamstress[6 * i + 2];
                    cmptStressBeam[1] = beamstress[6 * i]*beamstress[6 * i + 1] - beamstress[6 * i + 3]*beamstress[6 * i + 3]
                                     + beamstress[6 * i]*beamstress[6 * i + 2] - beamstress[6 * i + 5]*beamstress[6 * i + 5]
                                     + beamstress[6 * i + 1]*beamstress[6 * i + 2] - beamstress[6 * i + 4]*beamstress[6 * i + 4];
                    cmptStressBeam[2] = beamstress[6 * i]*beamstress[6 * i + 1]*beamstress[6 * i + 2]
                                     + beamstress[6 * i + 3]*beamstress[6 * i + 4]*beamstress[6 * i + 5]*2
                                     - beamstress[6 * i + 1]*beamstress[6 * i + 5]*beamstress[6 * i + 5]
                                     - beamstress[6 * i + 2]*beamstress[6 * i + 3]*beamstress[6 * i + 3]
                                     - beamstress[6 * i + 4]*beamstress[6 * i + 4]*beamstress[6 * i];
                    cmptStressBeam[3] = sqrt(cmptStressBeam[0]*cmptStressBeam[0] - cmptStressBeam[1]);
					*this << setw(Datalength) << cmptStressBeam[0]
					<< setw(Datalength) << cmptStressBeam[1]
					<< setw(Datalength) << cmptStressBeam[2]
					<< setw(Datalength) << cmptStressBeam[3];

                    for (unsigned DegOF = 0; DegOF < 6; DegOF++)
						*this << setw(Datalength) << beamstress[6 * i + DegOF];
                    

                    *this << endl;
                }
            }

			for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned NumEleNode = 0; NumEleNode < 8; NumEleNode++)
					*this << setw(Datalength) << Ele * 8 + NumEleNode + 1;
                
                *this << endl;
            }
            *this << endl;
			 
			break;

		case ElementTypes::T3: // 3T element

            *this << "ZONE T = \"Bridge\", N = " << NUME * 3 << ",E = " << NUME
                  << " ,F = FEPOINT , ET = TRIANGLE, C = RED" << endl;

            double stress3T[3];
            double PrePosition3T[9];
            double PostPosition3T[9];
            double cmptStress3T[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp[Ele];
                Element.ElementPostInfo(stress3T, Displacement, PrePosition3T, PostPosition3T);
                C3TMaterial material =
                    *dynamic_cast<C3TMaterial*>(Element.GetElementMaterial());

                for (unsigned nodeIndex = 0; nodeIndex < 3; nodeIndex++)
                {
                    for (unsigned dof = 0; dof < 3; dof++)
                        *this << setw(Datalength) << PostPosition3T[nodeIndex * 3 + dof];

                    cmptStress3T[0] = stress3T[6 * nodeIndex] + stress3T[6 * nodeIndex + 1];
                    cmptStress3T[1] = stress3T[6 * nodeIndex]*stress3T[6 * nodeIndex + 1] - stress3T[6 * nodeIndex + 2]*stress3T[6 * nodeIndex + 2];
					cmptStress3T[2] = sqrt(cmptStress3T[0] * cmptStress3T[0] - cmptStress3T[1]);
                    cmptStress3T[3] = 0.0;
                    *this << setw(Datalength) << cmptStress3T[0]
                          << setw(Datalength) << cmptStress3T[1]
                          << setw(Datalength) << cmptStress3T[2]
                          << setw(Datalength) << cmptStress3T[3];

                    *this << setw(Datalength) << stress3T[0] << setw(Datalength)
                          << stress3T[1] << setw(Datalength) << 0.0
                          << setw(Datalength) << stress3T[2] << setw(Datalength)
                          << 0.0 << setw(Datalength) << 0.0 << std::endl;
                }
            }
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned _ = 0; _ < 3; _++)
                    *this << setw(Datalength) << Ele * 3 + _ + 1;
                *this << std::endl;
            }
            *this << endl;

            break;

        case ElementTypes::H8: // 8H element
        {
            *this << "ZONE T= \"Bridge\", N = " << NUME * 8 << " ,E = " << NUME
                  << " ,F = FEPOINT , ET = BRICK, C = RED" << endl;


            double*  stressHex = new double[NUME*48];
            double*  PrePosition8H = new double[NUME*24];
            double*  Position8H = new double[NUME*24];
            double cmptStress8H[4];            


			for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp[Ele];
                Element.ElementPostInfo( &stressHex[48*Ele], Displacement, &PrePosition8H[24*Ele], &Position8H[24*Ele]);
            }

           for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned _ = 0; _ < 8; _++)
                {
                    *this << setw(Datalength)
                          << PrePosition8H[24*Ele+_*3 + 0]+coeff*(Position8H[24*Ele+_ * 3 + 0]-PrePosition8H[24*Ele+_*3 + 0])
                          << setw(Datalength)
                          << PrePosition8H[24*Ele+_*3 + 1]+coeff*(Position8H[24*Ele+_ * 3 + 1]-PrePosition8H[24*Ele+_*3 + 1])
                          << setw(Datalength)
                          << PrePosition8H[24*Ele+_*3 + 2]+coeff*(Position8H[24*Ele+_ * 3 + 2]-PrePosition8H[24*Ele+_*3 + 2]);

                    cmptStress8H[0] = stressHex[48*Ele+6 * _] + stressHex[48*Ele+6 * _ + 1] + stressHex[48*Ele+6 * _ + 2];
                    cmptStress8H[1] = stressHex[48*Ele+6 * _]*stressHex[48*Ele+6 * _ + 1] - stressHex[48*Ele+6 * _ + 3]*stressHex[48*Ele+6 * _ + 3]
                                     + stressHex[48*Ele+6 * _]*stressHex[48*Ele+6 * _ + 2] - stressHex[48*Ele+6 * _ + 5]*stressHex[48*Ele+6 * _ + 5]
                                     + stressHex[48*Ele+6 * _ + 1]*stressHex[48*Ele+6 * _ + 2] - stressHex[48*Ele+6 * _ + 4]*stressHex[48*Ele+6 * _ + 4];
                    cmptStress8H[2] = stressHex[48*Ele+6 * _]*stressHex[48*Ele+6 * _ + 1]*stressHex[48*Ele+6 * _ + 2]
                                     + stressHex[48*Ele+6 * _ + 3]*stressHex[48*Ele+6 * _ + 4]*stressHex[48*Ele+6 * _ + 5]*2
                                     - stressHex[48*Ele+6 * _ + 1]*stressHex[48*Ele+6 * _ + 5]*stressHex[48*Ele+6 * _ + 5]
                                     - stressHex[48*Ele+6 * _ + 2]*stressHex[48*Ele+6 * _ + 3]*stressHex[48*Ele+6 * _ + 3]
                                     - stressHex[48*Ele+6 * _ + 4]*stressHex[48*Ele+6 * _ + 4]*stressHex[48*Ele+6 * _];
                    cmptStress8H[3] = sqrt(cmptStress8H[0]*cmptStress8H[0] - cmptStress8H[1]);
			
                    *this << setw(Datalength) << cmptStress8H[0]
                          << setw(Datalength) << cmptStress8H[1]
                          << setw(Datalength) << cmptStress8H[2]
                          << setw(Datalength) << cmptStress8H[3];

					*this << setw(Datalength) << stressHex[48*Ele+_*6 + 0]
                          << setw(Datalength) << stressHex[48*Ele+_*6 + 1]
                          << setw(Datalength) << stressHex[48*Ele+_*6 + 2]
                          << setw(Datalength) << stressHex[48*Ele+_*6 + 3]
                          << setw(Datalength) << stressHex[48*Ele+_*6 + 4]
                          << setw(Datalength) << stressHex[48*Ele+_*6 + 5]
                          << endl;
                }
            }

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned int i = 1; i < 9; i++)
                {
                    *this << setw(Datalength) << Ele * 8 + i;
                }
                *this << endl;
            }

			delete stressHex;
			delete PrePosition8H;
			delete Position8H;
		}
            break;

		case ElementTypes::Plate:
            *this << "ZONE T = \"Bridge\", N = " << 8 * NUME << " E = " << NUME
                  << " F = FEPOINT , ET = BRICK, C = RED" << endl;

            double stresses4PE[48];
            double PrePositions4PE[24];
            double Positions4PE[24];
            double cmptStressPlate[4];

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                EleGrp[Ele].ElementPostInfo(stresses4PE, Displacement, PrePositions4PE,
                                                       Positions4PE);
                for (unsigned i = 0; i < 4; ++i)
                { // four gauss points
					*this << setw(Datalength)
						<< (1 - coeff) * PrePositions4PE[3 * i] + coeff * Positions4PE[3 * i]
						<< setw(Datalength)
						<< (1 - coeff) * PrePositions4PE[3 * i + 1] +
						coeff * Positions4PE[3 * i + 1]
						<< setw(Datalength)
						<< (1 - coeff) * PrePositions4PE[3 * i + 2] +
						coeff * Positions4PE[3 * i + 2];

					cmptStressPlate[0] = stresses4PE[6 * i] + stresses4PE[6 * i + 1] + stresses4PE[6 * i + 2];
					cmptStressPlate[1] = stresses4PE[6 * i] * stresses4PE[6 * i + 1] - stresses4PE[6 * i + 3] * stresses4PE[6 * i + 3]
						+ stresses4PE[6 * i] * stresses4PE[6 * i + 2] - stresses4PE[6 * i + 5] * stresses4PE[6 * i + 5]
						+ stresses4PE[6 * i + 1] * stresses4PE[6 * i + 2] - stresses4PE[6 * i + 4] * stresses4PE[6 * i + 4];
                    cmptStressPlate[2] = stresses4PE[6 * i]*stresses4PE[6 * i + 1]*stresses4PE[6 * i + 2]
                                     + stresses4PE[6 * i + 3]*stresses4PE[6 * i + 4]*stresses4PE[6 * i + 5]*2
                                     - stresses4PE[6 * i + 1]*stresses4PE[6 * i + 5]*stresses4PE[6 * i + 5]
                                     - stresses4PE[6 * i + 2]*stresses4PE[6 * i + 3]*stresses4PE[6 * i + 3]
                                     - stresses4PE[6 * i + 4]*stresses4PE[6 * i + 4]*stresses4PE[6 * i];
					cmptStressPlate[2] = sqrt(cmptStressPlate[0] * cmptStressPlate[0] - cmptStressPlate[1]);
					*this << setw(Datalength) << cmptStressPlate[0]
						<< setw(Datalength) << cmptStressPlate[1]
						<< setw(Datalength) << cmptStressPlate[2]
						<< setw(Datalength) << cmptStressPlate[3];
					*this << setw(Datalength) << stresses4PE[6 * i]
						<< setw(Datalength) << stresses4PE[6 * i + 1]
						<< setw(Datalength) << stresses4PE[6 * i + 2]
						<< setw(Datalength) << stresses4PE[6 * i + 3]
						<< setw(Datalength) << stresses4PE[6 * i + 4]
						<< setw(Datalength) << stresses4PE[6 * i + 5] << endl;
                }
            }
            *this << endl;
            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                for (unsigned int i = 0; i < 8; ++i)
                {
                    *this << setw(Datalength) << 8 * Ele + i + 1;
                }
                *this << std::endl;
            }
            break;

		case ElementTypes::Shell:
			*this << "ZONE T = \"Bridge\", N = " << 8 * NUME << " E = " << NUME
				<< " F = FEPOINT , ET = BRICK, C = RED" << endl;

			double stresses4SE[48];
			double PrePositions4SE[24];
			double Positions4SE[24];
			double cmptStressShell[4];

			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				EleGrp[Ele].ElementPostInfo(stresses4SE, Displacement, PrePositions4SE,
					Positions4SE);

				for (unsigned i = 0; i < 8; ++i)
				{
					// four gauss points;
					*this << setw(Datalength)
						<< (1 - coeff) * PrePositions4SE[3 * i] + coeff * Positions4SE[3 * i]
						<< setw(Datalength)
						<< (1 - coeff) * PrePositions4SE[3 * i + 1] +
						coeff * Positions4SE[3 * i + 1]
						<< setw(Datalength)
						<< (1 - coeff) * PrePositions4SE[3 * i + 2] +
						coeff * Positions4SE[3 * i + 2];

					cmptStressShell[0] = stresses4SE[6 * i] + stresses4SE[6 * i + 1] + stresses4SE[6 * i + 2];
					cmptStressShell[1] = stresses4SE[6 * i] * stresses4SE[6 * i + 1] - stresses4SE[6 * i + 3] * stresses4SE[6 * i + 3]
						+ stresses4SE[6 * i] * stresses4SE[6 * i + 2] - stresses4SE[6 * i + 5] * stresses4SE[6 * i + 5]
						+ stresses4SE[6 * i + 1] * stresses4SE[6 * i + 2] - stresses4SE[6 * i + 4] * stresses4SE[6 * i + 4];
					cmptStressShell[2] = stresses4SE[6 * i] * stresses4SE[6 * i + 1] * stresses4SE[6 * i + 2]
						+ stresses4SE[6 * i + 3] * stresses4SE[6 * i + 4] * stresses4SE[6 * i + 5] * 2
						- stresses4SE[6 * i + 1] * stresses4SE[6 * i + 5] * stresses4SE[6 * i + 5]
						- stresses4SE[6 * i + 2] * stresses4SE[6 * i + 3] * stresses4SE[6 * i + 3]
						- stresses4SE[6 * i + 4] * stresses4SE[6 * i + 4] * stresses4SE[6 * i];
					cmptStressShell[3] = sqrt(cmptStressShell[0] * cmptStressShell[0] - cmptStressShell[1]);
					*this << setw(Datalength) << cmptStressShell[0]
						<< setw(Datalength) << cmptStressShell[1]
						<< setw(Datalength) << cmptStressShell[2]
						<< setw(Datalength) << cmptStressShell[3];

					*this << setw(Datalength) << stresses4SE[6 * i]
						<< setw(Datalength) << stresses4SE[6 * i + 1]
						<< setw(Datalength) << stresses4SE[6 * i + 2]
						<< setw(Datalength) << stresses4SE[6 * i + 3]
						<< setw(Datalength) << stresses4SE[6 * i + 4]
						<< setw(Datalength) << stresses4SE[6 * i + 5] << endl;
				}
			}
			for (unsigned int Ele = 0; Ele < NUME; Ele++)
			{
				for (unsigned int i = 0; i < 8; ++i)
				{
					*this << setw(Datalength) << 8 * Ele + i + 1;
				}
				*this << std::endl;
			}

			break;

		default: // Invalid element type
				cerr << "*** Error *** Elment type  " << ElementType
					<< " has not been implemented for PostOutputter.\n\n";

		};

			


		
		}
};