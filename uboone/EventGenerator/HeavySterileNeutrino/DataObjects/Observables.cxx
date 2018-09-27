#include "Observables.h"

namespace hsngen {
	
	// Default empty constructor
	Observables::Observables() {}

	// Print output particle in the hepevt format
	void Observables::PrintHepEvt(int i)
	{
		// pdg1 = 11;
		// pdg2 = 211;
		printf("%i 2\n",i);
		printf("1 %i 0 0 0 0 %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",pdg1,P1[0],P1[1],P1[2],E1,mass1,xPos,yPos,zPos,time);
		printf("1 %i 0 0 0 0 %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",pdg2,P2[0],P2[1],P2[2],E2,mass2,xPos,yPos,zPos,time);
	} // END function PrintHepEvt

} // END namespace hsngen
