#include "FourMomentum.h"

namespace hsngen {

	// Constructor for FourMomentum
	FourMomentum::FourMomentum(double energy, std::vector<double> momentum)
	{
		Populate(energy, momentum);
	} //  END constructor FourMomentum


	// Constructor for empty FourMomentum
	// (to be filled later with "Populate" method)
	FourMomentum::FourMomentum()
	{
		E = 0.0;
		modp=0;
		p.push_back(0.0);
		p.push_back(0.0);
		p.push_back(0.0);
		mass = 0.0;
	} // END constructor FourMomentum


	// Assign values to FourMomentum from with energy and 3-momentum.
	int FourMomentum::Populate(double energy, std::vector<double> momentum)
	{
		E = energy;
		p = momentum;

		//Just check that the input is OK.
		if(p.size() != 3) { std::cout<<"ERROR: 3-momentum wrong size."<<std::endl; }

		//Define the total 3-momentum and the full Minkowski norm (and check it's an on-shell 4-momentum).
		modp = sqrt(p.at(0)*p.at(0) + p.at(1)*p.at(1) + p.at(2)*p.at(2));

		mass = E*E - modp*modp; //Using mass as a temporary variable here. True value created a few lines down.
		if(fabs(mass) < 1e-12){ mass = 0.0; }

		if(mass < 0.0 ){ std::cout<<"ERROR: 4-vector is spacelike. This isn't what we had agreed on!"<<std::endl; }
		else{ mass = sqrt(mass); }

		return 0;
	} // END function Populate


	// Print the content of a FourMomentum
	void FourMomentum::Print(std::string name)
	{
		std::cout<<"Fourvector '"<<name<<"'"<<" = ("<<E<<", "<<p.at(0)<<", "<<p.at(1)<<", "<<p.at(2)<<"),\t"<<"[Inv. Mass^2: "<<mass<<", Norm of 3-momentum: "<<modp<<"]"<<std::endl;

		return;
	} // END function Print


	// Return direction vector
	std::vector<double> FourMomentum::Direction()
	{
		std::vector<double> temp;

		if(modp==0)
		{
			std::cout<<"Cannot compute direction: 3-momentum vanishes."<<std::endl;
			temp.push_back(0.0);
			temp.push_back(0.0);
		}
		else{
			temp.push_back(acos(p.at(2)/modp)); //acos returns 0 to pi.
			double phi = fabs(atan(p.at(1)/p.at(0)));

			if(p.at(0)>=0.0 && p.at(1) >= 0.0) temp.push_back(phi);
			else if(p.at(0)>=0.0 && p.at(1) < 0.0) temp.push_back(-phi);
			else if(p.at(0)<0.0 && p.at(1) >= 0.0) temp.push_back(M_PI-phi);
			else if(p.at(0)<0.0 && p.at(1) < 0.0)	temp.push_back(-M_PI+phi);
			else printf("Something horrible happened when trying to calculate direction:\nx:%.2f y:%.2f", p.at(0),p.at(1));
		}

		return temp;
	} // END function Direction


	// Return gamma factor
	double FourMomentum::Gamma()
	{
		double temp;
		if(mass==0)
		{
			std::cout<<"ERROR: Trying to compute gamma factor for NULL four momentum."<<std::endl;
			temp = 1e-5;
		}
		else temp = E/mass;

		return temp;
	} // END function Gamma


	// Apply rotation boost from a parent, changing the FourMomentum attributes
	int FourMomentum::RotBoostFromParent(FourMomentum * parentP)
	{
		double costheta = cos(parentP->Direction().at(0));
		double sintheta = sin(parentP->Direction().at(0));
		double phi = parentP->Direction().at(1);

		double gamma = parentP->Gamma();
		double beta = sqrt(1.0-1.0/(gamma*gamma));

		double temp[4];
		temp[0]=E;
		temp[1]=p.at(0);
		temp[2]=p.at(1);
		temp[3]=p.at(2);

		double new_E;
		std::vector<double> new_p;

		new_E = gamma*temp[0] + gamma*beta*temp[3];
		new_p.push_back( gamma*beta*cos(phi)*sintheta*temp[0] + cos(phi)*costheta*temp[1] - sin(phi)*temp[2] + gamma*cos(phi)*sintheta*temp[3] );
		new_p.push_back( gamma*beta*sin(phi)*sintheta*temp[0] + sin(phi)*costheta*temp[1] + cos(phi)*temp[2] + gamma*sin(phi)*sintheta*temp[3] );
		new_p.push_back( gamma*beta*costheta*temp[0] - sintheta*temp[1] + gamma*costheta*temp[3] );

		FourMomentum::Populate(new_E, new_p);

		return 0;
	} // END function RotBoostFromParent

} // END namespace hsngen