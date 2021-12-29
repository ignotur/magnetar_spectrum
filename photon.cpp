#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "monte_carlo.h"
#include "planck.h"

using namespace std;

// Code to model spectra of the radiation emitted by a hot NS and scattered in the twisted magnetosphere
// The main scattering mechanism is the resonant cyclotron scattering
// The code is based on article by Nobili, Turolla & Zane (2008)
// Written by Dr. Andrei P. Igoshev ignotur@gmail.com

class magnetosphere {
	private:
		double p;       // Parameter describing the magnetic field twist
		double C;       // Eigenvalue
		int N;          // Angular resolution for the differential equation
		double * mu;    // Pointer for mesh for cos \theta  
		double * F;     //
		double Bpole;   // Strength of the magnetic field at the equator
		double beta;    // Average charge velocity in units of c
		double Rns;     // NS radius cm
		double echarge; // Elementary charge
	public:
		magnetosphere (double Bpole_par, double beta_par, double Rns); // Initialise the twisted magnetosphere configuration by reading file with solved equation
		double Ffun    (double theta);
		double dFfun   (double theta);
		double Br      (double r, double theta) {return -0.5*Bpole * pow( Rns/r , 2.0+p) * dFfun (theta);};
		double Btheta  (double r, double theta) {return  0.5*Bpole * pow(Rns/r  , 2.0+p) * p * Ffun (theta) / sin(theta);};
		double Bphi    (double r, double theta) {return Btheta (r, theta) * sqrt(C / (p*(p+1))) * pow(Ffun(theta), 1.0/p);};
		double B       (double r, double theta) {return sqrt(pow(Br(r, theta), 2.0) + pow(Btheta (r,theta), 2.0) + pow(Bphi(r,theta), 2.0));};
		double ne      (double r, double theta) {return (p+1) / (4.0 * M_PI * echarge) * (Bphi(r,theta) / Btheta (r,theta)) * B (r, theta) / r / beta ;};
};


class photon {
	private:
		double b;       // Beaming parameter
		double omega;   // Photon frequency 
		int s;          // Photon polarisation state: 0 - photon does not exist; 1 - ordinary mode; 2 - extraordinary mode
		double mu;      // Cosine of the angle between the initial photon direction and magnetic field
		double azimuth; // Azimuth angle for photon
	public:
		photon (double theta, double phi, double T, double beaming); // Constructor which corresponds to emission of the photon from the NS surface path with coordinates phi, theta and temperature T
		int    propagate ();
		int    get_polarisation_state () {return s;};
		double get_beaming_parameter  () {return b;};
		double get_mu                 () {return mu;};
		double get_omega              () {return omega;};
};

photon::photon (double theta, double phi, double T, double beaming) {
	b = beaming;

	srand (time (NULL));

	s = rand() % 2 + 1; // Draw the initial polarisation from the uniform distribution

	mu = pow (uniform (0.0, 1.0), 1.0 / b); // cos is uniform taking into account the beaming factor

	azimuth = 2.0 * M_PI * uniform (0.0, 1.0);

	omega = sed_planck (T);


}

magnetosphere::magnetosphere (double Bpole_par, double beta_par, double Rns_par) {

	Bpole   = Bpole_par;
	beta    = beta_par; 
	Rns     = Rns_par;
	echarge = 4.8032e-10;

	ifstream infile ("F_magnetosphere.txt");

	infile >> N;
	infile >> p; 
	infile >> C;

	mu = new double [N];
	F  = new double [N];

	for (int i =0; i < N; i++) {

		infile >> mu[i] >> F[i];

	}

	cout << "We initiliased twisted magnetosphere using file F_magnetosphere.txt" << endl;
	cout << "N = " << N << " \t p = " << p << " \t C = "<< C << endl;

};

double magnetosphere::Ffun (double theta) {

	double res = 0.0;
	double mu_val = cos(theta);

	for (int i=0; i < N; i++) {
		res += F[i] * cos(i*mu_val);
	}

	return res;
}

double magnetosphere::dFfun (double theta) {

	double res = 0.0;
	double mu_val = cos(theta);

	for (int i=0; i < N; i++) {
		res += - i * F[i] * sin(i*mu_val);
	}

	return res;
}


int main () {

	photon pht (0.0, 0.0, 1.0e6, 1.0);
	cout << "Polarisation state: " << pht.get_polarisation_state() << endl;
	cout << "Beaming parameter: "<<   pht.get_beaming_parameter () << endl;
	cout << "mu: " <<                 pht.get_mu () << endl; 
	cout << "omega: " <<              pht.get_omega () << endl;

	magnetosphere (1e14, 0.1, 1e6);



}
