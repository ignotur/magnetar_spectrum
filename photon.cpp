#include <iostream>
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

class magnetic_field {
	private:
		double p;      // Parameter describing the magnetic field twist
		double C;      // Eigenvalue
		int N = 20;    // Angular resolution for the differential equation
		double mu [N]; // Mesh for cos \theta  
		double F [N];  //
		double Bpole;  // Strength of the magnetic field at the equator
		double beta;   // Average charge velocity in units of c
	public:
		magnetic_field (double p, int N, double beta);
		double Br      (r, theta);
		double Btheta  (r, theta);
		double Bphi    (r, theta);
		double ne      (r, theta);
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


int main () {

	photon pht (0.0, 0.0, 1.0e6, 1.0);
	cout << "Polarisation state: " << pht.get_polarisation_state() << endl;
	cout << "Beaming parameter: "<<   pht.get_beaming_parameter () << endl;
	cout << "mu: " <<                 pht.get_mu () << endl; 
	cout << "omega: " <<              pht.get_omega () << endl;



}
