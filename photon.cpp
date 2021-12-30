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
		double p;         // Parameter describing the magnetic field twist
		double C;         // Eigenvalue
		int N;            // Angular resolution for the differential equation
		double * mu;      // Pointer for mesh for cos \theta  
		double * F;       //
		double Bpole;     // Strength of the magnetic field at the equator
		double beta_bulk; // Average charge velocity in units of c
		double Te;        // Electron and positron temperature
		double norm_f;    // Normalisation coeffitient for velocity distribution of charged particles
		double Rns;       // NS radius cm
		double echarge;   // Elementary charge
		double me;        // Electron rest mass
		double speed_of_light; // Speed of light
		double re;        // Classical electron radius
		double kB;        // Boltzman constant
	public:
		magnetosphere (double Bpole_par, double beta_par, double Te_par, double Rns); // Initialise the twisted magnetosphere configuration by reading file with solved equation
		double Ffun    (double theta);
		double dFfun   (double theta);
		double Br      (double r, double theta) {return -0.5*Bpole * pow( Rns/r , 2.0+p) * dFfun (theta);};
		double Btheta  (double r, double theta) {return  0.5*Bpole * pow(Rns/r  , 2.0+p) * p * Ffun (theta) / sin(theta);};
		double Bphi    (double r, double theta) {return Btheta (r, theta) * sqrt(C / (p*(p+1))) * pow(Ffun(theta), 1.0/p);};
		double B       (double r, double theta) {return sqrt(pow(Br(r, theta), 2.0) + pow(Btheta (r,theta), 2.0) + pow(Bphi(r,theta), 2.0));};
		double ne      (double r, double theta) {return (p+1) / (4.0 * M_PI * echarge) * (Bphi(r,theta) / Btheta (r,theta)) * B (r, theta) / r / beta ;};
		double omega_B (double r, double theta) {return echarge * B(r,theta) / (me * speed_of_light);};
		double get_Rns () {return Rns;};
		double f_beta  (double beta_v); // Velocity distribution of charged particles in magnetosphere
};


class photon {
	private:
		double b;        // Beaming parameter
		double omega;    // Photon frequency 
		int s;           // Photon polarisation state: 0 - photon does not exist; 1 - ordinary mode; 2 - extraordinary mode
		double mu;       // Cosine of the angle between the initial photon direction and magnetic field
		double azimuth;  // Azimuth angle for photon
		double pos [3];  // Instantenious location of photon in the Carthesian coordinate system [x,y,z]
		double pos_r [3];// Instantenious location of photon in the spherical coordinate system [r,theta,phi]
		double k [3];    // Instantenious k vector for the photon in the Carthesian coordinate system [kx, ky, kz]
		double kr[3];    // Instantenious k vector for the photon in the spherical coordinate system [kr, ktheta, kphi]
		double c;        // Speed of light
		magnetosphere mg // Object of the twisted magnetosphere class. Physically the magnetosphere through which the photon propagates

	public:
		photon (double theta, double phi, double T, double beaming); // Constructor which corresponds to emission of the photon from the NS surface path with coordinates phi, theta and temperature T
		int    propagate ();
		double propagate_one_step (); // Method to propagate photon one numerical timestep. It returns back the optical depth
		int    scatter ();            // Method to scatter the photon, i.e. it updates omega and k 
		int    get_polarisation_state () {return s;};
		double get_beaming_parameter  () {return b;};
		double get_mu                 () {return mu;};
		double get_omega              () {return omega;};
		double get_mu                 ();
		double beta_plus              ();
		double beta_minus             ();
};

photon::photon (double theta, double phi, double T, double beaming, magnetosphere mag_NS) {
	b = beaming;

	srand (time (NULL));

	s = rand() % 2 + 1; // Draw the initial polarisation from the uniform distribution

	mu = pow (uniform (0.0, 1.0), 1.0 / b); // cos is uniform taking into account the beaming factor

	azimuth = 2.0 * M_PI * uniform (0.0, 1.0);

	omega = sed_planck (T);

	c = 2.99792458e10; // cm/s

	mg = mag_NS;

	pos[0] = mag_NS.get_Rns * cos (theta) * cos (phi); // Check the coordniate system transformation
	pos[1] = mag_NS.get_Rns * cos (theta) * sin (phi);
	pos[2] = mag_NS.get_Rns * sin (theta);


	// two rotations in respect to normal: R(azimuth) * R(acos(mu)) * n
	k[0] = c * cos (theta) * cos (phi); // along the orthogonal direction to NS surface
	k[1] = c * cos (theta) * sin (phi); // tbd
	k[2] = c * sin (theta); 

	cout << "Be careful, position and k vector are not right at the moment" << endl;

}

magnetosphere::magnetosphere (double Bpole_par, double beta_par, double Te_par, double Rns_par) {

	Bpole     = Bpole_par;
	beta_bulk = beta_par; 
	Te        = Te_par;
	Rns       = Rns_par;
	echarge   = 4.8032e-10;
	me        = 9.1094e-28; // g
	speed_of_light = 2.99792458e10; // cm/s
	re        = 2.8179403227e-13; // cm

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

double magnetosphere::f_beta  (double beta_v) {

	double gamma_v, gamma_bulk, Theta_e, gamma_ap;
	double res;
	gamma_v    = 1.0 / (1.0 - beta_v*beta_v);
	gamma_bulk = 1.0 / (1.0 - beta_bulk*beta_bulk);
	Theta_e = kB * Te / (me * speed_of_light);
	gamma_ap = gamma_v * gamma_bulk * (1.0 - beta_v * beta_bulk);

	res = norm_f * exp (-gamma_ap / Theta_e); // eq. (4) in Nobili, Turolla & Zane, normalised numerically 
	
	return res;

}

// Compute the ange in respect to local magnetic field
double photon::get_mu () {

	double B_inst [3];
	double res;

	B_inst[0] = mg.Br     (pos_r[0], pos_r[1]);
	B_inst[1] = mg.Btheta (pos_r[0], pos_r[1]);
	B_inst[2] = mg.Bphi   (pos_r[0], pos_r[1]);

	res = kr[0] * B_inst[0] + kr[1] * B_inst[1] + kr[2] * B_inst[2]; // scalar product 

	return res;
}

// eq. (26) in Fernandez & Thompson (2007) article

double photon::beta_plus () {

	double mu_val, omega_c, res, omegac2;

	mu_val  = get_mu();
	omega_c = mg.omega_B (pos_r[0], pos_r[1]);

	omegac2 = pow(omega_c / omega, 2.0);

	if ((mu_val*mu_val + omegac2) < 1) 
		return -20; // error message - scattering is not possible

	res = 1.0 / (omegac2 + mu_val*mu_val) * (mu_val + omega_c/omega * sqrt(omegac2 + mu_val*mu_val - 1.0));

	return res;
}

// eq. (26) in Fernandez & Thompson (2007) article

double photon::beta_minus () {

	double mu_val, omega_c, res, omegac2;

	mu_val  = get_mu();
	omega_c = mg.omega_B (pos_r[0], pos_r[1]);

	omegac2 = pow(omega_c / omega, 2.0);

	if ((mu_val*mu_val + omegac2) < 1) 
		return -20; // error message - scattering is not possible

	res = 1.0 / (omegac2 + mu_val*mu_val) * (mu_val - omega_c/omega * sqrt(omegac2 + mu_val*mu_val - 1.0));

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
