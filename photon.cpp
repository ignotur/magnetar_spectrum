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
		double ne      (double r, double theta) {return (p+1) / (4.0 * M_PI * echarge) * (Bphi(r,theta) / Btheta (r,theta)) * B (r, theta) / r / beta_bulk;};
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
		magnetosphere *mg;// Object of the twisted magnetosphere class. Physically the magnetosphere through which the photon propagates

	public:
		photon (double theta, double phi, double T, double beaming, magnetosphere); // Constructor which corresponds to emission of the photon from the NS surface path with coordinates phi, theta and temperature T
		int    propagate ();
		double propagate_one_step (double delta_t); // Method to propagate photon one numerical timestep. It returns back the optical depth
		int    scatter ();            // Method to scatter the photon, i.e. it updates omega and k 
		int    get_polarisation_state () {return s;};
		double get_beaming_parameter  () {return b;};
		//double get_mu                 () {return mu;};
		double get_omega              () {return omega;};
		double get_mu                 ();
		double beta_plus              ();
		double beta_minus             ();
		double r                      () {return pos_r[0];};
		double theta                  () {return pos_r[1];};
		double phi                    () {return pos_r[2];};
		double x                      () {return pos[0];};
		double y                      () {return pos[1];};
		double z                      () {return pos[2];};
		void   print_pos              () {cout << "x = "<<pos[0]<<" y = "<<pos[1] <<" z = "<<pos[2] << "; r = "<<pos_r[0]<<" theta = "<<pos_r[1] << " phi = "<<pos_r[2] <<endl; };
		void   print_k                () {cout << "kx = "<<k[0]<<" ky = "<<k[1] <<" kz = "<<k[2]    << "; kr = " <<kr[0] << " k_theta = "<<kr[1] << " k_phi = "<<kr[2] << endl;}
		bool   is_inside_ns           () {if (pos_r[0] < mg->get_Rns()) return true; else return false;};
};

photon::photon (double theta, double phi, double T, double beaming, magnetosphere mag_NS) {

	double theta0, phi0;

	static bool called_ever;

	if (!called_ever) {
		called_ever = true;
		srand (time(NULL));
	}


	b = beaming;

	s = rand() % 2 + 1; // Draw the initial polarisation from the uniform distribution

	mu = pow (uniform (0.0, 1.0), 1.0 / b); // cos is uniform taking into account the beaming factor

	azimuth = 2.0 * M_PI * uniform (0.0, 1.0);

	omega = sed_planck (T);

	c = 2.99792458e10; // cm/s

	mg = &mag_NS;

	pos[0] = mag_NS.get_Rns() * sin (theta) * cos (phi); // Theta is computed from the pole down
	pos[1] = mag_NS.get_Rns() * sin (theta) * sin (phi);
	pos[2] = mag_NS.get_Rns() * cos (theta);

	pos_r[0] = mag_NS.get_Rns();
	pos_r[1] = theta;
	pos_r[2] = phi;

	kr[0] = c;
	kr[1] = acos(mu);
	kr[2] = azimuth;

	theta0 = acos(mu);
	phi0   = azimuth;

	k[0] = c * sin(theta) * cos(phi) + mag_NS.get_Rns() * cos(theta) * theta0 * cos (phi) - mag_NS.get_Rns() * sin(theta) * sin(phi) * phi0;
	k[1] = c * sin(theta) * sin(phi) + mag_NS.get_Rns() * cos(theta) * theta0 * sin (phi) + mag_NS.get_Rns() * sin(theta) * cos(phi) * phi0;
	k[2] = c * cos(theta)            - mag_NS.get_Rns() * sin(theta) * theta0; 


	// two rotations in respect to normal: R(azimuth) * R(acos(mu)) * n
	//k[0] = c * cos (theta) * cos (phi); // along the orthogonal direction to NS surface
	//k[1] = c * cos (theta) * sin (phi); // tbd
	//k[2] = c * sin (theta); 

	//cout << "Be careful, position and k vector might not be right at the moment" << endl;

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

	B_inst[0] = mg->Br     (pos_r[0], pos_r[1]);
	B_inst[1] = mg->Btheta (pos_r[0], pos_r[1]);
	B_inst[2] = mg->Bphi   (pos_r[0], pos_r[1]);

	res = kr[0] * B_inst[0] + kr[1] * B_inst[1] + kr[2] * B_inst[2]; // scalar product 

	return res;
}

// eq. (26) in Fernandez & Thompson (2007) article

double photon::beta_plus () {

	double mu_val, omega_c, res, omegac2;

	mu_val  = get_mu();
	omega_c = mg->omega_B (pos_r[0], pos_r[1]);

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
	omega_c = mg->omega_B (pos_r[0], pos_r[1]);

	omegac2 = pow(omega_c / omega, 2.0);

	if ((mu_val*mu_val + omegac2) < 1) 
		return -20; // error message - scattering is not possible

	res = 1.0 / (omegac2 + mu_val*mu_val) * (mu_val - omega_c/omega * sqrt(omegac2 + mu_val*mu_val - 1.0));

	return res;
}

double photon::propagate_one_step (double delta_t) {

	double k_new[3];
	double pos_new[3], r_new, theta_new, phi_new;
	double k_r_new[3];

	pos_new[0] = pos[0] + k[0] * delta_t;
	pos_new[1] = pos[1] + k[1] * delta_t;
	pos_new[2] = pos[2] + k[2] * delta_t;

	r_new     = sqrt(pos_new[0]* pos_new[0] + pos_new[1]*pos_new[1] + pos_new[2]*pos_new[2]);
	theta_new = atan2 (sqrt(pos_new[0]*pos_new[0] + pos_new[1]*pos_new[1])/r_new, pos_new[2] / r_new );
	phi_new   = atan2 (pos_new[1], pos_new[0]);

	k_r_new[0] = (pos_new[0] * k[0] + pos_new[1] * k[1] + pos_new[2] * k[2]) / r_new;
	k_r_new[1] = (k_r_new[0] * cos(theta_new) - k[2]) / (r_new * sin(theta_new));
	if ( abs(phi_new) > 0.1) 
		k_r_new[2] = (k_r_new[0] * sin(theta_new) * cos(phi_new) + r_new * cos(theta_new) * k_r_new[1] * cos(phi_new) - k[0]) / (r_new * sin(theta_new) * sin(phi_new));
	else
		k_r_new[2] = (k[1] - k_r_new[0] * sin(theta_new) * sin(phi_new) + r_new * cos(theta_new) * k_r_new[1] * sin(phi_new)) / (r_new * sin(theta_new) * cos(phi_new));



	pos[0] = pos_new[0];
	pos[1] = pos_new[1];
	pos[2] = pos_new[2];

	pos_r[0] = r_new;
	pos_r[1] =theta_new;
	pos_r[2] = phi_new;

	kr[0] = k_r_new[0];
	kr[1] = k_r_new[1];
	kr[2] = k_r_new[2];

	return 0;

}



int main () {

	double dr, dtheta, dphi, dt;
	photon * pht1;

	magnetosphere mg (1e14, 0.1, 1e6, 1e6);

	photon pht (M_PI/2.0, 0.0, 1.0e6, 1.0, mg);
	cout << "Polarisation state: " << pht.get_polarisation_state() << endl;
	cout << "Beaming parameter: "<<   pht.get_beaming_parameter () << endl;

	dt = 1e-3;

	pht.print_pos();
	pht.print_k ();

	dr     = - pht.r();
	dtheta = - pht.theta();
	dphi   = - pht.phi ();

	pht.propagate_one_step(dt);
	pht.print_pos();

	dr     += pht.r();
	dtheta += pht.theta();
	dphi   += pht.phi();

	cout << "kr = " << dr / dt << " ktheta = " << dtheta / dt << " kphi = " << dphi / dt << endl;

	pht.print_k();



	pht.propagate_one_step(dt);
	pht.print_pos();
	pht.propagate_one_step(dt);
	pht.print_pos();
	pht.propagate_one_step(dt);
	pht.print_pos();
	pht.propagate_one_step(dt);
	pht.print_pos();

	cout << "Is it inside NS? "<< pht.is_inside_ns() << endl;

	ofstream ofile ("photon_list.txt");

	for (int i = 0; i < 100; i++) {
		pht1 = new photon (M_PI/2.0, 0.0, 1.0e6, 1.0, mg);
		pht1->propagate_one_step(dt);
		ofile << pht1->x() << "\t" << pht1->y() << "\t" << pht1->z() << endl;
		delete pht1;
	}
	ofile.close();


	//cout << "mu: " <<                 pht.get_mu () << endl; 
	//cout << "omega: " <<              pht.get_omega () << endl;






}
