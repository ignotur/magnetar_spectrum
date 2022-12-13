#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "monte_carlo.h"
#include "planck.h"
#include "photon.h"

using namespace std;

// Code to model spectra of the radiation emitted by a hot NS and scattered in the twisted magnetosphere
// The main scattering mechanism is the resonant cyclotron scattering
// The code is based on article by Nobili, Turolla & Zane (2008)
// Written by Dr. Andrei P. Igoshev ignotur@gmail.com


photon::photon (double theta, double phi, double T, double beaming, magnetosphere *mag_NS) {

	double theta0, phi0;

	static bool called_ever;

	if (!called_ever) {
		called_ever = true;
		srand (time(NULL));
	}


	b = beaming;

	s = rand() % 2 + 1; // Draw the initial polarisation from the uniform distribution

	mu = pow (uniform (0.0, 1.0), 1.0 / b); // cos is uniform taking into account the beaming factor // assumed to be an angle between normal and photon direction, could instead be an angle between local field direction and photon emission - still need to check

	azimuth = 2.0 * M_PI * uniform (0.0, 1.0);

	omega = sed_planck (T);

	c  = 2.99792458e10; // cm/s
	re = 2.8179403227e-13; // cm 

	max_number_of_propagation_steps = 5000;

	mg = mag_NS;

	pos_em[0] = mag_NS->get_Rns() * sin (theta) * cos (phi); // Theta is computed from the pole down
	pos_em[1] = mag_NS->get_Rns() * sin (theta) * sin (phi);
	pos_em[2] = mag_NS->get_Rns() * cos (theta);

	pos[0] = pos_em[0];
	pos[1] = pos_em[1];
	pos[2] = pos_em[2];


	pos_r[0] = mag_NS->get_Rns();
	pos_r[1] = theta;
	pos_r[2] = phi;

	//kr[0] = c;
	//kr[1] = acos(mu);
	//kr[2] = azimuth;

	theta0 = acos(mu);
	phi0   = azimuth;

	//k[0] = c * sin(theta) * cos(phi) + mag_NS.get_Rns() * cos(theta) * theta0 * cos (phi) - mag_NS.get_Rns() * sin(theta) * sin(phi) * phi0;
	//k[1] = c * sin(theta) * sin(phi) + mag_NS.get_Rns() * cos(theta) * theta0 * sin (phi) + mag_NS.get_Rns() * sin(theta) * cos(phi) * phi0;
	//k[2] = c * cos(theta)            - mag_NS.get_Rns() * sin(theta) * theta0; 
	//
	
	
	k[0] = sin(phi)*sin(phi0)*sin(theta0)*cos(theta) + sin(theta) *cos(theta0) + sin(theta0) *cos(phi) * cos(phi0) *cos(theta);
	k[1] = -sin(phi) *sin(theta0) *cos(phi0) + sin(phi0) *sin(theta0) * cos(phi);
	k[2] = -sin(phi) *sin(phi0) *sin(theta) * sin(theta0) - sin(theta) *sin(theta0)*cos(phi)*cos(phi0) +cos(theta) *cos(theta0);

	//cout << "Be careful, position and k vector might not be right at the moment" << endl;
	//
	//cout << "At this moment p = "<<mg->get_p()<<endl;

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
	kB        = 1.3807e-16;       // cm^2 g s^-2 K^-1

	ifstream infile ("F_magnetosphere.txt");

	infile >> N;
	infile >> p1; 
	infile >> C;

	mu = new double [N];
	F  = new double [N];

	for (int i =0; i < N; i++) {

		infile >> mu[i] >> F[i];

	}

	norm_f = 1.0;

	normalise_f_beta ();

	cout << "We initiliased twisted magnetosphere using file F_magnetosphere.txt" << endl;
	cout << "N = " << N << " \t p = " << (double) p1 << " \t C = "<< C << endl;

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
	Theta_e = kB * Te / (me * speed_of_light * speed_of_light);
	gamma_ap = gamma_v * gamma_bulk * (1.0 - beta_v * beta_bulk);

	//cout << "Theta_e = "<< Theta_e << "\t" << "gamma_ap = " << gamma_ap << endl;

	res = norm_f * exp (-gamma_ap / Theta_e); // eq. (4) in Nobili, Turolla & Zane, normalised numerically 
	
	return res;

}

void   magnetosphere::normalise_f_beta () {

	double k1, k2, k3, k4, res, dbeta;
	int N;

	N = 1000.0;

	dbeta = 2.0 / N;

	res = 0.0;

	for (int i = 0; i < N-1; i++) {

		k1 = f_beta (-1.0 + dbeta * i);
		k2 = f_beta (-1.0 + dbeta * (i + 0.5));
		k4 = f_beta (-1.0 + dbeta * (i + 1.0));
		res = res + (k1 + 4.0 * k2 + k4) / 6.0 * dbeta;
		//cout << i << "\t"<<-1.0 + dbeta * i << "\t" << k1 << endl;

	}

	norm_f = 1.0 / res;

	cout << "Our new normalisation factor is: "<<norm_f << endl;

	//return res;

}

// Compute the ange in respect to local magnetic field
double photon::get_mu () {

	double B_inst [3];
	double res;
	double module_kr;
	double module_B;

	B_inst[0] = mg->Br     (pos_r[0], pos_r[1]);
	B_inst[1] = mg->Btheta (pos_r[0], pos_r[1]);
	B_inst[2] = mg->Bphi   (pos_r[0], pos_r[1]);

	res = kr[0] * B_inst[0] + kr[1] * B_inst[1] + kr[2] * B_inst[2]; // scalar product 

	module_kr = sqrt(kr[0]*kr[0] + kr[1]*kr[1] + kr[2]*kr[2]);
	module_B  = sqrt(B_inst[0]*B_inst[0] + B_inst[1]*B_inst[1] + B_inst[2]*B_inst[2]);

	res = res / module_kr / module_B;

	cout << "Intermediate steps of mu calculations: vec B = "<<B_inst[0]<<" "<<B_inst[1]<<" "<<B_inst[2]<<"\t vec kr = "<<kr[0]<<" "<<kr[1]<<" "<<kr[2]<<endl;
	cout << "Photon position: "<<pos_r[0] << "\t" << pos_r[1] << "\t" << pos_r[2] << endl;
	cout << "Magnetic field (Carthesian) "<< mg->Bx(pos_r[0], pos_r[1], 0) << endl;
	cout << "Check if class is initialised correctly: "<< mg->get_Rns ()<<endl;
  
	cout << "Inside get_mu"<<endl;
	cout << "Standard calculations: "<<res<<endl;
	cout << "Carthesian calculations: "<<(k[0]*mg->Bx(pos_r[0], pos_r[1], pos_r[2]) + k[1]*mg->By(pos_r[0], pos_r[1], pos_r[2]) + k[2]*mg->Bz(pos_r[0], pos_r[1], pos_r[2])) / mg->B(pos_r[0], pos_r[1]) / sqrt(k[0]*k[0] + k[1]*k[1]+k[2]*k[2])<<endl;
	cout << "Module B = "<<mg->B(pos_r[0], pos_r[1]) <<"\t alternatively "<<module_B<<endl;
	cout << "Module k = "<<sqrt(k[0]*k[0] + k[1]*k[1]+k[2]*k[2]) << "\t alternatively "<<module_kr<<endl;

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

// part of eq. (16 and 17). It corresponds to individual summands: first index - initial polarisation, second index - plus or minus
double photon::Sf (int ind1, int ind2) {

	double mu_v;
	double beta_k;
	double res;

	if ((ind1 < 0) || (ind1 > 2)) {
		cout << "Index i in function Sf is wrong "<<ind1 << endl;
		exit(2);
	}
	if ((ind2 < 0) || (ind2 > 2)) {
		cout << "Index j in function Sf is wrong "<<ind2 << endl;
		exit(2);
	}

	mu_v = get_mu();


	if (ind1 == 1) {
		if (ind2 == 1) {
			beta_k = beta_plus();
			res = abs(mu_v - beta_k) / (1.0 - mu_v * beta_k) * mg->f_beta  (beta_k);
		}
		if (ind2 == 2) {
			beta_k = beta_minus();
			res = abs(mu_v - beta_k) / (1.0 - mu_v * beta_k) * mg->f_beta  (beta_k);
		}
	}
	if (ind1 == 2) {
		if (ind2 == 1){
			beta_k = beta_plus();
			res = (1.0 - mu_v * beta_k) / abs(mu_v - beta_k) * mg->f_beta  (beta_k);
		}
		if (ind2 == 2){
			beta_k = beta_minus();
			res = (1.0 - mu_v * beta_k) / abs(mu_v - beta_k) * mg->f_beta  (beta_k);
		}
	}
	

return res;
}

double photon::propagate_one_step (double delta_t) {

	double k_new[3];
	double pos_new[3], r_new, theta_new, phi_new;
	double k_r_new[3];
	double beta_plus_v, beta_minus_v;
	double omega_B;
	double dtau;
	double mu_v, c;

	c  = 2.99792458e10; // cm/s

	pos_new[0] = pos[0] + k[0] * delta_t * c;
	pos_new[1] = pos[1] + k[1] * delta_t * c;
	pos_new[2] = pos[2] + k[2] * delta_t * c;

	cout << "New pos: "<< pos_new[0] / 1e6 <<"\t" << pos_new[1] / 1e6 << "\t" << pos_new[2] / 1e6 << endl;  

	r_new     = sqrt(pos_new[0]* pos_new[0] + pos_new[1]*pos_new[1] + pos_new[2]*pos_new[2]);
	theta_new = atan2 (sqrt(pos_new[0]*pos_new[0] + pos_new[1]*pos_new[1]), pos_new[2] );
	phi_new   = atan2 (pos_new[1], pos_new[0]);

	cout << "New r: "<< r_new  <<endl;

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

	beta_plus_v  = beta_plus();
	beta_minus_v = beta_minus();
	cout << "Test (1)" << endl;


	omega_B = mg->omega_B (pos_r[0], pos_r[1]);
	mu_v = get_mu();


	if ((beta_plus_v != -20) && (beta_minus_v != -20)) {

		if (s == 1) {
//			dtau = 2.0*M_PI*M_PI*re*c*mg->ne (pos_r[0], pos_r[1]) *omega_B / (omega*omega) * ( abs(mu_v - beta_plus_v) / (1.0 - mu_v*beta_plus_v) * mg->f_beta  (beta_plus_v) + abs(mu_v - beta_minus_v) / (1.0 - mu_v*beta_minus_v) * mg->f_beta  (beta_minus_v) * c * delta_t  );

			dtau = c*delta_t * 2.0*M_PI*M_PI*re*c*mg->ne (pos_r[0], pos_r[1]) *omega_B / (omega*omega) * ( Sf (1,1) + Sf (1,2) );

		}
		else if (s == 2) {
			dtau = c*delta_t * 2.0*M_PI*M_PI*re*c*mg->ne (pos_r[0], pos_r[1]) *omega_B / (omega*omega) * ( Sf (2,1) + Sf (2,2) );
		}

	}
	else 
		dtau = 0.0;
  	
        //cout <<"pos = "<<pos_r[0]<<"\t"<< pos_r[1] <<" in respect to NS radius = "<<pos_r[0] / mg->get_Rns() <<endl;
	//cout <<"Bphi = "<<mg->Bphi (pos_r[0], pos_r[1]) << "\t Btheta = "<< mg->Btheta(pos_r[0], pos_r[1]) <<endl;
	//cout <<"What is p? = "<<mg->get_p() <<endl;
        //cout <<"Magnetospheric properties: "<<mg->ne (pos_r[0], pos_r[1]) << endl;
	//cout <<"Comparison of optical depths: "<< 2.0*M_PI*M_PI*re*c*mg->ne (pos_r[0], pos_r[1]) *omega_B / (omega*omega) * c * delta_t << "\t" << dtau << endl;
	//cout <<"Remaining factors are: f(beta_k) = "<< mg->f_beta  (beta_minus_v) <<"\t" << mg->f_beta  (beta_plus_v) <<endl;
	//cout <<"mu_v - beta_plus_v = "<<abs(mu_v - beta_plus_v) << endl;
	//cout <<"mu_v = "<<mu_v << endl;
	//cout <<"Distance from NS = "<<pos_r[0]/1e6 <<" mu_v = "<<mu_v<<"\t k0*delta_t = "<<k[0]*delta_t/1e6<<endl;

	return dtau;

}

void photon::update_kr() {
	double theta1, phi1;

	theta1 = pos_r[1];
	phi1   = pos_r[2];

	kr[0] =  k[0] * sin(theta1) * cos(phi1) + k[1] * sin(theta1) * sin(phi1)  + k[2] * cos(theta1);
	kr[1] =  k[0] * cos(theta1) * cos(phi1) + k[1] * cos(theta1) * sin(phi1)  - k[2] * sin(theta1);
	kr[2] = -k[0] * sin(phi1)               + k[1] * cos(phi1);

	//cout << "mod k = "<<sqrt(k[0]*k[0] + k[1]*k[1]+k[2]*k[2]) << " mod k2 = "<<sqrt(kr[0]*kr[0]+kr[1]*kr[1]+kr[2]*kr[2])<<endl;

	//exit(2);


}

// Scattering event - photon gets new polarisation, frequency and k-vector
int photon::scatter () {

	double u1, u2, u3, u4;
	double phi_new, theta_new, s_new, new_mu;
	double beta_plus_v, beta_minus_v;
	double beta_val;
	double Bx, By, Bz;
	double mu_k, mu_b, mu_val;
	double k_new [3];
	double alpha11;
	double omega_new;

	u1 = uniform (0.0, 1.0);
	u2 = uniform (0.0, 1.0);
	u3 = uniform (0.0, 1.0);
	u4 = uniform (0.0, 1.0);
 
	if ((s == 1) && (u1 > 0.25))
		s_new = 2;
	if ((s == 2) && (u1 > 0.75))
		s_new = 1;

	if (((s == 1) && (s_new == 1)) || ((s==2) && (s_new==1))) {
		theta_new = 2.0 * u4 - 1.0;
		theta_new = pow(theta_new, 1.0 / 3.0);
		//theta_new = acos(theta_new);
	}
	if (((s == 1) && (s_new == 2)) || ((s==2) && (s_new==2))) {
		theta_new = 2.0 * u4 - 1.0;
		//theta_new = acos(theta_new);
	}

	phi_new = 2.0 * M_PI * u3;

	mu_val = get_mu ();
	beta_plus_v  = beta_plus();
	beta_minus_v = beta_minus();

	if ((s == 1) && (u2 < Sf(1,1) / (Sf(1,1) + Sf(1,2)) ))
		beta_val = beta_plus_v;
	else if (s == 1) 
		beta_val = beta_minus_v;
	if ((s == 2) && (u2 < Sf(2,1) / (Sf(2,1) + Sf(2,2)) ))
		beta_val = beta_plus_v;
	else if (s == 2)
		beta_val = beta_minus_v;

	new_mu = (theta_new + beta_val) / (1.0 + beta_val * theta_new);

	Bx = mg->Bx (pos_r[0], pos_r[1], pos_r[2]) / mg->B (pos_r[0], pos_r[1]); // Direction of the magnetic field
	By = mg->By (pos_r[0], pos_r[1], pos_r[2]) / mg->B (pos_r[0], pos_r[1]);
	Bz = mg->Bz (pos_r[0], pos_r[1], pos_r[2]) / mg->B (pos_r[0], pos_r[1]);
	
	if (abs(sqrt(Bx*Bx+By*By + Bz*Bz)-1.0) > 1e-3) {
		cout << "Something strange has happened with magnetic field transformation" <<endl;
		exit(2);
	}


	mu_b = Bz;

	mu_k = mu_b*new_mu + sqrt((1.0 - mu_b*mu_b)*(1.0-new_mu*new_mu)) * cos(phi_new);

	k_new[2] = mu_k;

	alpha11 = sqrt((1.0-mu_b*mu_b)*(1.0 - new_mu*new_mu)) * sin(phi_new);

	k_new[1] = (new_mu - k_new[2] * Bz - alpha11 * Bx/By ) / (By + Bx*Bx / By);

	k_new[0] = (alpha11 + k_new[1] * Bx) / By;

	//cout << "Check for new k_new; |k_new| = "<<sqrt(k_new[0]*k_new[0] + k_new[1]*k_new[1]+ k_new[2]*k_new[2]) <<endl;
	//
	if (abs(sqrt(k_new[0]*k_new[0] + k_new[1]*k_new[1]+ k_new[2]*k_new[2]) - 1.0) > 1e-3) {
		cout << "Something is strange with new k_new at the scattering" <<endl;
		exit(3);
	}

	k_new[0] *= c;
	k_new[1] *= c;
	k_new[2] *= c;


	//cout << "Before the scattering: k = "<<k[0]<<" "<<k[1]<<" "<<k[2]<<" k_new = "<<k_new[0]<<" "<<k_new[1] <<" "<<k_new[2]<<endl;


	omega_new = pow(1.0 / (1.0 - beta_val*beta_val), 2.0) * omega * (1.0 - beta_val * mu_val) * (1.0 + beta_val * new_mu);

	//cout << "Omega before scattering "<< omega << " after "<<omega_new <<endl;

	//cout << "Scattering" << endl;

	k[0] = k_new[0];
	k[1] = k_new[1];
	k[2] = k_new[2];

	update_kr();

	omega = omega_new;

	s = s_new;

	//cout << "Check if our update of mu makes sence - mu_new = "<<new_mu << " compute mu = "<<get_mu()<<endl;
	//cout << "Directly compute mu = "<<(Bx*k_new[0] + By*k_new[1] + Bz*k_new[2]) / c <<endl;
	//exit(0);

	return 0;
}


// In this function we integrate the photon motion
int photon::propagate (bool verbose) {

	double u, tau, dt;
	int cnt;
	bool success;

	u = uniform (0.0, 1.0);

	cnt = 0;
	tau = 0.0;
	dt = 1e-5; // i.e. distance of ~3 km 
	success = true;


	 
	ofstream ofile ("photon_history.txt");

	while (true) { 

		tau += propagate_one_step (dt);
		cnt += 1;

		if (cnt > max_number_of_propagation_steps) {
			success = false;
			break;
		}

		if (tau > - log(u)) {  // Scattering occured, it means that we integrate the optical depth again
			scatter();
			tau = 0.0;

		}


		if (is_inside_ns ()) { // Photon was scattered back into NS, we disregard it
			success = false;
			break;
		}

		if (dist_from_ns_center () > 1000 * mg->get_Rns()) { // Photon left the magnetosphere and has to be recorded
			success = true;
			break;
		}

		//cout << cnt <<"\t" << tau <<endl;
		if (verbose)
			ofile << cnt << "\t" <<x() / 1e6 << "\t" << y() / 1e6 << "\t" << z() / 1e6 << endl; // print current photon location in units of Rns

	}

	return success;
}


