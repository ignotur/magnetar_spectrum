#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "monte_carlo.h"
#include "planck.h"
#include "photon.h"
#include <cassert>

int main () {

        photon * pht1;
        bool success;
        double phi_surf, theta_surf;
        double T_surf;
        int N_phot;
        double E_init;
	double sum_fbeta;
	double B, Bx, By, Bz;
	double dt;
	double dtauv;

        T_surf = 2e6;

        N_phot = 1;

        magnetosphere mg (1e14, 0.1, 348135750.1848, 1.e6);

	B = mg.B (3e6, 1.52);
	Bx = mg.Bx (3e6, 1.52, 0.0);
	By = mg.By (3e6, 1.52, 0.0);
	Bz = mg.Bz (3e6, 1.52, 0.0);

	cout << "* Test of magnetic field transformation from spherical coordinate system to the Carthesian coordinate system" << endl;
	//cout << "B = " << B <<" other B = "<< sqrt(Bx*Bx + By*By + Bz*Bz) << " diff: " << abs(B - sqrt(Bx*Bx + By*By + Bz*Bz)) <<endl;

	assert (abs(B - sqrt(Bx*Bx + By*By + Bz*Bz)) < 1e-6);
	cout << "* Succes" << endl;

        phi_surf   = 2.0 * M_PI * uniform(0.0, 1.0);
        theta_surf = acos(1 - 2 * uniform(0.0, 1.0));

        pht1 = new photon (theta_surf, phi_surf, T_surf, 1.0, &mg);

        cout << "* Test of |k| in the Carthesian coordinate system" << endl;
        //cout << "Check if |k|==1: "<<pow(pht1->kx(), 2.0) + pow(pht1->ky(), 2.0) + pow(pht1->kz(), 2.0) << endl;
	assert (abs(1 - (pow(pht1->kx(), 2.0) + pow(pht1->ky(), 2.0) + pow(pht1->kz(), 2.0))) < 1e-6 );
	assert (abs(1 - (pow(pht1->k_r(), 2.0) + pow(pht1->k_theta(), 2.0) + pow(pht1->k_phi(), 2.0))) < 1e-6 );

	cout << "* Succes" << endl;


	cout << "mu is: " << pht1->get_mu() << endl;

	// Let us test the normalisation of f_beta distribution here
	//
	//
	cout << "* Test of f_beta numerical normalisation" << endl;

	sum_fbeta =  mg.f_beta(-1.0) + mg.f_beta(-0.9) + mg.f_beta(-0.8) + mg.f_beta(-0.7) + mg.f_beta(-0.6);
	sum_fbeta += mg.f_beta(-0.5) + mg.f_beta(-0.4) + mg.f_beta(-0.3) + mg.f_beta(-0.2) + mg.f_beta(-0.1);
	sum_fbeta += mg.f_beta(0.0)  + mg.f_beta(0.1)  + mg.f_beta(0.2) + mg.f_beta(0.3) + mg.f_beta(0.4);
	sum_fbeta += mg.f_beta(0.5)  + mg.f_beta(0.6)  + mg.f_beta(0.7) + mg.f_beta(0.8) + mg.f_beta(0.9);
	sum_fbeta += mg.f_beta(1.0);

	sum_fbeta = sum_fbeta * 0.1;
	//cout << " sum_fbeta = " << sum_fbeta << endl;
	
	assert (abs(sum_fbeta - 1) < 0.1);

	cout << "* Succes" << endl;
	


        ofstream ofile ("propagation_test.txt");

	dt = 1e-5;

	for (int i=0; i < 200; i++) {

		dtauv = pht1->propagate_one_step (dt);
		ofile << pht1->x() << "\t" << pht1->y() << "\t" << pht1->z() <<"\t" << pht1->r() << "\t" <<  pht1->theta() << "\t" <<  pht1->phi() <<"\t" <<dtauv <<"\t" << pht1->get_mu() <<"\t"<< mg.omega_B (pht1->r(), pht1->theta()) / pht1->get_omega() << endl;

		//cout << "Step: "<<i << endl;

		//cout << (pow(pht1->k_r(), 2.0) + pow(pht1->k_theta(), 2.0) + pow(pht1->k_phi(), 2.0)) << endl;
		//cout << pow(pht1->kx(), 2.0) + pow(pht1->ky(), 2.0) + pow(pht1->kz(), 2.0) << endl;

	        assert (abs(1 - (pow(pht1->kx(), 2.0) + pow(pht1->ky(), 2.0) + pow(pht1->kz(), 2.0))) < 1e-6 );
       		assert (abs(1 - (pow(pht1->k_r(), 2.0) + pow(pht1->k_theta(), 2.0) + pow(pht1->k_phi(), 2.0))) < 1e-6 );

		assert (abs(pht1->get_mu()) <= 1.0);

		//cout << "mu = " << pht1->get_mu() << "\t "<< (pht1->kx()*mg.Bx(pht1->r(), pht1->theta(), pht1->phi()) + pht1->ky()*mg.By(pht1->r(), pht1->theta(), pht1->phi()) + pht1->kz()*mg.Bz(pht1->r(), pht1->theta(), pht1->phi())) / mg.B(pht1->r(), pht1->theta())   << endl;

		assert (abs(pht1->get_mu() - (pht1->kx()*mg.Bx(pht1->r(), pht1->theta(), pht1->phi()) + pht1->ky()*mg.By(pht1->r(), pht1->theta(), pht1->phi()) + pht1->kz()*mg.Bz(pht1->r(), pht1->theta(), pht1->phi()) ) / mg.B(pht1->r(), pht1->theta())  ) < 1e-6 );




	}





        ofile.close();

}

