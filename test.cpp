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
	double B, Bx, By, Bz;
	double dt;

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
	cout << "* Succes" << endl;


	pht1->get_mu();


        ofstream ofile ("propagation_test.txt");

	dt = 1e-5;

	for (int i=0; i < 1000; i++) {

		ofile << pht1->x() << "\t" << pht1->y() << "\t" << pht1->z() <<"\t" << pht1->r() << "\t" <<  pht1->theta() << "\t" <<  pht1->phi() << endl;
		pht1->propagate_one_step (dt);


	}





        ofile.close();

}

