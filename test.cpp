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




        ofstream ofile ("photon_list.txt");

        for (int i = 0; i < N_phot; i++) {

                phi_surf   = 2.0 * M_PI * uniform(0.0, 1.0);
                theta_surf = acos(1 - 2 * uniform(0.0, 1.0));

                pht1 = new photon (theta_surf, phi_surf, T_surf, 1.0, mg);

                E_init = pht1->get_E_keV();

                cout << "Test" << endl;
                cout << "Check if |k|==1: "<<pow(pht1->kx(), 2.0) + pow(pht1->ky(), 2.0) + pow(pht1->kz(), 2.0) << endl;

                success = pht1->propagate(false);


                if ((i%100) == 0)
                        cout <<i<<endl;

                if (success)
                        ofile << pht1->theta() <<"\t"<<pht1->phi() <<"\t"<<E_init<<"\t"<<pht1->get_E_keV()<<endl;

                delete pht1;

        }
        ofile.close();

}

