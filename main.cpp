#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "monte_carlo.h"
#include "planck.h"
#include "photon.h"


int main () {

        photon * pht1;
        bool success;
        double phi_surf, theta_surf;
        double T_surf;
        int N_phot;
        double E_init;

        T_surf = 2e6;

        N_phot = 1;

        magnetosphere mg (1e14, 0.1, 348135750.1848, 1.e6);

        ofstream ofile ("photon_list.txt");
        //ofstream oinit ('')
        //

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

