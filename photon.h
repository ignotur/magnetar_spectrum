
using namespace std;

class magnetosphere {
        private:
                int idk;
                double p1;         // Parameter describing the magnetic field twist
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
                double get_p   () {return p1;};
                double Br      (double r, double theta) {return -0.5*Bpole * pow( Rns/r , 2.0+p1) * dFfun (theta);};  // Eqs. (2)  Fernandez & Thompson (2007)
                double Btheta  (double r, double theta) {return  0.5*Bpole * pow( Rns/r , 2.0+p1) * p1 * Ffun (theta) / sin(theta);};
                double Bphi    (double r, double theta) {return Btheta (r, theta) * sqrt(C / (p1*(p1+1))) * pow(Ffun(theta), 1.0/p1);};
                double B       (double r, double theta) {return sqrt(pow(Br(r, theta), 2.0) + pow(Btheta (r,theta), 2.0) + pow(Bphi(r,theta), 2.0));};
                double Bx      (double r, double theta, double phi) {return Br (r, theta) * sin(theta) * cos(phi) + Btheta (r,theta) * cos(theta)*cos(phi) - Bphi (r,theta) * sin(phi); };
                double By      (double r, double theta, double phi) {return Br (r, theta) * sin(theta) * sin(phi) + Btheta (r,theta) * cos(theta)*sin(phi) + Bphi (r,theta) * cos(phi); };
                double Bz      (double r, double theta, double phi) {return Br (r, theta) * cos(theta)  - Btheta (r,theta) * sin(theta); };
                double ne      (double r, double theta) {return (p1+1) / (4.0 * M_PI * echarge) * (Bphi(r,theta) / Btheta (r,theta)) * B (r, theta) / r / beta_bulk;};
                double omega_B (double r, double theta) {return echarge * B(r,theta) / (me * speed_of_light);};
                double get_Rns () {return Rns;};
                double f_beta  (double beta_v); // Velocity distribution of charged particles in magnetosphere
                void   normalise_f_beta ();     // Compute and store numerical coeffitient for f_beta distribution
};

class photon {
        private:
                double b;         // Beaming parameter
                double omega_em;   // Emission frequency
                double omega;     // Instantenious photon frequency
                int s;            // Instantenious Photon polarisation state: 0 - photon does not exist; 1 - ordinary mode; 2 - extraordinary mode
                double mu;        // Cosine of the angle between the initial photon direction and magnetic field
                double azimuth;   // Azimuth angle for photon
                double pos_em[3];  // Emission position
                double pos [3];   // Instantenious location of photon in the Carthesian coordinate system [x,y,z]
                double pos_r [3]; // Instantenious location of photon in the spherical coordinate system [r,theta,phi]
                double k [3];     // Instantenious k vector for the photon in the Carthesian coordinate system [kx, ky, kz] - this is a unit vector in direction of the photon propagation
                double kr[3];     // Instantenious k vector for the photon in the spherical coordinate system [kr, ktheta, kphi]
                double c;         // Speed of light
                double re;        // Classical electron radius
                int max_number_of_propagation_steps; // maximum number of integration steps before the photon get discarded
                magnetosphere *mg;// Object of the twisted magnetosphere class. Physically the magnetosphere through which the photon propagates

        public:
                photon (double theta, double phi, double T, double beaming, magnetosphere); // Constructor which corresponds to emission of the photon from the NS surface path with coordinates phi, theta and temperature T
                int    propagate (bool); // propagate all the way from NS atmosphere to the last scattering, flag is to print photon history for debugging
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
                double kx                     () {return k[0];};
                double ky                     () {return k[1];};
                double kz                     () {return k[2];};
                void   print_pos              () {cout << "x = "<<pos[0]<<" y = "<<pos[1] <<" z = "<<pos[2] << "; r = "<<pos_r[0]<<" theta = "<<pos_r[1] << " phi = "<<pos_r[2] <<endl; };
                void   print_k                () {cout << "kx = "<<k[0]<<" ky = "<<k[1] <<" kz = "<<k[2]    << "; kr = " <<kr[0] << " k_theta = "<<kr[1] << " k_phi = "<<kr[2] << endl;}
                bool   is_inside_ns           () {if (pos_r[0] < mg->get_Rns()) return true; else return false;};
                double dist_from_ns_center    () {return pos_r[0];};
                double Sf                     (int, int);
                void   update_kr              ();
                double get_E_keV              () {return 6.582119569e-16*omega / 1.0e3;};
};

