/*
 * Parma.cc
 * Created 2021 by Sameer Erramilli
 * to interface the Parma4 CosmicSpectra program (FORTRAN)
 * with the MATHUSLA simulation (C++).
 */

//#include <G4ThreeVector.hh>
//#include <G4Event.hh>

#include <Parma.hh>
#include <iostream>
using namespace std;

extern"C"{
    /* initialization should give ip, ie, ia, nebin, nabin; 
     * ehigh (ARRAY), etable (ARRAY), 
     * ahigh (ARRAY), atable (ARRAY), e, mass, flux (ARRAY).
     * so input is 7 doubles in an array, plus FIVE more arrays of doubles
     * with particular dimensions.
     * NOTE THAT FORTRAN TREATS ARRAYS DIFFERENTLY FROM C!!!!!! 
     * In Fortran, arrays are stored in COLUMN-MAJOR ORDER;
     * In C, arrays are stored in ROW-MAJOR ORDER.
     */
    // Order: [nebin, nabin, ip, seed], [ie, ia], mass, ehigh[nebin+1], etable[nebin+1],
    //          ahigh[nabin+1], atable[nebin+1,nabin+1], flux[nebin,nabin + 1].
    //
    // The above dimensions are given AS THEY WOULD BE IN C, so they are REVERSED in Fortran.
    // e.g. flux in fortran has dimension (0:nabin, nebin).
    void initializegen_(int [], int [], double *, double [], double [],
                        double [], double [][Parma4::nabin+1], double [][Parma4::nabin+1]);

    // Order: [nebin, nabin, ip, seed], [ie, ia], mass, ehigh[nebin+1], etable[nebin+1],
    //          ahigh[nabin+1], atable[nebin+1,nabin+1], flux[nebin,nabin + 1],
    //          [e+mass, zenith, azimuthal]
    // particle info takes 3 doubles (total energy, zenith, azimuthal)
    // note that of all the inputs taken from initialization, only flux gets modified in getparticle.
    void getparticle_(int [], int [], double *, double [], double [],
                        double [], double [][Parma4::nabin+1], double [][Parma4::nabin+1],
                        double []);
}

namespace Parma4 {
//Taking a break. TODO: I want to have nebin, nabin specific to the instance of the class, but it's complaining about the array sizing then.
//I also need to flesh out what I want to be global to a class and what I want to be specific to a given function.
//I don't want any of nebin, nabin, mass, etc. to be set until init();, since then I can later make a way to specify those quantities externally.

Parma::Parma() {
    isInit = false;
}

Parma::Parma(int particleIP, int genSeed) {
    ip = particleIP;
    seed = genSeed;
    isInit = false;
}

Parma::~Parma() {};

bool Parma::init() { //Default initializer.

    isInit = false;

    if(ip == 0) { //if ip has not been specified yet
        ip = 29;
    }
    if(seed == 0) { //if seed has not been specified.
        cout << "Using default seed as specified in Parma.hh.\n";
        cout << "To specify a seed,  use the Parma(particleIP, genSeed) constructor.\n";
        seed = defaultSeed;
    }

    int iparams[] = { nebin, nabin, ip, seed };

    initializegen_(iparams, ivars, &mass, ehigh, etable, ahigh, atable, flux);

    isInit = true;
    return true;
}

bool Parma::init(int particleIP, int genSeed) { //Initialize w/ custom particle and seed.

    isInit = false;

    ip = particleIP;
    seed = genSeed;

    int iparams[] = { nebin, nabin, ip, seed };

    initializegen_(iparams, ivars, &mass, ehigh, etable, ahigh, atable, flux);

    isInit = true;
    return true;
}

double * Parma::getParticle() { //Returns {ip, total energy, zenith, azimuthal}
    double *particle = new double[4];

    int iparams[] = { nebin, nabin, ip, seed };

    getparticle_(iparams, ivars, &mass, ehigh, etable, ahigh, atable, flux, particle);

    return particle;
}

}