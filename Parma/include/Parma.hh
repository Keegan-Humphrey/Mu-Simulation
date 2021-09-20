#ifndef PARMA4_HH
#define PARMA4_HH

//#include <G4ThreeVector.hh>
//#include <G4Event.hh>

namespace Parma4 {

const int nebin = 1000;
const int nabin = 1000;
const int defaultSeed = 5918;

class Parma {
    public:
        //Constructor
        Parma();
        Parma(int particleIP, int genSeed);

        //Destructor
        ~Parma();

        //TODO: read from string
        // bool readString(std::string config);

        //initialize Parma instance
        bool init();
        bool init(int particleIP, int genSeed);

        //gen next event
        double * getParticle();

        int getIP() { return ip; }
        int getSeed() { return seed; }

        void setIP(int n) { ip = n; }

    private:
        bool    isInit;

        int ip;
        int seed;
        int ivars[2]; //ie, ia
        double  mass;
        double  ehigh[nebin+1], etable[nebin+1], ahigh[nabin+1];
        double  atable[nebin+1][nabin+1], flux[nebin][nabin+1];
        // double  ehigh[nebin+1], etable[nebin+1], ahigh[nabin+1],
        //         atable[nebin+1][nabin+1], flux[nebin][nabin+1];
        

};

} /* namespace Parma4 */ //////////////////////////////////////////////////////////////////////
#endif