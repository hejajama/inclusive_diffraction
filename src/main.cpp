/*
 * Inclusive diffraction
 * Main reference: 0805.4071
 * 
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016
 * 
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include <gsl/gsl_rng.h>

#include <tools/tools.hpp>
#include <gsl/gsl_errno.h>

#include <dipole.hpp>
#include <smooth_ws_nuke.hpp>
#include <ipsat_nucleons.hpp>
#include <ipsat_proton.hpp>
#include <vector.hpp>
#include <ipglasma.hpp>
#include <nucleons.hpp>

#include "inclusive_diffraction.hpp"

using namespace std;

gsl_rng *global_rng;
string InfoStr();
DipoleAmplitude* amp;

using namespace std;

int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler); // Do not let code to crash if an error occurs, we should
    // have error handling everywhere!
    
    double Qsqr=0;
    double xbj=0; // x for F2
    double t=0.1;
    int A=1;
    DGLAPDist *gd=0;  // Initialized and used if we have nucleus consisting of ipsatnucleons (old ipsat)
    
   
    
    
    cout << "# Inclusive Diffraction by H. Mäntysaari <mantysaari@bnl.gov>, 2016" << endl;
    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    
    if (string(argv[1])=="-help")
    {
        
        return 0;
    }
    
    for (int i=1; i<argc; i++)
    {
       
        if (string(argv[i])=="-dipole")
        {
            A = StrToInt(argv[i+1]);
            if (A==1)
            {
                if (string(argv[i+2])=="ipsatproton")
                {
                    amp = new Ipsat_Proton;
                    ((Ipsat_Proton*)amp)->SetProtonWidth(StrToReal(argv[i+3]));
                    ((Ipsat_Proton*)amp)->SetQuarkWidth(StrToReal(argv[i+4]));
                  
                        ((Ipsat_Proton*)amp)->SetShape(GAUSSIAN);
                        if (argc > i+5)
                        {
                            if (string(argv[i+5])=="fluxtube")
                            {
                                ((Ipsat_Proton*)amp)->SetStructure(CENTER_TUBES);
                                ((Ipsat_Proton*)amp)->SetFluxTubeNormalization(StrToReal(argv[i+5]));
                            }
                            else if (string(argv[i+5]).substr(0,1)!="-")
                            {
                                cerr << "Unknown ipsatproton option " << argv[i+4] << endl;
                                exit(1);
                            }
                        }
                    
                }
                else if (string(argv[i+2])=="ipglasma")
                    amp = new IPGlasma(argv[i+3]);
                else
                {
                    cerr << "Unknown dipole " << argv[i+1] << endl;
                    return -1;
                }
            }
            else
            {
                // Construct nucleus
                std::vector<DipoleAmplitude* > nucleons;
                for (int j=0; j<A; j++)
                {
                    if (string(argv[i+2])=="ipsatproton")
                    {
                        if (j==0)
                            gd = new DGLAPDist;
                        //Ipsat_Proton *nucleon = new Ipsat_Proton(gd);
                        Ipsat_Proton *nucleon = new Ipsat_Proton();
                        nucleon->SetProtonWidth(StrToReal(argv[i+3]));
                        nucleon->SetQuarkWidth(StrToReal(argv[i+4]));
                       
                        
                            nucleon->SetShape(GAUSSIAN);
                            if (argc > i+5)
                            {
                                if (string(argv[i+5])=="fluxtube")
                                {
                                    nucleon->SetStructure(CENTER_TUBES);
                                    nucleon->SetFluxTubeNormalization(StrToReal(argv[i+5]));
                                }
                                else if (string(argv[i+5]).substr(0,1)!="-")
                                {
                                    cerr << "Unknown ipsatproton option " << argv[i+5] << endl;
                                    exit(1);
                                }
                            }
                        
                        nucleons.push_back(nucleon);

                    }
                }
                amp = new Nucleons(nucleons);
            }
            
        }
       
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unknown parameter " << argv[i] << endl;
            exit(1);
        }
        
    }

    
    
    // Initialize global random number generator
    gsl_rng_env_setup();
    global_rng = gsl_rng_alloc(gsl_rng_default);

    amp->InitializeTarget();
    
    
    cout << "# " << InfoStr() << endl;
    InitializeWSDistribution(197);
        

    
    InclusiveDiffraction diffraction(amp);
    
    
    for (double beta=0.05; beta<0.95; beta+=0.05)
    {
		//double f_qq_t = diffraction.DiffractiveStructureFunction_qq_T(0.001, beta, 5);
		//double f_qq_l = diffraction.DiffractiveStructureFunction_qq_L(0.001, beta, 5);
		//cout << beta << " " << 0.001*f_qq_t << " " << 0.001*f_qq_l << endl;
		double gbw = diffraction.DiffractiveStructureFunction_qqg_GBW_T(0.001, beta, 5);
		//double ms = diffraction.DiffractiveStructureFunction_qqg_MS_T(0.001, beta, 5);
		cout << beta << " " << gbw << endl;
		//cout << beta << " " << f_qq_t << " " << f_qq_l << endl;
	}
    
    
    
    
    gsl_rng_free(global_rng);


    
    delete amp;
    
    if (gd != 0)
        delete gd;
    
}


string InfoStr()
{
    stringstream info;
    
    
    info << amp->InfoStr();

    
    
    return info.str();

}

