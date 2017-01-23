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
    double xpom=0; // x for F2
    double t=0.1;
    int A=1;
    DGLAPDist *gd=0;  // Initialized and used if we have nucleus consisting of ipsatnucleons (old ipsat)
    
    bool ms=false;
    bool smallb=false;
   
    
    
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
                amp = new Smooth_ws_nuke(A);
            }
            
        }
        else if (string(argv[i])=="-Q2")
            Qsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-xpom")
            xpom = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-ms")
            ms=true;
        else if (string(argv[i])=="-smallb")
            smallb=true;
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
    
    
    cout << "# " << InfoStr() ;
    InitializeWSDistribution(197);
        

    
    InclusiveDiffraction diffraction(amp);
    
    cout << "#Q^2=" << Qsqr << " GeV^2, xpom=" << xpom << endl;
    
    if (ms)
    {
        double ms = diffraction.DiffractiveStructureFunction_qqg_MS_T(xpom, 0, Qsqr);
        cout << "#MS result" << endl;
        cout << ms << endl;
        return 0;
    }
    
    if (smallb)
    {
        cout << "# GBW result at small beta" << endl;
        for (double beta=0.0003; beta<0.1; beta*=2)
        {
            double gbw = diffraction.DiffractiveStructureFunction_qqg_GBW_T(xpom, beta, Qsqr);
            cout << beta << " " << gbw << endl;
            return 0;
        }
    }
    
    for (double beta=0.01; beta<=0.98; beta+=0.04)
    {
		double f_qq_t = diffraction.DiffractiveStructureFunction_qq_T(xpom, beta, 5);
		double f_qq_l = diffraction.DiffractiveStructureFunction_qq_L(xpom, beta, 5);
		double gbw = diffraction.DiffractiveStructureFunction_qqg_GBW_T(xpom, beta, Qsqr);
        cout << beta << " " << f_qq_t << " " << f_qq_l << " " << gbw << endl;
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

