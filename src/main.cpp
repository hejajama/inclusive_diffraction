/*
 * Inclusive diffraction
 * Main reference: 0805.4071
 * 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2016-2019
 * 
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream> 
#include <iomanip>

#include <gsl/gsl_rng.h>

//#include <tools/tools.hpp>
#include <gsl/gsl_errno.h>
#include <subnucleon_config.hpp>

#include <dipole.hpp>
#include <smooth_ws_nuke.hpp>
#include <ipsat_proton.hpp>
#include <vector.hpp>
#include <ipglasma.hpp>
#include <nucleons.hpp>

#include "inclusive_diffraction.hpp"
#include "gitsha1.h"

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
    double beta = -1;   // if beta>0, use fixed beta and compute as a function of xpom
    double t=0.1;
    int A=1;
    DGLAPDist *gd=0;  // Initialized and used if we have nucleus consisting of ipsatnucleons (old ipsat)
    
    bool ms=false;
    bool smallb=false;
    double gbw = false;
    bool total=false;
    bool charm=false;
   
    
    cout << "# Inclusive Diffraction by H. Mäntysaari <heikki.mantysaari@jyu.fi>, 2016-2019" << endl;
    cout << "# Git version " << g_GIT_SHA1 << " local repo " << g_GIT_LOCAL_CHANGES << " main build " << __DATE__  << " " << __TIME__ << endl;
    cout << "# Command: ";
    for (int i=1; i<argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    
    if (string(argv[1])=="-help")
    {
        
        return 0;
    }
    
    Ipsat_version ipsatv;
    
    for (int i=1; i<argc; i++)
    {
       
        if (string(argv[i])=="-dipole")
        {
            A = StrToInt(argv[i+1]);
            if (A==1)
            {
                if (string(argv[i+2])=="ipsatproton" or string(argv[i+2])=="ipnonsatproton" or string(argv[i+2])=="ipsat12proton")
                {
                    if(string(argv[i+2])=="ipsatproton")
                        ipsatv = MZSAT;
                    else if(string(argv[i+2])=="ipnonsatproton")
                        ipsatv = MZNONSAT;
                    else if(string(argv[i+2])=="ipsat12proton")
                        ipsatv = IPSAT12;
                    else
                    {
                        cerr << "WTF ipsat: " << argv[i+2] << endl;
                        exit(1);
                    }
                    amp = new Ipsat_Proton(ipsatv);
                    
                    ((Ipsat_Proton*)amp)->SetProtonWidth(StrToReal(argv[i+3]));
                    ((Ipsat_Proton*)amp)->SetQuarkWidth(StrToReal(argv[i+4]));
                  
                    ((Ipsat_Proton*)amp)->SetShape(GAUSSIAN);
                    
                }
                else if (string(argv[i+2])=="ipglasma")
                    amp = new IPGlasma(argv[i+3]);
                else
                {
                    cerr << "Unknown dipole " << argv[i+2] << endl;
                    return -1;
                }
            }
            else // A>1
            {
                if(string(argv[i+2])=="ipsatproton")
                    ipsatv = MZSAT;
                else if(string(argv[i+2])=="ipnonsatproton")
                    ipsatv = MZNONSAT;
                else
                {
                    cerr << "Unknown dipole " << argv[i+2] << endl;
                    exit(1);
                }
                amp = new Smooth_ws_nuke(A, ipsatv);
                ((Smooth_ws_nuke*)amp)->SetSmoothApproximation(true);
            }
            
            
            
        }
        else if (string(argv[i])=="-Q2")
            Qsqr = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-xpom")
            xpom = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-beta")
            beta = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-ms")
            ms=true;
        else if (string(argv[i])=="-gbw")
            gbw=true;
        else if (string(argv[i])=="-smallb")
            smallb=true;
        else if (string(argv[i])=="-qq")
            gbw=false;
        else if (string(argv[i])=="-total")
            total=true;
        else if (string(argv[i])=="-charm")
            charm=true;
        
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
    
    cout <<"# Dipole amplitude initialized, N(r=1 GeV^-1, xp=0.001,b=0)=" << amp->Amplitude(0.001, Vec(-0.5,0), Vec(0.5,0)) << endl;
    
    
    cout << "# " << InfoStr() ;
    
        

    
    InclusiveDiffraction diffraction(amp, ipsatv, charm, A);
    
    if (charm)
        cout << "# Including only charm quarks" << endl;
    else
        cout << "# Quarks: u,d,s,c" << endl;
    
    if (total)
    {
        cout <<"# Q^2=" << Qsqr << endl;
        cout << "#xp diffractive inclusive ratio" << endl;
        for (double xp = 1e-12; xp < 0.02; xp*=1.2*1.2*1.2)
        {
            double inclusive = diffraction.TotalInclusive_qq(xp,Qsqr);
            
            double diffractive =diffraction.TotalDiffractive_qq(xp,Qsqr);
            
            
            cout << xp << " " << diffractive << " " << inclusive << " " << diffractive/inclusive << endl;
        }
        
        
        return 0;
    }
    
    
    
    if (beta > 0)
    {
        cout <<"# Q^2=" << Qsqr <<", beta=" << beta << endl;
        cout << "# xp transverse longitudinal" << endl;
        for (double xp = 1e-12; xp < 0.02; xp*=1.2*1.2*1.2)
        {
            double f_qq_t = diffraction.DiffractiveStructureFunction_qq_T(xp, beta, Qsqr);
            double f_qq_l = diffraction.DiffractiveStructureFunction_qq_L(xp, beta, Qsqr);
            cout << xp << " " << f_qq_t << " " << f_qq_l << endl;
        }
        
        
        return 0;
    }
    
    
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
    
    
    for (double beta=0.02; beta<1; beta+=0.04)
    {
        if (gbw)
        {
            double gbw = diffraction.DiffractiveStructureFunction_qqg_GBW_T(xpom, beta, Qsqr);
            cout << beta << " " <<gbw << endl;
        }
        else
        {
            double f_qq_t = diffraction.DiffractiveStructureFunction_qq_T(xpom, beta, Qsqr);
            double f_qq_l = diffraction.DiffractiveStructureFunction_qq_L(xpom, beta, Qsqr);
            cout << beta << " " << f_qq_t << " " << f_qq_l << endl;
        }
        
       
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

