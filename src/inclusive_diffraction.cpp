/*
 * Inclusive diffraction
 * Main reference: 0805.4071
 * 
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016
 * 
 */
 
 // Note, this does not support IPglasma, as we do not keep track of the 
 // dipole orientation
 // Monte carlo ipsat NOT FULLY supported, as we have 2d vector to describe b
 
#include "inclusive_diffraction.hpp"
#include <gsl/gsl_integration.h>
#include <vector.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <cmath>

using namespace std;
 

const double MAXR = 100;
const int INTERVALS = 5;
const double ACCURACY = 0.01;

 
InclusiveDiffraction::InclusiveDiffraction(DipoleAmplitude* amp)
{
	amplitude=amp;
	m_f.push_back(0.14); m_f.push_back(0.14); m_f.push_back(0.14); m_f.push_back(1.4);
    e_f.push_back(2.0/3.0); e_f.push_back(1.0/3.0); e_f.push_back(1.0/3.0), e_f.push_back(2.0/3.0);
}
    
 // qq component 
double inthelperf_zint_t(double z, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	
    int flavor=par->flavor;
	
	double phi1 = par->diffraction->Qq_component_n(par->xpom, par->qsqr, par->Mxsqr, z, 1, flavor);
	double phi0 = par->diffraction->Qq_component_n(par->xpom, par->qsqr, par->Mxsqr, z, 0, flavor);
	
	double mf=par->diffraction->QuarkMass(flavor);
	double eps = sqrt(z*(1.0-z)*par->qsqr + mf*mf);
	
	return z*(1.0-z)*(eps*eps*(z*z + pow(1.0-z, 2.0))*phi1 + mf*mf*phi0);
	
}

double inthelperf_zint_l(double z, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	
    int flavor=par->flavor;
	
	double phi0 = par->diffraction->Qq_component_n(par->xpom, par->qsqr, par->Mxsqr, z, 0, flavor);
	
	double mf=par->diffraction->QuarkMass(flavor);
	double eps = sqrt(z*(1.0-z)*par->qsqr + mf*mf);
	
	return pow(z*(1.0-z),3.0)*phi0;
	
}
	
double InclusiveDiffraction::DiffractiveStructureFunction_qq_T(double xpom, double beta, double qsqr)
{
	inthelper_inclusive par;
	par.diffraction = this;
	par.xpom=xpom;
	par.beta=beta;
	par.qsqr=qsqr;
	double mxsqr = qsqr / beta - qsqr;
	par.Mxsqr = mxsqr;
    double sum=0;
    for (unsigned int flavor=0; flavor<NumberOfQuarks(); flavor++)
    {
        
        gsl_function f;
        f.params = &par;
        f.function = inthelperf_zint_t;
        
        par.flavor=flavor;
        
        double z0 = (1.0 - sqrt(1.0 - 4.0*m_f[flavor]*m_f[flavor]/mxsqr))/2.0;
        
        f.function = inthelperf_zint_t;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
        double result,error;
        int status = gsl_integration_qag(&f, z0, 0.5, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        cout << "# Transverse, flavor " << flavor << "contribution w.o. quark charge " << result*3.0*qsqr*qsqr/(16.0*pow(M_PI,3.0)*beta) << endl;
        //cout << "zint from " << z0 << " to 1/2: " << result << " relerr " << error/result << endl;
        
        if (status)
            cerr << "#z failed, result " << result << " relerror " << error  << endl;
        
        gsl_integration_workspace_free(w);
        
        sum += result*e_f[flavor]*e_f[flavor];
    }
    
    
    return 3.0*qsqr*qsqr/(16.0*pow(M_PI,3.0)*beta) *  sum;
    
}

double InclusiveDiffraction::DiffractiveStructureFunction_qq_L(double xpom, double beta, double qsqr)
{
	inthelper_inclusive par;
	par.diffraction = this;
	par.xpom=xpom;
	par.beta=beta;
	par.qsqr=qsqr;
	double mxsqr = qsqr / beta - qsqr;
	par.Mxsqr = mxsqr;
	
	
    
    double sum=0;
    for (unsigned int flavor=0; flavor<NumberOfQuarks(); flavor++)
    {
        gsl_function f;
        f.params = &par;
        f.function = inthelperf_zint_t;
        
        par.flavor=flavor;
    
        double z0 = (1.0 - sqrt(1.0 - 4.0*m_f[flavor]*m_f[flavor]/mxsqr))/2.0;
        
        f.function = inthelperf_zint_l;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
        double result,error;
        int status = gsl_integration_qag(&f, z0, 0.5, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
        
        //cout << "zint from " << z0 << " to 1/2: " << result << " relerr " << error/result << endl;
        
        if (status)
            cerr << "#z failed, result " << result << " relerror " << error  << endl;
        
        gsl_integration_workspace_free(w);
        
        sum += result * e_f[flavor]*e_f[flavor];
    }
    
    return 3.0*qsqr*qsqr*qsqr/(4.0*pow(M_PI,3.0)*beta) *  sum;
}


/*
 * Helper function to calculate Qq component
 * Eq. 7
 */
double inthelperf_Qq_component_n(double r, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	double z = par->z;
	double Q2 = par->qsqr;
	double Mxsqr = par->Mxsqr;
	
	// DipoleAmplitude requires positions of the quark and the antiquark
	// We assume here that it uses IPsat, so dipole orientation does not matter
	Vec q1 = par->b;
	Vec q2 = par->b;
	q1.SetX(q1.GetX()-r/2.0);
	q2.SetX(q2.GetX()+r/2.0); // Dipole center is now b
	double dsigma = 2.0 * par->amp->Amplitude(par->xpom, q1, q2);	// 2 as this is sigma_qq, not N

	double mf=par->diffraction->QuarkMass(par->flavor);
	double eps = sqrt(z*(1.0-z)*Q2 + mf*mf);
	double Kn = gsl_sf_bessel_Kn(par->bessel_component, eps*r);
	double k = sqrt(z*(1.0-z)*Mxsqr-mf*mf);
	double Jn = gsl_sf_bessel_Jn(par->bessel_component, k*r);
	double result = r*Jn*Kn*dsigma;
	
	
	return result;
}


double inthelperf_Qq_component_n_b_theta(double b_theta, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->b_theta = b_theta;
	par->b.SetX(par->b_len*cos(par->b_theta));
	par->b.SetY(par->b_len*sin(par->b_theta));
	gsl_function f;
    f.params = par;
    
    f.function = inthelperf_Qq_component_n;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
    
    //if (status)
    //    cerr << "#r failed, result " << result << " relerror " << error << " b " << par->b_theta << endl;
    
    gsl_integration_workspace_free(w);
    
    return result*result;
}


double inthelperf_Qq_component_n_b(double b, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->b_len = b;
	gsl_function f;
    f.params = par;
    
    // Skip b_theta integral as we have symmetrical system now
    par->b_theta = 0;
    return 2.0*M_PI*b*inthelperf_Qq_component_n_b_theta(0, par);
    /*
    f.function = inthelperf_Qq_component_n_b_theta;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#thetaint failed, result " << result << " relerror " << error << " b " << b << endl;
    
    gsl_integration_workspace_free(w);
    
    return result;
    */
}




double InclusiveDiffraction::Qq_component_n(double xpom, double qsqr, double Mxsqr, double z, int n, int flavor)
{
	inthelper_inclusive par;
	par.amp=amplitude;
	par.xpom=xpom;
	par.diffraction = this;
	par.qsqr = qsqr;
	par.Mxsqr = Mxsqr;
	par.flavor = flavor;
	par.bessel_component = n;
	par.z=z;
	
	gsl_function f;
    f.params = &par;
    
    f.function = inthelperf_Qq_component_n_b;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
    
    if (status)
        cerr << "#thetaint failed, result " << result << " relerror " << error << endl;
    
    gsl_integration_workspace_free(w);
    
    //cout << "Helper with n=" << n << " z " << z << " mxsqr " << Mxsqr << " res " << result << " err " << error << endl;
	
	return result;
	
}



int InclusiveDiffraction::NumberOfQuarks()
{
	return m_f.size();
}

double InclusiveDiffraction::QuarkMass(int i)
{
	return m_f[i];
}

double InclusiveDiffraction::QuarkCharge(int i)
{
    return e_f[i];
}
