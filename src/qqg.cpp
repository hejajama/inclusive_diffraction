/*
 * Inclusive diffraction
 * Main reference: 0805.4071
 * 
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016-2025
 * 
 */
 
 // Note, this does not support IPglasma, as we do not keep track of the 
 // dipole orientation
 // Monte carlo ipsat is NOT supported, as we have 2d vector to describe b
 
#include "inclusive_diffraction.hpp"
#include <gsl/gsl_integration.h>
#include <vector.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <cmath>
#include <iostream>
#include <virtual_photon.hpp>
using namespace std;

const int INTERVALS=16;
const double ACCURACY=0.001;
const double MAXR=100;

gsl_integration_workspace *gsl_wp_kint;
gsl_integration_workspace *gsl_wp_bint;
gsl_integration_workspace *gsl_wp_rint;
gsl_integration_workspace *gsl_wp_thetaint;

/*
 * Helper function to calculate Qq component
 * Eq. 10
 */
double inthelperf_Qqg_component_r(double r, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	double z = par->z;
	double Q2 = par->qsqr;
	double Mxsqr = par->Mxsqr;
	double ksqr = par->ksqr;
	double k = sqrt(ksqr);
	
	// DipoleAmplitude requires positions of the quark and the antiquark
	// We assume here that it uses IPsat, so dipole orientation does not matter
	Vec q1 = par->b;
	Vec q2 = par->b;
	q1.SetX(q1.GetX()-r/2.0);
	q2.SetX(q2.GetX()+r/2.0); // Dipole center is now b
	
	double dsigma = 2.0 * par->amp->Amplitude(par->xpom, q1, q2);	// 2 as this is sigma_qq, not N
    double dsigma_adj;
    if (par->ipsat == MZNONSAT)
        dsigma_adj = 2.0 * dsigma; // N_adj = 2N - N^2 \approx 2N
    else
        dsigma_adj = 2.0 * ( 1.0 - pow( 1.0 - 0.5*dsigma, 2.0 ) );

	double K2 = gsl_sf_bessel_Kn(2, sqrt(z)*k*r);

	double J2 = gsl_sf_bessel_Jn(2, sqrt(1.0-z)*k*r);

	double result = r*J2*K2*dsigma_adj;

	return result;
}

double inthelperf_Qqg_component_k(double ksqr, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->ksqr=ksqr;
	gsl_function f;
    f.params = par;
    
    f.function = inthelperf_Qqg_component_r;
    //gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    double result,error;
    int status = gsl_integration_qag(&f, 0, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, gsl_wp_rint, &result, &error);
    
    // cout << "rint at k " << par->k << " res " << result << " error " << error << endl;
     
    //if (status)
    //    cerr << "#r failed, result " << result << " relerror " << error << " b " << par->b_theta << endl;
    
    //gsl_integration_workspace_free(w);
    
    return result*result * ksqr*ksqr*log( par->qsqr/ksqr) * ( pow(1.0-par->beta/par->z, 2.0) + pow(par->beta/par->z, 2.0) );
    
}

double inthelperf_Qqg_component_n_b_theta(double b_theta, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->b_theta = b_theta;
	par->b.SetX(par->b_len*cos(par->b_theta));
	par->b.SetY(par->b_len*sin(par->b_theta));
	gsl_function f;
    f.params = par;
    
    
    
    f.function = inthelperf_Qqg_component_k;
    //gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    double result,error;
    int status = gsl_integration_qag(&f, 0.001, par->qsqr, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, gsl_wp_kint, &result, &error);
    
    
    
    //cout << "Kint at b " << par->b << " res " << result << " error " << error << endl;
    //if (status)
    //    cerr << "#r failed, result " << result << " relerror " << error << " b " << par->b_theta << endl;
    
    //gsl_integration_workspace_free(w);
    
    return result;
}


double inthelperf_Qqg_component_n_b(double b, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->b_len = b;
	gsl_function f;
    f.params = par;
    
    // Skip b_theta integral as we have symmetrical system now
    par->b_theta = 0;
    return 2.0*M_PI*b*inthelperf_Qqg_component_n_b_theta(0, par);
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



    
// qqg, two approximations
// GBW

double inthelperf_zint_gbw(double z, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	
	par->z=z;
	
    if (par->ipsat == MZNONSAT or par->ipsat==MZSAT)
    {
        // Do b integral analytically
        // Factorize N(r,b) = N(r) e^(-b^2/(2B))
        // b integral gives \int d^2 b [e^(-b^2/(2B))]^2 = pi * B_p
        
        // This is possible as we assume here that in ipnonsat model the adjoint dipole
        // N_adj = 2*N_fund
        
        // Effectively evaluate N(r,b) at b=0

        // Note that b dependence does not actually factor out from IPsat!
        double B_p=4.0;
        par->b_len = 0;
        return M_PI*B_p*inthelperf_Qqg_component_n_b_theta(0, par);
    }
    
    else
    {
        gsl_function f;
        f.params = par;
        
        f.function = inthelperf_Qqg_component_n_b;
        //gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
        double result,error;
        int status = gsl_integration_qag(&f, 0, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, gsl_wp_bint, &result, &error);
        
        if (status)
            cerr << "#bint failed, result " << result << " relerror " << error << endl;
        
        //gsl_integration_workspace_free(w);
        
        return result;
    }
}

double InclusiveDiffraction::DiffractiveStructureFunction_qqg_GBW_T(double xpom, double beta, double qsqr)
{
	inthelper_inclusive par;
	par.diffraction = this;
	par.amp=amplitude;
	par.xpom=xpom;
	par.beta=beta;
	par.qsqr=qsqr;
    par.ipsat=ipsat;
	
	double mxsqr = qsqr / beta - qsqr;
	par.Mxsqr = mxsqr;
	

    // Init memory
    gsl_wp_bint = gsl_integration_workspace_alloc(INTERVALS);
    gsl_wp_kint = gsl_integration_workspace_alloc(INTERVALS);
    gsl_wp_rint = gsl_integration_workspace_alloc(INTERVALS);
    
    gsl_function f;
    f.params = &par;
    f.function = inthelperf_zint_gbw;


    par.flavor=0;   // Not used
    
    double z0 = beta;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
    double result,error;
    int status = gsl_integration_qag(&f, z0, 1.0, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
    
    //cout << "zint from " << z0 << " to 1/2: " << result << " relerr " << error/result << endl;
    
    if (status)
        cerr << "#z from " << z0 << " to 1 failed, result " << result << " relerror " << error/result  << endl;
    
    gsl_integration_workspace_free(w);
    gsl_integration_workspace_free(gsl_wp_bint);
    gsl_integration_workspace_free(gsl_wp_rint);
    gsl_integration_workspace_free(gsl_wp_kint);
    
    cout << "# GBW integral,  beta " << beta << " result " << result << " relerror " << error/result << endl;
    
    
    // e_f^2 sum
    double efsum=0;
    for (unsigned int i=0; i<NumberOfQuarks(); i++)
    {
        efsum += e_f[i]*e_f[i];
    }
    
    
    
    return ALPHAs*beta / (8.0*pow(M_PI,4.0)) *  result*efsum;
}


double inthelperf_MS_r(double r, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->r = r;
	
	
    VirtualPhoton photon; // default uds quarks
    double wf_sqr = photon.PsiSqr_T_intz(par->qsqr, r);
        
    
	return r*wf_sqr*par->diffraction->A_bint(r, par->xpom);
}

double InclusiveDiffraction::DiffractiveStructureFunction_qqg_MS_T(double xpom, double beta, double qsqr)
{
	gsl_wp_thetaint = gsl_integration_workspace_alloc(INTERVALS);
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
	
	inthelper_inclusive par;
	par.diffraction=this;
	par.amp=amplitude;
	par.xpom=xpom;
	par.qsqr = qsqr;
	par.beta=beta;
	
	gsl_function f;
	f.params=&par;
	f.function = inthelperf_MS_r;
	double result,error;
    int status = gsl_integration_qag(&f, 1e-10, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);
    
	
	double CF = (3.0*3.0-1.0)/(2.0*3.0);
	double ALPHA_em = 1.0/137.0;
	
	gsl_integration_workspace_free(gsl_wp_thetaint);
	gsl_integration_workspace_free(w);
	
	
	return CF * ALPHAs * qsqr / (4.0*M_PI*M_PI*ALPHA_em) * result;
	
	
	
	
	
}
    
// Interpolation formula for qqg
double InclusiveDiffraction::DiffractiveStructureFunction_qqg_T(double xpom, double beta, double qsqr)
{
	
	return 0;
}

// Eq. (13)

double inthelperf_A_r2_theta(double theta, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	
	// r, r2, r-r2 dipoles, with given impact parameter b
	// As we assume that the target dipole amplitude does not depend on angles, 
	// we can set vec b = (b,0)

    // Note that when using a factorized b dependence, this is evaluated with b=0
	Vec q1(0,0);
    Vec q2(par->r,0); 
    Vec q3(par->r2*cos(theta), par->r2*sin(theta));

    q1 = q1+par->b;
    q2 = q2+par->b;
    q3 = q3+par->b;

    Vec r = q1-q2;
    Vec rp = q1-q3;
    Vec r_m_rp = q2-q3;

	
    double nr2 = par->amp->Amplitude(par->xpom, q1,q3);	
    double nr1 = par->amp->Amplitude(par->xpom, q2, q3);
    double nr0 = par->amp->Amplitude(par->xpom, q1, q2);
	
	double kernel = r.LenSqr() / ( rp.LenSqr()*r_m_rp.LenSqr() + 1e-40);
    
    // nonlinear term, ipnonsat: consistent to drop (?)
    double nonlinear =nr1*nr2;
    if (par->ipsat == MZNONSAT)
        nonlinear = 0;
    
	
	double dipole = pow(nr1 + nr2 - nr0 - nonlinear,2.0);
	
	
	//cout << kernel << " " << dipole <<" " << par->r << " " << par->r2 <<endl;
	
	return par->r2*kernel*dipole;
	
	
	
}
double inthelperf_A_r2(double r2, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	par->r2=r2;
	
	gsl_function f;
	f.function = inthelperf_A_r2_theta;
	f.params = par;
	double result,error;
    int status = gsl_integration_qag(&f, 0, 2.0*M_PI, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, gsl_wp_thetaint, &result, &error);

    //if (status and result>1e-20)
    //    cerr << "#r2theta failed, result " << result << " relerror " << error  << endl;
	
	
	return result;
}


double InclusiveDiffraction::A(double r, double xpom, Vec b)
{
	inthelper_inclusive par;
	par.amp=amplitude;
	par.diffraction=this;
	par.xpom=xpom;
	par.b=b;
	par.b_len = b.Len();
	par.r=r;
	
	
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
	gsl_function f;
	f.function = inthelperf_A_r2;
	f.params = &par;
	double result,error;
    int status = gsl_integration_qag(&f, 1e-10, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);

    //if (status and abs(result)>1e-20)
    //    cerr << "#r2int failed, result " << result << " relerror " << error  << endl;
	
	gsl_integration_workspace_free(w);
	
	
	return result;
	
}

double inthelperf_A_b(double b, void* p)
{
	inthelper_inclusive* par = (inthelper_inclusive*)p;
	Vec bvec(b,0);
	return 2.0*M_PI*b*par->diffraction->A(par->r, par->xpom, bvec);
}

double InclusiveDiffraction::A_bint(double r, double xpom)
{
	inthelper_inclusive par;
	par.amp=amplitude;
	par.diffraction=this;
	par.xpom=xpom;
	par.r=r;

    // Factorized b dependence. NOTE: Currently the code assumes that for ipsat also, 
    // although T(b) dependence does not actually factor out 
    if (ipsat == MZSAT or ipsat==MZNONSAT)
    {
        Vec bvec(0,0);
        
        // ¨ do b integral for all N(r,b) by pulling out e^(-b^2/(2B))
        // Effectively we evaluate A at b=0 and multiply the result by \int d^2 [e^(-b^2/(2B))]^2 = pi*B
        return M_PI*4.0*A(r, xpom, bvec);
        
    }
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTERVALS);
	gsl_function f;
	f.function = inthelperf_A_b;
	f.params = &par;
	double result,error;
    int status = gsl_integration_qag(&f, 0, MAXR, 0, ACCURACY, INTERVALS, GSL_INTEG_GAUSS51, w, &result, &error);

    if (status and result>0.0001)
        cerr << "#A_bint failed, result " << result << " relerror " << error  << endl;
	
	gsl_integration_workspace_free(w);
	
	//cout << "A_bint at r=" << r <<" gives " << result << " pm " << error << endl;
	
	return result;
}
