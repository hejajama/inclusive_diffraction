#ifndef INCLUSIVE_DIFFRACTION_HPP
#define INCLUSIVE_DIFFRACTION_HPP
/*
 * Inclusive diffraction
 * Main reference: 0805.4071
 * 
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2016
 * 
 */
 
 


#include "dipole.hpp"
#include <vector.hpp>
#include <vector>

const double ALPHAs=0.3;

class InclusiveDiffraction
{
public:
    InclusiveDiffraction(DipoleAmplitude* amp);
    
    // qq component 
    double DiffractiveStructureFunction_qq_T(double xpom, double beta, double qsqr);
    double DiffractiveStructureFunction_qq_L(double xpom, double beta, double qsqr);
    
    // qqg, two approximations
    double  DiffractiveStructureFunction_qqg_GBW_T(double xpom, double beta, double qsqr);
    double  DiffractiveStructureFunction_qqg_MS_T(double xpom, double beta, double qsqr);
    
    // Interpolation formula for qqg
    double  DiffractiveStructureFunction_qqg_T(double xpom, double beta, double qsqr);
    
    double Qq_component_n(double xpom, double qsqr, double Mxsqr, double z, int n, int flavor=0	); // Eq. 7
    
    double A(double r, double xpom, Vec b);	// Eq. 13
    double A_bint(double r, double xpom);	// Eq. 13 integrated over d^2b
     
    double QuarkMass(int i);
    int NumberOfQuarks();
    
private:
    DipoleAmplitude* amplitude;
    
    
    std::vector<double> m_f;	// Quark amsses
    
};

struct inthelper_inclusive
{
	DipoleAmplitude* amp;
	InclusiveDiffraction* diffraction;
	double xpom;
	double beta;
	double qsqr;
	double r;
	double Mxsqr;
	double z;
	int flavor;	// quark flavor sum
	int bessel_component;
	Vec b;
	double b_len;
	double b_theta;
	double ksqr;
	double r2;
};

#endif /* smooth_ws_nuke_hpp */

