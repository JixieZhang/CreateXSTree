//2009 version of P. Bosted's electron XS model
//
//Fits to inclusive inelastic electron scattering (see M.E. Christy and P.E. Bosted, 
//arXiv:0711.0159 (submitted to PRC), for proton fit, and P.E. Bosted and M.E. Christy, 
//arXiv:0711.0159 [also Phys. Rev. C 77, 065206 (2008)], for deuteron and neutron fit.
// This 2007 code also includes the quasi-elastic and inelastic model for nuclei shown 
//in Phys. Rev. C. 78, 015202 (2008). (arXiv:0712.2438) (Ratios of 15N/12C and 4He/12C 
//inclusive electroproduction cross sections in the nucleon resonance region, 
//P.E. Bosted, R. Fersch {\it et al.}). 
//The 2009 fit is the same for proton and deuteron as the 2007 fit, but for A>2 the 
//fit has been greatly improved by the addition of new data from JLab, and an additional 
//25 parameters for A-dependence. This code works well up to copper, but may have problems 
//for heavier nuclei. It works well for 4He, but the model for 3He still needs tweaking. 
//The reference for the 2009 $A>2$ fit is: http://arxiv.org/abs/1203.2262. 
//
/* F1F221Wrapper
 * Unit of xs is ub/MeV-sr.
 * F1F209.f is valid for all W<5 GeV and all Q2<11 GeV2.
 * Inside the code of F1F209.f, I saw that this fit does not work for Q2>11. 
 * https://github.com/jixie/CreateXSTree/blob/main/XSModel/Bosted/bosted.f#L1010
 * For A>2, it used function fitemc() , which is a fit to "EMC effect". 
 * fitemc() marks a fit to good_fit only if 0.0085< Xbj <0.7.
 * https://github.com/jixie/CreateXSTree/blob/main/XSModel/Bosted/bosted.f#L1696
 *
 * In F1F221.f, Q2 and W2 have been extended up to 32.
 * For A > 2, it does not use fitemc() any longer. It uses the parameterization  based on
 * 12C, 27Al, 56Fe and 64Cu.
 * It claims that this fit scales relatively well for all nuclei with 10 < A < 80.
 *
 * Use -fno-leading-underscore -fno-second-underscore when compiling F1F221.f
 */
#include <math.h>
#include <iostream>

namespace PBosted 
{
	static const bool kUseF1F2Version21 = true;

	extern "C" 
	{
		void bosted_(double* Z, double* A, double* Ei, double* Ep, double* ang, double* xs, double *Tb, double *Ta);
		
		void f1f2in21_(double* Z, double* A, double* Q2, double* W2, double* F1, double* F2);
		void f1f2qe21_(double* Z, double* A, double* Q2, double* W2, double* F1, double* F2);
	}


	double GetXS21(double Z, double A, double Ei, double Ef, double theta)
	{
		const double M = 0.93825;
		theta = fabs(theta);  //in case the theta is negative 
		double halftheta = theta/2;

		double nu = Ei - Ef;
		double Q2 = 4. * Ei * Ef * pow(sin(halftheta), 2.0);
		double w2 = M * M + 2. * M * nu - Q2;

		double F1, F2;
		double xs1, xs2;

		//original code
      	//xs2=(2./137.*Ep/q2*cos(theta/2))**2. !mott
      	//xs2=xs2*(2/M*F1*tan(abs(theta)/2)**2+F2/nu)
      	//xs2=xs2*389.379

		f1f2in21_(&Z, &A, &Q2, &w2, &F1, &F2);
		xs1 = pow(2. / 137. * Ef / Q2 * cos(halftheta), 2.0); // mott
		xs1 = xs1 * (2. / M * F1 * pow(tan(halftheta), 2.0) + F2 / nu);
		xs1 = xs1 * 389.379;

		f1f2qe21_(&Z, &A, &Q2, &w2, &F1, &F2);
		xs2 = pow(2. / 137. * Ef / Q2 * cos(halftheta), 2.0); // mott
		xs2 = xs2 * (2. / M * F1 * pow(tan(halftheta), 2.0) + F2 / nu);
		xs2 = xs2 * 389.379;

		//double deg = atan(1.0)/45.;
		//std::cout<<" Ei="<<Ei<<" Theta="<<theta/deg<<" Ef="<<Ef<<"  XS_in="<<xs1<<"  XS_qe="<<xs2<<std::endl; 
		return (xs1 + xs2) / 1000.; // ub/MeV-sr
	}


	//input
	//       Z,N: proton and neutron number of the nucleus.	;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;
	double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb, double Ta)
	{
		double NZ, NA;
		NZ = Z;
		NA = Z+N;
		double XS;
#ifdef WIN32
		//this fortran routine does not work in windows, do not know why
		return 1.0;
#else
		if(kUseF1F2Version21) {
			XS = GetXS21(NZ, NA, Ei, Ef, theta);
		}
		else {
			bosted_(&NZ, &NA, &Ei, &Ef, &theta, &XS, &Tb, &Ta);  //rad corr does not work in this version, no time to check it yet
		}
#endif
		return XS;
	}
}

