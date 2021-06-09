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
/* add F1F221Wrapper
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


#ifndef _PBosted_H
#define _PBosted_H

namespace PBosted 
{

	//input
	//Z,N: proton and neutron number of the nucleus.;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;
	//bVersion21: flag to tell whether to use F1F221 or F1F209;

	double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb=-0.001, double Ta=-0.001, bool bVersion21=true);
	double GetXS21(double Z, double A, double Ei, double Ef, double theta);  //use F1F221.f
}
#endif

