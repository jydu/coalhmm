//
// File: IlsCoalHmmTransitionMatrix09.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 04 16:20 2008
// NB: this replaces former file with only one recombination rate.
//

//This file is part of the CoalHMM program and library.
//
//CoalHMM is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//any later version.
//
//CoalHMM is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with CoalHMM.  If not, see <http://www.gnu.org/licenses/>.

#include "IlsCoalHmmTransitionMatrix09.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;
using namespace bpp;

const short IlsCoalHmmTransitionMatrix09::ONE_RHO = 1;
const short IlsCoalHmmTransitionMatrix09::TWO_RHOS_12 = 2;
const short IlsCoalHmmTransitionMatrix09::THREE_RHOS = 3;

IlsCoalHmmTransitionMatrix09::IlsCoalHmmTransitionMatrix09(const ThreeSpeciesCoalHmmStateAlphabet* hiddenAlphabet, double rho1, double rho2, double rho3, double lowerRho, double upperRho):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet), rhoOption_(THREE_RHOS), rho1_(rho1), rho2_(rho2), rho3_(rho3),
  s_(), u_(), v1_(), v2_()
{
  //We set a maximum bound for practical use.
  lowerRho = std::max(0., lowerRho);
  upperRho = std::max(lowerRho + 0.1, upperRho);
  addParameter_(new bpp::Parameter("coal.rho1", rho1, new bpp::IntervalConstraint(lowerRho, upperRho, true, true), true));
  addParameter_(new bpp::Parameter("coal.rho2", rho2, new bpp::IntervalConstraint(lowerRho, upperRho, true, true), true));
  addParameter_(new bpp::Parameter("coal.rho3", rho3, new bpp::IntervalConstraint(lowerRho, upperRho, true, true), true));
  unsigned int n = getNumberOfStates();
  transitions_.resize(n, n);
  freqs_.resize(n);
  actualizeTransitionMatrix();
}

IlsCoalHmmTransitionMatrix09::IlsCoalHmmTransitionMatrix09(const ThreeSpeciesCoalHmmStateAlphabet* hiddenAlphabet, double rho12, double rho3, short rhoOption, double lowerRho, double upperRho):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet), rhoOption_(rhoOption), rho1_(), rho2_(), rho3_(),
  s_(), u_(), v1_(), v2_()
{
  if (rhoOption == TWO_RHOS_12)
  {
    rho1_ = rho12;
    rho2_ = rho12;
    rho3_ = rho3;
  }
  else throw bpp::Exception("IlsCoalHmmTransitionMatrix08 (constructor 2). Bad rho option.");
  //We set a maximum bound for practical use.
  addParameter_(new bpp::Parameter("coal.rho12", rho12, new bpp::IntervalConstraint(lowerRho, upperRho, true, true), true));
  addParameter_(new bpp::Parameter("coal.rho3" , rho3, new bpp::IntervalConstraint(lowerRho, upperRho, true, true), true));
  unsigned int n = getNumberOfStates();
  transitions_.resize(n, n);
  freqs_.resize(n);
  actualizeTransitionMatrix();
}

IlsCoalHmmTransitionMatrix09::IlsCoalHmmTransitionMatrix09(const ThreeSpeciesCoalHmmStateAlphabet* hiddenAlphabet, double rho, double lowerRho, double upperRho):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet), rhoOption_(ONE_RHO), rho1_(rho), rho2_(rho), rho3_(rho),
  s_(), u_(), v1_(), v2_()
{
  //We set a maximum bound for practical use.
  addParameter_(new bpp::Parameter("coal.rho", rho, new bpp::IntervalConstraint(lowerRho, upperRho, true, true), true));
  unsigned int n = getNumberOfStates();
  transitions_.resize(n, n);
  freqs_.resize(n);
  actualizeTransitionMatrix();
}


double phc1(double tau2, double thtHC)
{
  return 1-exp(-tau2/thtHC);
}
 
double phc2(double tau2, double thtHC)
{
  return exp(-tau2/thtHC)/3.;
}
 
double phg(double tau2, double thtHC)
{
  return exp(-tau2/thtHC)/3.; 
}
 
double pcg(double tau2, double thtHC)
{
  return exp(-tau2/thtHC)/3.;
}

double expRho(double rho)
{
  return exp(rho);
}

double expRhoTau(double rho, double tau)
{
  return exp(rho * tau);
}

double expTauOverTheta(double tau, double theta)
{
  return exp(tau/theta);
}

double IlsCoalHmmTransitionMatrix09::threeS(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC)
{
  double x = 1./expTauOverTheta(tau2, thtHC);
  double y = (1./expRhoTau(rhoH, tau1+tau2)) * (1./expRhoTau(rhoM, tau1+tau2));
  double z = expTauOverTheta(tau2, thtHC);
  double w = expRhoTau(rhoH, tau2) * expRhoTau(rhoM, tau2);
  double u = expTauOverTheta(tau2, thtHC) - 1;
  double v = thtHC*(rhoH + rhoM) - 1;
  return x * (1 + (y * (z - w))/(u*v));
}

double IlsCoalHmmTransitionMatrix09::u(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC)
{
  double x = 1./expTauOverTheta(tau2, thtHC); 
  double y = (1. / expRhoTau(rhoH, tau1+tau2)) * (1./expRhoTau(rhoM, tau1+tau2));
  double z = (1./expTauOverTheta(tau2, thtHC)) * (1./expRhoTau(rhoH, tau1)) * (1./expRhoTau(rhoM, tau1));
  double u = thtHC*(rhoH + rhoM) - 1;
  return 1 - x + y/u - z/u;
}

double e001a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = (1./expTauOverTheta(tau2, thtHC)) * (1./expRhoTau(rhoH, tau1 + tau2)) * (1./expRhoTau(rhoM, tau1 + tau2));
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = 3./(2. * (3. + thtHCG*(rhoH + rhoM)));
  double w = 1./(2. * (5. + thtHCG*(rhoH + rhoM)));
  return (1./3.) * x * y * (z - w); 
} 

double e001b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double z = thtHCG*rhoG;
  double u = 6. + thtHCG*(rhoH + rhoM);
    
  double v1 = 3. + thtHCG*(rhoH + rhoM); 
  double v2 = 5. + thtHCG*(rhoH + rhoM);
  double v3 = 3. + thtHCG*(rhoG + rhoH + rhoM);

  double nr = x * y1 * y2 *y3 * z * u;
  double dr = 3. * v1 * v2 * v3;

  return (nr/dr);
} 

double e001(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e001a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e001b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y);
}
    
double hrec(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = pow(expTauOverTheta(tau2, thtHC),-2.);
  double y = 1 - 1./expRhoTau(rhoH, tau1);
  double z = (1./expRhoTau(rhoH, tau1)) * pow(expTauOverTheta(tau2, thtHC),-2);
  double w = (expTauOverTheta(tau2, thtHC) * 1./expRhoTau(rhoH, tau2)) - 1;
  double u = thtHC * rhoH;
  return (x*y - (z*w*u)/(u-1));
}

double mrec(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = pow(expTauOverTheta(tau2, thtHC), -2);
  double y = 1 - 1./expRhoTau(rhoM, tau1);
  double z = (1./expRhoTau(rhoM, tau1)) * pow(expTauOverTheta(tau2, thtHC), -2);
  double w = (expTauOverTheta(tau2, thtHC) * 1./expRhoTau(rhoM, tau2)) - 1;
  double u = thtHC * rhoM;
  return (x*y - (z*w*u)/(u-1));
}

double e100x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = (1./expRhoTau(rhoG, tau1 + tau2)) * (1./expRhoTau(rhoM, tau1 + tau2)) ;
  double y = 3. + thtHCG*(rhoM + rhoG);
  double z = 5. + thtHCG*(rhoM + rhoG);
  return (x/(y*z));
}

double e100a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e100x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e100b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z = thtHCG * rhoH;

  double u1 = 3. + thtHCG*(rhoG + rhoM);
  double u2 = 5. + thtHCG*(rhoG + rhoM); 
  double u3 = 3. + thtHCG*(rhoG + rhoH + rhoM);
    
  double nr = x * y1 * y2 * y3 * z;
  double dr = u1 * u2 * u3;

  return (nr/dr);
}

double e100(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e100a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e100b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y);
} 

double e010x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y = 5. + thtHCG*(rhoH + rhoG);
  return ((x1 * x2)/(3.*y));
}

double e010a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e010x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e010b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2); 
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z = thtHCG * rhoM;

  double u1 = 5. + thtHCG*(rhoG + rhoH);
  double u2 = 3. + thtHCG*(rhoG + rhoH + rhoM);
    
  double nr = x * y1 * y2 * y3 * z;
  double dr = 3. * u1 * u2;

  return (nr/dr);
}

double e010(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e010a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e010b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y);
}

double hcrec(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = pow(expTauOverTheta(tau2, thtHC), -2);

  double y1 = (1 - 1./expRhoTau(rhoH, tau1)) * (1 - 1./expRhoTau(rhoM, tau1));
  double y2 = (1./expRhoTau(rhoH, tau1)) * (1 - 1./expRhoTau(rhoH, tau2)) * (1 - 1./expRhoTau(rhoM, tau1));
  double y3 = (1./expRhoTau(rhoH, tau1)) * (1 - 1./expRhoTau(rhoH, tau1)) * (1 - 1./expRhoTau(rhoM, tau2));

  double z1 = (1./expRhoTau(rhoM, tau2)) * pow(expTauOverTheta(tau2, thtHC), -2) * thtHC;
  double z2 = (expTauOverTheta(tau2, thtHC) * 1./expRhoTau(rhoH, tau2)) - 1.;

  double w1 = pow(expTauOverTheta(tau2, thtHC), -2)  * thtHC;
  double w2 = 1. - (expTauOverTheta(tau2, thtHC) * 1./(expRhoTau(rhoH, tau2)) * (1./expRhoTau(rhoM, tau2)));

  double u1 = (1./expRhoTau(rhoH, tau2)) * pow(expTauOverTheta(tau2, thtHC), -2) * thtHC;
  double u2 = (expTauOverTheta(tau2, thtHC) * 1./expRhoTau(rhoM, tau2)) - 1.;

  double d1 = thtHC * rhoH - 1.;
  double d2 = thtHC * (rhoH + rhoM) - 1.;

  double t1 = x1 * (y1 + y2 + y3);
  double t2 = rhoH * ((z1 * z2)/d1 + (w1 * w2)/d2);
  double t3 = rhoM * ((u1 * u2)/d1 + (w1 * w2)/d2);

  return (t1+t2+t3);
}

double e110x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoG, tau1 + tau2);
  double y = 3. + thtHCG * rhoG;
  return (x/(3.*y));
}

double e110a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hcrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e110x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double e110y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{   
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double y = thtHCG*rhoM;
  double z = 6. + thtHCG*(rhoG + rhoM);
    
  double u1 = 3. + thtHCG*rhoG;
  double u2 = 3. + thtHCG*(rhoG + rhoM);
  double u3 = 5. + thtHCG*(rhoG + rhoM);

  double nr = x1 * x2 * y * z;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double e110b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e110y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e110z(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoH, tau1 + tau2);

  double y = thtHCG*rhoH;
    
  double z1 = 3. + thtHCG * rhoG;
  double z2 = 5. + thtHCG * (rhoG + rhoH);

  double nr = x1 * x2 * y;
  double dr = 3. * z1 * z2;

  return (nr/dr);
}

double e110c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  double y = e110z(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e110d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z1 = thtHCG * rhoH;
  double z2 = thtHCG * rhoM;

  double w1 = 6. + thtHCG*rhoM; 
  double w2 = 13. + thtHCG*rhoM; 
  double w3 = 19. + thtHCG*rhoH + 3.*thtHCG*rhoM;

  double nr = x * y1 * y2 * y3 * z1 * z2 * (45. + thtHCG * (2. * thtHCG * (rhoG * rhoG) + rhoH * w1 + rhoM * w2 + rhoG * w3));
    
  double u1 = 3. + thtHCG*rhoG;
  double u2 = 5. + thtHCG*(rhoG + rhoH);
  double u3 = 3. + thtHCG*(rhoG + rhoM);
  double u4 = 5. + thtHCG*(rhoG + rhoM);
  double u5 = 3. + thtHCG*(rhoG + rhoH + rhoM);

  double dr = 3. * u1 * u2 * u3 * u4 * u5;
    
  return (nr/dr);
}

double e110(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e110a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e110b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = e110c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = e110d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+w);
}

double e101x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoM, tau1 + tau2);
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = 3. + thtHCG * rhoM;
  return ((x*y)/(3.*z));
}

double e101a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e101x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double e101y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double y = thtHCG*rhoG;
  double z = 6. + thtHCG*(rhoG + rhoM);
    
  double u1 = 3. + thtHCG*rhoM;
  double u2 = 3. + thtHCG*(rhoG + rhoM);
  double u3 = 5. + thtHCG*(rhoG + rhoM);

  double nr = x1 * x2 * y * z;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double e101b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e101y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double e101c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double z = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double w = thtHCG*rhoH;
  double u = 6. + thtHCG*(rhoH + rhoM);
    
  double v1 = 3. + thtHCG*rhoM;
  double v2 = 3. + thtHCG*(rhoH + rhoM);
  double v3 = 5. + thtHCG*(rhoH + rhoM);

  double nr = x * y1 * y2 * z * w * u;
  double dr = 3. * v1 * v2 * v3;

  return (nr/dr);
} 

double e101d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z1 = thtHCG*rhoH;
  double z2 = thtHCG*rhoM;
  double z3 = thtHCG*rhoG;

  double w1 = (z1 * z1) * (6. + z2);
  double w2 = (z3 * z3) * (6. + z1 + z2);
  double w3 = z1 * (63. + 28.*z2 + 3.*(z2*z2));
  double w4 = 2 * (90. + 63.*z2 + 14.*(z2*z2) + z2*z2*z2);
  double w5 = z3 * (63. + z1*z1 + 28. * z2 + 3.*(z2*z2) + (4. * z1) * (4. + z2));
    
  double nr = x * y1 * y2 * y3 * z1 * z3 * (w1 + w2 + w3 + w4 + w5);

  double u1 = 3. + z2;
  double u2 = 3. + z1 + z2 + z3;
  double u3 = 15. + z3*z3 + 8.*z2 + z2*z2 + 2.*z3*(4. + z2);
  double u4 = 15. + z3*z3 + 8.*z2 + z2*z2 + 2.*z1*(4. + z2);
    
  double dr = 3. * u1 * u2 * u3 * u4;

  return (nr/dr);
}

double e101(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e101a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e101b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = e101c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = e101d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+w);
}

double e011x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoH, tau1 + tau2);
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = 3. + thtHCG * rhoH;
  return ((x*y)/(3.*z));
}

double e011a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e011x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double e011y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoH, tau1 + tau2);
    
  double y = thtHCG*rhoG;
    
  double z1 = 3. + thtHCG*rhoH;
  double z2 = 5. + thtHCG*rhoG + thtHCG*rhoH;

  double nr = x1 * x2 * y;
  double dr = 3. * z1 * z2;

  return (nr/dr);
}

double e011b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e011y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e011c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double z = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double w = thtHCG*rhoM;
  double u = 6. + thtHCG*(rhoH + rhoM);
    
  double v1 = 3. + thtHCG*rhoH;
  double v2 = 3. + thtHCG*(rhoH + rhoM);
  double v3 = 5. + thtHCG*(rhoH + rhoM);

  double nr = x * y1 * y2 * z * w * u;
  double dr = 3. * v1 * v2 * v3;

  return (nr/dr);
} 

double e011d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2); 
  double y3 = 1./expRhoTau(rhoG, tau1 + tau2);

  double z1 = thtHCG*rhoH;
  double z2 = thtHCG*rhoM;
  double z3 = thtHCG*rhoG;

  double w = 45. + thtHCG * (19.*rhoH + 13.*rhoM + thtHCG*(rhoH+rhoM)*(2.*rhoH+rhoM) + rhoG*(6.+thtHCG*(rhoH+rhoM)));
    
  double nr = x * y1 * y2 * y3 * z2 * z3 * w;

  double u1 = 3. + z1;
  double u2 = 5. + thtHCG*(rhoH + rhoG);
  double u3 = 3. + thtHCG*(rhoH + rhoM);
  double u4 = 5. + thtHCG*(rhoH + rhoM);
  double u5 = 3. + thtHCG*(rhoH + rhoM + rhoG);
    
  double dr = 3. * u1 * u2 * u3 * u4 * u5;

  return (nr/dr);
}

double e011(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e011a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e011b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = e011c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = e011d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+w);
}

double e111x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoG, tau1 + tau2);
  double y = 3. + thtHCG * rhoG;
  return ((1./9.) * (1. - (3.*x)/y));
}

double e111a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hcrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double e111y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoM, tau1 + tau2);
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = thtHCG * rhoM;
  double w = 3. + thtHCG * rhoM;
  return ((x*y*z)/(9.*w));
}

double e111b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e111z(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoH, tau1 + tau2);
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = thtHCG * rhoH;
  double w = 3. + thtHCG * rhoH;
  return ((x*y*z)/(9.*w));
}

double e111c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111z(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double e111u(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoM, tau1 + tau2); 

  double y1 = thtHCG*rhoG;
  double y2 = thtHCG*rhoM;

  double z = 6. + thtHCG*(rhoM + rhoG);

  double nr = x1 * x2 * y1 * y2 * (z*z);

  double w1 = 3. + thtHCG*rhoG;
  double w2 = 3. + thtHCG*rhoM;
  double w3 = 3. + thtHCG*(rhoM + rhoG);
  double w4 = 5. + thtHCG*(rhoM + rhoG);

  double dr = 9. * w1 * w2 * w3 * w4;

  return (nr/dr);
}

double e111d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111u(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return(x*y);
}

double e111v(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{    
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoH, tau1 + tau2);

  double y1 = thtHCG*rhoG;
  double y2 = thtHCG*rhoH;

  double z = 6 + thtHCG*(rhoH + rhoG);

  double nr = x1 * x2 * y1 * y2 * z;

  double w1 = 3. + thtHCG*rhoG;
  double w2 = 3. + thtHCG*rhoH;
  double w3 = 5. + thtHCG*(rhoH + rhoG);

  double dr = 9. * w1 * w2 * w3;
  return (nr/dr);
}

double e111e(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111v(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double e111f(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z = 1. - 1./expRhoTau(rhoG, tau1 + tau2);

  double w1 = thtHCG*rhoH;
  double w2 = thtHCG*rhoM;

  double u = 6. + thtHCG*(rhoH + rhoM);

  double nr = x * y1 * y2 * z * w1 * w2 * (u*u);

  double v1 = 3. + thtHCG*rhoH;
  double v2 = 3. + thtHCG*rhoM;
  double v3 = 3. + thtHCG*(rhoH + rhoM);
  double v4 = 5. + thtHCG*(rhoH + rhoM);

  double dr = 9. * v1 * v2 * v3 * v4;

  return (nr/dr);
}

double e111g1(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. + thtHCG*(rhoM + rhoG);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoG;
  double z2 = 3. + thtHCG*(rhoM + rhoG);
  double z3 = 5. + thtHCG*(rhoM + rhoG);
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double e111g2(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. + thtHCG*(rhoM + rhoG);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoM;
  double z2 = 3. + thtHCG*(rhoM + rhoG);
  double z3 = 5. + thtHCG*(rhoM + rhoG);
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double e111g3(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{    
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. +thtHCG*(rhoH + rhoM);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoM;
  double z2 = 3. + thtHCG*(rhoH + rhoM);
  double z3 = 5. + thtHCG*(rhoH + rhoM); 
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double e111g4(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{  
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. +thtHCG*(rhoH + rhoM);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoH;
  double z2 = 3. + thtHCG*(rhoH + rhoM);
  double z3 = 5. + thtHCG*(rhoH + rhoM); 
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double e111g5(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{    
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double nr = x1 * x2 * x3;

  double z1 = 3. +thtHCG*rhoH;
  double z2 = 5. + thtHCG*(rhoH + rhoG);
  double z3 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3;

  return (nr/dr);
}

double e111g6(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double nr = x1 * x2 * x3; 

  double z1 = 3. +thtHCG*rhoG;
  double z2 = 5. + thtHCG*(rhoH + rhoG);
  double z3 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3;

  return (nr/dr);
}

double e111w(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2); 
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2); 
  double y3 = 1./expRhoTau(rhoG, tau1 + tau2);
   
  return ((1./3.) * x * y1 * y2 * y3);
}


double e111g(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double w = e111w(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);

  double x = e111g1(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  double y = e111g2(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = e111g3(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double t = e111g4(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double u = e111g5(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double v = e111g6(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (w*(x+y+z+t+u+v));
}

double e111(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = e111a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = e111c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double t = e111d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double u = e111e(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double v = e111f(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = e111g(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+t+u+v+w);
}

double IlsCoalHmmTransitionMatrix09::pcgtohg(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double dr = phg(tau2, thtHC);
  double x001 = e001(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x100 = e100(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x010 = e010(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x110 = e110(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x101 = e101(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x011 = e011(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x111 = e111(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  //cout << "e001\t" << x001 << endl;
  //cout << "e100\t" << x100 << endl;
  //cout << "e010\t" << x010 << endl;
  //cout << "e110\t" << x110 << endl;
  //cout << "e101\t" << x101 << endl;
  //cout << "e011\t" << x011 << endl;
  //cout << "e111\t" << x111 << endl;
  return ((x001+x100+x010+x110+x101+x011+x111)/dr);
}

double k001a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = (1./expTauOverTheta(tau2, thtHC)) * (1./expRhoTau(rhoH, tau1 + tau2)) * (1./expRhoTau(rhoM, tau1 + tau2));
  double y = 1. - (1./expRhoTau(rhoG, tau1 + tau2));
  double z = 5. + thtHCG*(rhoH + rhoM);
  return ((x*y)/(3.*z));
} 

double k001b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2); 
  double y3 = 1./expRhoTau(rhoG, tau1 + tau2);

  double z = thtHCG*rhoG;

  double w1 = 5. + thtHCG*(rhoH + rhoM);
  double w2 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double nr = x * y1 * y2 * y3 * z;
  double dr = 3. * w1 * w2;

  return (nr/dr);
} 

double k001(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k001a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k001b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y);
}

double k100x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = (1./expRhoTau(rhoG, tau1 + tau2)) * (1./expRhoTau(rhoM, tau1 + tau2));
  double y = 3. + thtHCG*(rhoM + rhoG);
  double z = 5 + thtHCG*(rhoM + rhoG);
  return (x/(y*z));
}

double k100a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k100x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k100b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z = thtHCG * rhoH;

  double u1 = 3. + thtHCG*(rhoG + rhoM);
  double u2 = 5. + thtHCG*(rhoG + rhoM);
  double u3 = 3. + thtHCG*(rhoG + rhoH + rhoM);
    
  double nr = x * y1 * y2 * y3 * z;
  double dr = u1 * u2 * u3;

  return (nr/dr);
}

double k100(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k100a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k100b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y);
} 

double k010x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = (1./expRhoTau(rhoH, tau1 + tau2)) * (1./expRhoTau(rhoG, tau1 + tau2));
  double y = 3./(2. * (3. + thtHCG*(rhoH + rhoG)));
  double z = 1./(2. * (5. + thtHCG*(rhoH + rhoG)));
  return ((1./3.) * x * (y - z));
}

double k010a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k010x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k010b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);
        
  double z = thtHCG * rhoM;
  double w = 6. + thtHCG*(rhoH + rhoG);

  double u1 = 3. + thtHCG*(rhoG + rhoH);
  double u2 = 5. + thtHCG*(rhoG + rhoH);
  double u3 = 3. + thtHCG*(rhoG + rhoH + rhoM);
    
  double nr = x * y1 * y2 * y3 * z * w;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double k010(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k010a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k010b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y);
}

double k110x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoG, tau1 + tau2);
  double y = 3. + thtHCG * rhoG;
  return (x/(3.*y));
}

double k110a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hcrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k110x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k110y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double y = thtHCG*rhoM;
  double z = 6. + thtHCG*(rhoG + rhoM);
    
  double u1 = 3. + thtHCG*rhoG;
  double u2 = 3. + thtHCG*(rhoG + rhoM);
  double u3 = 5. + thtHCG*(rhoG + rhoM);

  double nr = x1 * x2 * y * z;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double k110b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k110y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double k110z(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoG, tau1 + tau2); 
    
  double y = thtHCG*rhoH;
  double z = 6. + thtHCG*(rhoG + rhoH);
    
  double u1 = 3. + thtHCG*rhoG;
  double u2 = 3. + thtHCG*(rhoG + rhoH); 
  double u3 = 5. + thtHCG*(rhoG + rhoH);

  double nr = x1 * x2 * y * z;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double k110c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k110z(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k110d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2); 
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z1 = thtHCG*rhoH;
  double z2 = thtHCG*rhoM;
  double z3 = thtHCG*rhoG;

  double w1 = 2. * pow(z3, 3.);
  double w2 = (z1*z1) * (6. + z2);
  double w3 = z1 * (7 + z2) * (9. + z2);
  double w4 = (3. * z2) * (21. + 2. * z2);
  double w5 = (z3*z3) * (28. + 3.*z1 + 3.*z2);
  double w6 = z3 * (126. + z1*z1 + 4.*z1*(7+z2) + z2*(28.+z2));

  double nr = x * y1 * y2 * y3 * z1 * z2 * (180. + w1 + w2 + w3 + w4 + w5 + w6);

  double u1 = 3. + z3;
  double u2 = 3. + z1 + z3;
  double u3 = 5. + z1 + z3;
  double u4 = 3. + z2 + z3;
  double u5 = 5. + z2 + z3;
  double u6 = 3. + z1 + z2 + z3;

  double dr = 3. * u1 * u2 * u3 * u4 * u5 * u6;

  return (nr/dr);
}

double k110(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k110a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k110b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = k110c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = k110d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+w);
}

double k101x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoM, tau1 + tau2); 
  double y = 1 - 1./expRhoTau(rhoG, tau1 + tau2); 
  double z = 3. + thtHCG * rhoM;
  return ((x*y)/(3.*z));
}

double k101a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k101x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k101y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double y = thtHCG*rhoG;
  double z = 6. + thtHCG*(rhoG + rhoM);
    
  double u1 = 3. + thtHCG*rhoM;
  double u2 = 3. + thtHCG*(rhoG + rhoM);
  double u3 = 5. + thtHCG*(rhoG + rhoM);

  double nr = x1 * x2 * y * z;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double k101b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k101y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k101c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);     

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2); 
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);
    
  double z = 1 - 1./expRhoTau(rhoG, tau1 + tau2);
  double w = thtHCG*rhoH;
    
  double u1 = 3. + thtHCG*rhoM;
  double u2 = 5. + thtHCG*(rhoH + rhoM);

  double nr = x * y1 * y2 * z * w;
  double dr = 3. * u1 * u2;

  return (nr/dr);
}

double k101d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoG, tau1 + tau2); 
  double y2 = 1./expRhoTau(rhoH, tau1 + tau2); 
  double y3 = 1./expRhoTau(rhoM, tau1 + tau2); 

  double z1 = thtHCG*rhoH;
  double z2 = thtHCG*rhoM;
  double z3 = thtHCG*rhoG;

  double w1 = z3*z3;
  double w2 = z1 * (6. + z2);
  double w3 = z2 * (19. + 2.*z2);
  double w4 = z3 * (13. + z1 + 3.*z2);

  double nr = x * y1 * y2 * y3 * z1 * z3 * (45. + w1 + w2 + w3 + w4);

  double u1 = 3. + z2;
  double u2 = 3. + z2 + z3;
  double u3 = 5. + z2 + z3;
  double u4 = 5. + z1 + z2;
  double u5 = 3. + z1 + z2 + z3;

  double dr = 3. * u1 * u2 * u3 * u4 * u5;

  return (nr/dr);
}

double k101(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k101a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k101b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = k101c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = k101d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+w);
}

double k011x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoH, tau1 + tau2);
  double y = 1 - 1./expRhoTau(rhoG, tau1 + tau2); 
  double z = 3. + thtHCG * rhoH;
  return ((x*y)/(3.*z));
}

double k011a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k011x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double k011y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoH, tau1 + tau2);
    
  double y = thtHCG*rhoG;
  double z = 6. + thtHCG*(rhoG + rhoH);
    
  double u1 = 3. + thtHCG*rhoH;
  double u2 = 3. + thtHCG*(rhoG + rhoH);
  double u3 = 5. + thtHCG*(rhoG + rhoH);

  double nr = x1 * x2 * y * z;
  double dr = 3. * u1 * u2 * u3;

  return (nr/dr);
}

double k011b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k011y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return( x*y);
}

double k011c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);
   
  double z = 1 - 1./expRhoTau(rhoG, tau1 + tau2);
  double w = thtHCG*rhoM;
    
  double u1 = 3. + thtHCG*rhoH;
  double u2 = 5. + thtHCG*(rhoH + rhoM);

  double nr = x * y1 * y2 * z * w;
  double dr = 3. * u1 * u2;

  return (nr/dr);
} 

double k011d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoG, tau1 + tau2);

  double z1 = thtHCG*rhoH;
  double z2 = thtHCG*rhoM;
  double z3 = thtHCG*rhoG;

  double w1 = pow(z3, 3.);
  double w2 = 2.*(z3*z3)*(9.+2.*z1+z2);
  double w3 = z3 * (5. * pow(5.+z1,2.) + z2*(27.+5.*z1) + z2*z2);
  double w4 = (54. + 2.*z1*(11.+z1) + z2*(6.+z1)) * (5. + z1 + z2);
    
  double nr = x * y1 * y2 * y3 * z2 * z3 * (w1 + w2 + w3 + w4);

  double u1 = 3. + z1;
  double u2 = 3. + z1 + z3;
  double u3 = 5. + z1 + z3;
  double u4 = 5. + z1 + z2;
  double u5 = 3. + z1 + z2 + z3;
  double u6 = 5. + z1 + z2 + z3;
    
  double dr = 3. * u1 * u2 * u3 * u4 * u5 * u6;

  return (nr/dr);
}

double k011(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k011a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k011b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = k011c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = k011d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+w);
}

double k111x(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoG, tau1 + tau2);
  double y = 3. + thtHCG * rhoG;
  return ((1./9.) * (1. - (3.*x)/y));
}

double k111a(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hcrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k111x(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k111y(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoM, tau1 + tau2);
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = thtHCG * rhoM;
  double w = 3. + thtHCG * rhoM;
  return ((x*y*z)/(9.*w));
}

double k111b(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k111y(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k111z(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expRhoTau(rhoH, tau1 + tau2);
  double y = 1. - 1./expRhoTau(rhoG, tau1 + tau2);
  double z = thtHCG * rhoH;
  double w = 3. + thtHCG * rhoH;
  return ((x*y*z)/(9.*w));
}

double k111c(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = e111z(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double k111u(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{

  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoM, tau1 + tau2);

  double y1 = thtHCG*rhoG;
  double y2 = thtHCG*rhoM;

  double z = 6. + thtHCG*(rhoM + rhoG);

  double nr = x1 * x2 * y1 * y2 * (z*z);

  double w1 = 3. + thtHCG*rhoG;
  double w2 = 3. + thtHCG*rhoM;
  double w3 = 3. + thtHCG*(rhoM + rhoG);
  double w4 = 5. + thtHCG*(rhoM + rhoG);

  double dr = 9. * w1 * w2 * w3 * w4;

  return (nr/dr);
}

double k111d(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = hrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k111u(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG); 
  return (x*y);
}

double k111v(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{ 
  double x1 = 1./expRhoTau(rhoG, tau1 + tau2);
  double x2 = 1./expRhoTau(rhoH, tau1 + tau2);

  double y1 = thtHCG*rhoG;
  double y2 = thtHCG*rhoH;

  double z = 6 + thtHCG*(rhoH + rhoG);

  double nr = x1 * x2 * y1 * y2 * (z*z);

  double w1 = 3. + thtHCG*rhoG;
  double w2 = 3. + thtHCG*rhoH;
  double w3 = 3. + thtHCG*(rhoH + rhoG);
  double w4 = 5. + thtHCG*(rhoH + rhoG);

  double dr = 9. * w1 * w2 * w3 * w4;
  return (nr/dr);
}

double k111e(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = mrec(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k111v(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x*y);
}

double k111f(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{   
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);

  double z = 1. - 1./expRhoTau(rhoG, tau1 + tau2);

  double w1 = thtHCG*rhoH;
  double w2 = thtHCG*rhoM;

  double u = 6. + thtHCG*(rhoH + rhoM);

  double nr = x * y1 * y2 * z * w1 * w2 * u;

  double v1 = 3. + thtHCG*rhoH;
  double v2 = 3. + thtHCG*rhoM;
  double v3 = 5. + thtHCG*(rhoH + rhoM);

  double dr = 9. * v1 * v2 * v3;

  return (nr/dr);
}

double k111g1(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{   
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. + thtHCG*(rhoM + rhoG);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoG;
  double z2 = 3. + thtHCG*(rhoM + rhoG);
  double z3 = 5. + thtHCG*(rhoM + rhoG);
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double k111g2(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. + thtHCG*(rhoM + rhoG);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoM;
  double z2 = 3. + thtHCG*(rhoM + rhoG);
  double z3 = 5. + thtHCG*(rhoM + rhoG);
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double k111g3(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
    
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. +thtHCG*(rhoH + rhoG);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoG;
  double z2 = 3. + thtHCG*(rhoH + rhoG);
  double z3 = 5. + thtHCG*(rhoH + rhoG);
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return (nr/dr);
}

double k111g4(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
    
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double y = 6. +thtHCG*(rhoH + rhoG);
    
  double nr = x1 * x2 * x3 * y;

  double z1 = 3. +thtHCG*rhoH;
  double z2 = 3. + thtHCG*(rhoH + rhoG);
  double z3 = 5. + thtHCG*(rhoH + rhoG);
  double z4 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3 * z4;

  return(nr/dr);
}

double k111g5(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{    
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double nr = x1 * x2 * x3;

  double z1 = 3. +thtHCG*rhoH;
  double z2 = 5. + thtHCG*(rhoH + rhoM);
  double z3 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3;

  return (nr/dr);
}

double k111g6(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{  
  double x1 = thtHCG*rhoH;
  double x2 = thtHCG*rhoM;
  double x3 = thtHCG*rhoG;

  double nr = x1 * x2 * x3;

  double z1 = 3. +thtHCG*rhoM;
  double z2 = 5. + thtHCG*(rhoH + rhoM);
  double z3 = 3. + thtHCG*(rhoH + rhoM + rhoG);

  double dr = 3. * z1 * z2 * z3;

  return (nr/dr);
}

double k111w(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = 1./expTauOverTheta(tau2, thtHC);

  double y1 = 1./expRhoTau(rhoH, tau1 + tau2);
  double y2 = 1./expRhoTau(rhoM, tau1 + tau2);
  double y3 = 1./expRhoTau(rhoG, tau1 + tau2);
   
  return ((1./3.) * x * y1 * y2 * y3);
}

double k111g(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double w = k111w(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);

  double x = k111g1(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k111g2(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = k111g3(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double t = k111g4(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double u = k111g5(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double v = k111g6(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (w*(x+y+z+t+u+v));
}

double k111(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double x = k111a(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double y = k111b(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double z = k111c(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double t = k111d(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double u = k111e(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double v = k111f(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double w = k111g(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  return (x+y+z+t+u+v+w);
}

double IlsCoalHmmTransitionMatrix09::pcgtohc2(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG)
{
  double dr = phc2(tau2, thtHC);

  double x001 = k001(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x100 = k100(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x010 = k010(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x110 = k110(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x101 = k101(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x011 = k011(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);
  double x111 = k111(rhoH, rhoM, rhoG, tau1, tau2, thtHC, thtHCG);

  //cerr << x001 << "+" << x100 << "+" << x010 << "+" << x110 << "+" << x101 << "+" << x011 << "+" << x111 << "/" << dr << endl;
  return ((x001+x100+x010+x110+x101+x011+x111)/dr);
}

void IlsCoalHmmTransitionMatrix09::actualizeTransitionMatrix()
{
  double tau1, tau2, theta1, theta2;
  tau1 = dynamic_cast<const ThreeSpeciesCoalHmmStateAlphabet *>(getHmmStateAlphabet())->getTau1(); 
  tau2 = dynamic_cast<const ThreeSpeciesCoalHmmStateAlphabet *>(getHmmStateAlphabet())->getTau2(); 
  theta1 = dynamic_cast<const ThreeSpeciesCoalHmmStateAlphabet *>(getHmmStateAlphabet())->getTheta1(); 
  theta2 = dynamic_cast<const ThreeSpeciesCoalHmmStateAlphabet *>(getHmmStateAlphabet())->getTheta2();
  if (rhoOption_ == ONE_RHO)
  {
    rho1_ = getParameter_(0).getValue();
    rho2_ = getParameter_(0).getValue();
    rho3_ = getParameter_(0).getValue();
  }
  else if (rhoOption_ == TWO_RHOS_12)
  {
    rho1_ = getParameter_(0).getValue();
    rho2_ = getParameter_(0).getValue();
    rho3_ = getParameter_(1).getValue();
  }
  else
  {
    rho1_ = getParameter_(0).getValue();
    rho2_ = getParameter_(1).getValue();
    rho3_ = getParameter_(2).getValue();
  }
  s_  = threeS(rho1_, rho2_, rho3_, tau1, tau2, theta1)/3.;
  u_  = u(rho1_, rho2_, rho3_, tau1, tau2, theta1);
  v1_ = pcgtohc2(rho1_, rho2_, rho3_, tau1, tau2, theta1, theta2);
  v2_ = pcgtohg(rho1_, rho2_, rho3_, tau1, tau2, theta1, theta2);


  if (isnan(s_) || s_ < 0.)
  {
    cerr << "DEBUG: s = " << s_ << endl;
    cerr << "rho1 = " << rho1_ << endl;
    cerr << "rho2 = " << rho2_ << endl;
    cerr << "rho3 = " << rho3_ << endl;
    cerr << "tau1  = " << tau1 << endl;
    cerr << "tau2  = " << tau2 << endl;
    cerr << "theta1  = " << theta1 << endl;
    cerr << "theta2  = " << theta2 << endl;
    s_ = 0.;
  }
  if (isnan(u_) || u_ < 0.)
  {
    cerr << "DEBUG: u = " << u_ << endl;
    cerr << "rho1 = " << rho1_ << endl;
    cerr << "rho2 = " << rho2_ << endl;
    cerr << "rho3 = " << rho3_ << endl;
    cerr << "tau1  = " << tau1 << endl;
    cerr << "tau2  = " << tau2 << endl;
    cerr << "theta1  = " << theta1 << endl;
    cerr << "theta2  = " << theta2 << endl;
    u_ = 0.;
  }
  if (isnan(v1_) || v1_ < 0.)
  {
    cerr << "DEBUG: v1 = " << v1_ << endl;
    cerr << "rho1 = " << rho1_ << endl;
    cerr << "rho2 = " << rho2_ << endl;
    cerr << "rho3 = " << rho3_ << endl;
    cerr << "tau1  = " << tau1 << endl;
    cerr << "tau2  = " << tau2 << endl;
    cerr << "theta1  = " << theta1 << endl;
    cerr << "theta2  = " << theta2 << endl;
    v1_ = 0.;
  }
  if (isnan(v2_) || v2_ < 0.)
  {
    cerr << "DEBUG: v2 = " << v2_ << endl;
    cerr << "rho1    = " << rho1_ << endl;
    cerr << "rho2    = " << rho2_ << endl;
    cerr << "rho3    = " << rho3_ << endl;
    cerr << "tau1    = " << tau1 << endl;
    cerr << "tau2    = " << tau2 << endl;
    cerr << "theta1  = " << theta1 << endl;
    cerr << "theta2  = " << theta2 << endl;
    v2_ = 0.;
  }
  //cout << t1 << "\t" << t2 << "\t" << theta1 << "\t" << theta2 << "\t" << _rhoH << "\t" << _rhoM << "\t" << _rhoG << endl;
  //cout << s_ << "\t" << u_ << "\t" << v1_ << "\t" << v2_ << endl;
  //cout << e001a(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e001b(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e001(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << hrec(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << mrec(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e100x(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e100a(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e100b(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e100(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e010x(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e010a(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e010b(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << e010(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;
  //cout << hcrec(_rhoH, _rhoM, _rhoG, t1, t2, theta1, theta2) << endl;

  unsigned int n = getNumberOfStates();
  transitions_(0, 0) = 1. - static_cast<double>(n - 1) * s_;
  for(unsigned int i = 1; i < n; i++)
    transitions_(0, i) = s_;

  for(unsigned int i = 1; i < n; i++)
    transitions_(i, 0) = u_;

  transitions_(1, 1) =  1. - u_ - 2. * v1_;
  transitions_(2, 2) =  1. - u_ - v1_ - v2_;
  transitions_(3, 3) =  1. - u_ - v1_ - v2_;

  transitions_(1, 2) = v1_;
  transitions_(1, 3) = v1_;
  
  transitions_(2, 1) = v1_;
  transitions_(2, 3) = v2_;

  transitions_(3, 1) = v1_;
  transitions_(3, 2) = v2_;
  //MatrixTools::print(transitions_);


  double phi = 1. / (1. + 3. * s_ / u_);
  freqs_[0] = phi;
  for (unsigned int i = 1; i < n; i++)
    freqs_[i] = (1. - phi) / static_cast<double>(n - 1);
}

void IlsCoalHmmTransitionMatrix09::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  if (rhoOption_ == ONE_RHO)
  {
    (out << "rho = " << rho1_).endLine();
  }
  else if (rhoOption_ == TWO_RHOS_12)
  {
    (out << "rho12 = " << rho1_).endLine();
    (out << "rho3  = " << rho3_).endLine();
  }
  else
  {
    (out << "rho1 = " << rho1_).endLine();
    (out << "rho2 = " << rho2_).endLine();
    (out << "rho3 = " << rho3_).endLine();
  }
  (out << "s = " << s_).endLine();
  (out << "u = " << u_).endLine();
  (out << "v1 = " << v1_).endLine();
  (out << "v2 = " << v2_).endLine();
}

