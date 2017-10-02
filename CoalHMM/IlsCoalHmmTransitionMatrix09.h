//
// File: BasicCoalHmmTransitionMatrix09.h
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

#ifndef _ILSCOALHMMTRANSITIONMATRIX09_H_
#define _ILSCOALHMMTRANSITIONMATRIX09_H_

#include "CoalHmmTransitionMatrix.h"
#include "ThreeSpeciesCoalHmmStateAlphabet.h"

class IlsCoalHmmTransitionMatrix09:
  public AbstractCoalHmmTransitionMatrix
{
  private:
    short rhoOption_;
    double rho1_, rho2_, rho3_, s_, u_, v1_, v2_;
    
  public:
    static const short ONE_RHO;
    static const short TWO_RHOS_12;
    static const short THREE_RHOS;
    
  public:
    IlsCoalHmmTransitionMatrix09(const ThreeSpeciesCoalHmmStateAlphabet* hiddenAlphabet, double rho1, double rho2, double rho3, double lowerRho, double upperRho);
    IlsCoalHmmTransitionMatrix09(const ThreeSpeciesCoalHmmStateAlphabet* hiddenAlphabet, double rho12, double rho3, short rhoOption, double lowerRho, double upperRho);
    IlsCoalHmmTransitionMatrix09(const ThreeSpeciesCoalHmmStateAlphabet* hiddenAlphabet, double rho, double lowerRho, double upperRho);
    
    IlsCoalHmmTransitionMatrix09* clone() const { return new IlsCoalHmmTransitionMatrix09(*this); }
    
  public:
    void fireParameterChanged(const bpp::ParameterList& pl)
    {
      actualizeTransitionMatrix();
    }
    
    void printUserFriendlyParameters(bpp::OutputStream& out) const;
    
  protected:
    void actualizeTransitionMatrix();
    
  public:
    
    static double threeS(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC);
    static double threeS(double rho, double tau1, double tau2, double thtHC)
    {
        return threeS(rho, rho, rho, tau1, tau2, thtHC);
    }
    static double u(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC);
    static double u(double rho, double tau1, double tau2, double thtHC)
    {
        return u(rho, rho, rho, tau1, tau2, thtHC);
    }
    static double pcgtohg(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG);
    static double pcgtohg(double rho, double tau1, double tau2, double thtHC, double thtHCG)
    {
        return pcgtohg(rho, rho, rho, tau1, tau2, thtHC, thtHCG);
    }
    static double pcgtohc2(double rhoH, double rhoM, double rhoG, double tau1, double tau2, double thtHC, double thtHCG);
    static double pcgtohc2(double rho, double tau1, double tau2, double thtHC, double thtHCG)
    {
        return pcgtohc2(rho, rho, rho, tau1, tau2, thtHC, thtHCG);
    }
    
};

#endif //_ILSCOALHMMTRANSITIONMATRIX09_H_

