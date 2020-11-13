//
// File: IlsCoalHmmTransitionMatrix07.cpp
// Created by: Julien Dutheil
// Created on: Mon Oct 29 09:08 2007
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

#include "IlsCoalHmmTransitionMatrix07.h"

IlsCoalHmmTransitionMatrix07::IlsCoalHmmTransitionMatrix07(const CoalHmmStateAlphabet* hiddenAlphabet, double s, double u, double v1, double v2, bool onev):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet), sigma_(), mu_(), nu_(), lambda_(), onev_(onev)
{
  if(onev)
  {
    v2 = v1;
    lambda_ = 0.5;
  }
  else
  {
    lambda_ = v1 / (v1 + v2);
  }
  sigma_  = 3. * s;
  mu_     = u + v1 + v2;
  nu_     = (v1 + v2) / mu_;
  addParameter_(new bpp::Parameter("sigma" , sigma_, bpp::Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new bpp::Parameter("nu"    , nu_   , bpp::Parameter::PROP_CONSTRAINT_EX));
  addParameter_(new bpp::Parameter("mu"    , mu_   , bpp::Parameter::PROP_CONSTRAINT_EX));
  if (!onev) 
    addParameter_(new bpp::Parameter("lambda", lambda_, bpp::Parameter::PROP_CONSTRAINT_EX));
  transitions_.resize(4, 4);
  freqs_.resize(4);
  actualizeTransitionMatrix();
}

void IlsCoalHmmTransitionMatrix07::actualizeTransitionMatrix()
{
  transitions_(0, 0) = 1. - sigma_;
  for(unsigned int i = 1; i < 4; i++)
    transitions_(0, i) = sigma_ / (3.);

  for(unsigned int i = 1; i < 4; i++)
    transitions_(i, 0) = (1. - nu_) * mu_;

  transitions_(1, 1) = 1. - mu_ * (1. - nu_ * (1. - 2. * lambda_));
  transitions_(2, 2) = 1. - mu_;
  transitions_(3, 3) = 1. - mu_;

  transitions_(1, 2) = nu_ * mu_ * lambda_;
  transitions_(1, 3) = nu_ * mu_ * lambda_;
  transitions_(2, 1) = nu_ * mu_ * lambda_;
  transitions_(3, 1) = nu_ * mu_ * lambda_;
  
  transitions_(2, 3) = nu_ * mu_ * (1. - lambda_);
  transitions_(3, 2) = nu_ * mu_ * (1. - lambda_);

//  for(unsigned int i = 1; i < 4; i++)
//  {
//    _transitions(i, i) = 1. - mu_;
//    for(unsigned int j = 1; j < 4; j++)
//      if(j != i) _transitions(i, j) = nu_ * mu_ / ((double)_nbStates - 2.);
//  }

  double phi = 1. / (1. + sigma_ / ((1. - nu_) * mu_));
  freqs_[0] = phi;
  for(unsigned int i = 1; i < 4; i++)
    freqs_[i] = (1. - phi) / 3.;
}

void IlsCoalHmmTransitionMatrix07::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  double s  = sigma_ / ((double)getNumberOfStates() - 1.);
  double u  = (1. - nu_) * mu_;
  (out << "s  = " << s).endLine();
  (out << "u  = " << u).endLine();
  if (onev_)
  {
    double v = nu_ * mu_ / 2.;
    (out << "v = " << v).endLine();
  }
  else
  {
    double v1 = nu_ * mu_ * lambda_;
    double v2 = nu_ * mu_ * (1. - lambda_);
    (out << "v1 = " << v1).endLine();
    (out << "v2 = " << v2).endLine();
  }
}

