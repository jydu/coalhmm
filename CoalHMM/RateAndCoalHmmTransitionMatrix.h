//
// File: RateAndCoalHmmTransitionMatrix.h
// Created by: Julien Dutheil
// Created on: Wed Nov 21 11:23 2007
//

// This file is part of the CoalHMM program and library.
//
// CoalHMM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.
//
// CoalHMM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CoalHMM.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _RATEANDCOALHMMTRANSITIONMATRIX_H_
#define _RATEANDCOALHMMTRANSITIONMATRIX_H_

#include "CoalHmmTransitionMatrix.h"
#include "RateAndCoalHmmStateAlphabet.h"

class RateAndCoalHmmTransitionMatrix :
  public AbstractCoalHmmTransitionMatrix
{
protected:
  CoalHmmTransitionMatrix* coalModel_;
  double omega_;
  unsigned int nbCoalStates_;
  unsigned int nbRates_;

public:
  RateAndCoalHmmTransitionMatrix(
    const RateAndCoalHmmStateAlphabet* hiddenAlphabet,
    CoalHmmTransitionMatrix* coalModel,
    double w) throw (bpp::Exception);

  RateAndCoalHmmTransitionMatrix(const RateAndCoalHmmTransitionMatrix& tm):
    AbstractCoalHmmTransitionMatrix(tm),
    coalModel_(tm.coalModel_),
    omega_(tm.omega_),
    nbCoalStates_(tm.nbCoalStates_),
    nbRates_(tm.nbRates_)
  {}

  RateAndCoalHmmTransitionMatrix& operator=(const RateAndCoalHmmTransitionMatrix& tm)
  {
    AbstractCoalHmmTransitionMatrix::operator=(tm);
    coalModel_ = tm.coalModel_;
    omega_ = tm.omega_;
    nbCoalStates_ = tm.nbCoalStates_;
    nbRates_ = tm.nbRates_;
    return *this;
  }

  virtual ~RateAndCoalHmmTransitionMatrix()
  {}

  RateAndCoalHmmTransitionMatrix* clone() const { return new RateAndCoalHmmTransitionMatrix(*this); }

public:
  void fireParameterChanged(const bpp::ParameterList& pl)
  {
    coalModel_->matchParametersValues(getParameters());
    omega_ = getParameter_(0).getValue();
    actualizeTransitionMatrix();
  }

  void printUserFriendlyParameters(bpp::OutputStream& out) const;

protected:
  void actualizeTransitionMatrix();
};

#endif // _RATEANDCOALHMMTRANSITIONMATRIX_H_

