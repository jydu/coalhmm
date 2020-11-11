//
// File: RateAndCoalHmmTransitionMatrix.cpp
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

#include "RateAndCoalHmmTransitionMatrix.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;
using namespace bpp;

RateAndCoalHmmTransitionMatrix::RateAndCoalHmmTransitionMatrix(
  const RateAndCoalHmmStateAlphabet* hiddenAlphabet,
  CoalHmmTransitionMatrix* coalModel,
  double w) :
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet),
  coalModel_(coalModel), omega_(),
  nbCoalStates_(coalModel->getNumberOfStates()),
  nbRates_(hiddenAlphabet->getRateDistribution().getNumberOfCategories())
{
  if (nbRates_ == 1)
    throw Exception("RateAndCoalHmmTransitionMatrix: At least two rate categories are needed!");
  omega_        = w * (static_cast<double>(nbRates_ - 1));
  addParameter_(new Parameter("omega", omega_, Parameter::PROP_CONSTRAINT_EX));
  addParameters_(coalModel_->getParameters());
  size_t n = getNumberOfStates();
  transitions_.resize(n, n);
  freqs_.resize(n);
  actualizeTransitionMatrix();
}

void RateAndCoalHmmTransitionMatrix::actualizeTransitionMatrix()
{
  for (size_t ir = 0; ir < nbRates_; ir++)
  {
    for (size_t is = 0; is < nbCoalStates_; is++)
    {
      for (size_t jr = 0; jr < nbRates_; jr++)
      {
        for (size_t js = 0; js < nbCoalStates_; js++)
        {
          if (ir == jr)
          // Same rate class.
          // The corresponding submatrix is equal to (1-omega)*the coal model matrix.
          {
            transitions_(ir * nbCoalStates_ + is, jr * nbCoalStates_ + js) =
              coalModel_->Pij(is, js) * (1. - omega_);
          }
          else if (is == js)
          // Not the same rate, but the same state.
          // The corresponding probability is then omega / (n - 1),
          // n being the number of rates.
          {
            transitions_(ir * nbCoalStates_ + is, jr * nbCoalStates_ + js) = omega_ / static_cast<double>(nbRates_ - 1);
          }
          else
          // Not the same rate, not the same state: probability of change is 0.
          {
            transitions_(ir * nbCoalStates_ + is, jr * nbCoalStates_ + js) = 0.;
          }
        }
      }
      // Equilibrium frequency:
      freqs_[ir * nbCoalStates_ + is] = coalModel_->getEquilibriumFrequencies()[is] / static_cast<double>(nbRates_);
    }
  }
  // MatrixTools::print(transitions_);
}

void RateAndCoalHmmTransitionMatrix::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  coalModel_->printUserFriendlyParameters(out);
  double w = omega_ / static_cast<double>(nbRates_ - 1);
  cout << "omega = " << w << endl;
}

