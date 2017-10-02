//
// File: NonClockAverageCoalHmmStateAlphabet.cpp
// Created by: Julien Dutheil
// Created on: Tue May 01 11:31 2010
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

#include "NonClockAverageCoalHmmStateAlphabet.h"

using namespace bpp;
using namespace std;

NonClockAverageCoalHmmStateAlphabet::NonClockAverageCoalHmmStateAlphabet(AverageCoalHmmStateAlphabet* alphabet, const string& species, double error):
  AbstractParametrizable(""),
  alphabet_(alphabet),
  speciesWithError_(species),
  error_(error),
  brlenParameters_()
{
  trees_.resize(alphabet_->getNumberOfStates());
  brlenParameters_.resize(alphabet_->getNumberOfStates());
  for (unsigned int i = 0; i < alphabet_->getNumberOfStates(); ++i) {
    trees_[i] = alphabet_->getState(i).clone();
    brlenParameters_[i] = alphabet_->getBranchLengthParametersForState(i);
  }
  species_ = alphabet->getSpeciesNames();

  addParameters_(alphabet_->getParameters());
  addParameter_(new Parameter(speciesWithError_ + "_error", error_, &Parameter::R_PLUS, false));
  fireParameterChanged(getParameters());
}

void NonClockAverageCoalHmmStateAlphabet::fireParameterChanged(const ParameterList& params)
{
  alphabet_->matchParametersValues(params);
  error_ = getParameterValue(speciesWithError_ + "_error");

  //Update topologies and branch lengths:
  for (unsigned int i = 0; i < alphabet_->getNumberOfStates(); ++i) {
    *trees_[i] = alphabet_->getState(i);
    Node* node = trees_[i]->getNode(speciesWithError_);
    node->setDistanceToFather(node->getDistanceToFather() + error_);
    brlenParameters_[i].matchParametersValues(alphabet_->getBranchLengthParametersForState(i));
    brlenParameters_[i].setParameterValue(
        "BrLen" + TextTools::toString(node->getId()), node->getDistanceToFather()); 
  }
}

void NonClockAverageCoalHmmStateAlphabet::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  alphabet_->printUserFriendlyParameters(out);
  (out << speciesWithError_ << "_error = " << error_).endLine();
}


