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

NonClockAverageCoalHmmStateAlphabet::NonClockAverageCoalHmmStateAlphabet(
		AverageCoalHmmStateAlphabet* alphabet,
	       	const string& species, double error):
  AbstractParametrizable(""),
  alphabet_(alphabet),
  speciesWithError1_(species),
  speciesWithError2_(""),
  error1_(error),
  error2_(0),
  hasSecondError_(false),
  brlenParameters_()
{
  trees_.resize(alphabet_->getNumberOfStates());
  brlenParameters_.resize(alphabet_->getNumberOfStates());
  for (size_t i = 0; i < alphabet_->getNumberOfStates(); ++i) {
    trees_[i] = alphabet_->getState(i).clone();
    brlenParameters_[i] = alphabet_->getBranchLengthParametersForState(i);
  }
  species_ = alphabet->getSpeciesNames();

  addParameters_(alphabet_->getParameters());
  addParameter_(new Parameter(speciesWithError1_ + "_error", error1_, Parameter::R_PLUS));
  fireParameterChanged(getParameters());
}

NonClockAverageCoalHmmStateAlphabet::NonClockAverageCoalHmmStateAlphabet(
		AverageCoalHmmStateAlphabet* alphabet,
	       	const string& species1,
	       	const string& species2,
	       	double error1,
	       	double error2):
  AbstractParametrizable(""),
  alphabet_(alphabet),
  speciesWithError1_(species1),
  speciesWithError2_(species2),
  error1_(error1),
  error2_(error2),
  hasSecondError_(true),
  brlenParameters_()
{
  trees_.resize(alphabet_->getNumberOfStates());
  brlenParameters_.resize(alphabet_->getNumberOfStates());
  for (size_t i = 0; i < alphabet_->getNumberOfStates(); ++i) {
    trees_[i] = alphabet_->getState(i).clone();
    brlenParameters_[i] = alphabet_->getBranchLengthParametersForState(i);
  }
  species_ = alphabet->getSpeciesNames();

  addParameters_(alphabet_->getParameters());
  addParameter_(new Parameter(speciesWithError1_ + "_error", error1_, Parameter::R_PLUS));
  addParameter_(new Parameter(speciesWithError2_ + "_error", error2_, Parameter::R_PLUS));
  fireParameterChanged(getParameters());
}

void NonClockAverageCoalHmmStateAlphabet::fireParameterChanged(const ParameterList& params)
{
  alphabet_->matchParametersValues(params);
  error1_ = getParameterValue(speciesWithError1_ + "_error");
  if (hasSecondError_){
    error2_ = getParameterValue(speciesWithError2_ + "_error");
  }

  //Update topologies and branch lengths:
  for (size_t i = 0; i < alphabet_->getNumberOfStates(); ++i) {
    brlenParameters_[i].matchParametersValues(alphabet_->getBranchLengthParametersForState(i));
    *trees_[i] = alphabet_->getState(i);
    Node* node = trees_[i]->getNode(speciesWithError1_);
    node->setDistanceToFather(node->getDistanceToFather() + error1_);
    brlenParameters_[i].setParameterValue(
        "BrLen" + TextTools::toString(node->getId()), node->getDistanceToFather()); 
    if (hasSecondError_) {
      Node* node2 = trees_[i]->getNode(speciesWithError2_);
      node2->setDistanceToFather(node2->getDistanceToFather() + error2_);
      brlenParameters_[i].setParameterValue(
          "BrLen" + TextTools::toString(node2->getId()), node2->getDistanceToFather()); 
    }
  }
}

void NonClockAverageCoalHmmStateAlphabet::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  alphabet_->printUserFriendlyParameters(out);
  (out << speciesWithError1_ << "_error = " << error1_).endLine();
  if (hasSecondError_) {
    (out << speciesWithError2_ << "_error = " << error2_).endLine();
  }
}


