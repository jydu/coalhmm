//
// File: RateAndCoalHmmStateAlphabet.cpp
// Created by: Julien Dutheil
// Created on: Wed Nov 21 11:23 2007
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

#include "RateAndCoalHmmStateAlphabet.h"

RateAndCoalHmmStateAlphabet::RateAndCoalHmmStateAlphabet(
    AverageCoalHmmStateAlphabet* hAlpha,
    bpp::DiscreteDistribution* rDist):
  bpp::AbstractParametrizable(""),
  coalHmmAlphabet_(hAlpha),
  rDist_(rDist),
  brlenParameters_(),
  nbCoalStates_(hAlpha->getNumberOfStates()),
  nbRates_(rDist->getNumberOfCategories()),
  nbStates_(nbCoalStates_ * nbRates_)
{
  addParameters_(hAlpha->getParameters());
  addParameters_(rDist->getIndependentParameters());
  // Trees and branch lengths:
  trees_.resize(nbStates_);
  species_ = hAlpha->getSpeciesNames();
  for (unsigned int i = 0; i < nbStates_; i++)
  {
    unsigned int s = (unsigned int)(i % nbCoalStates_);
    trees_[i] = coalHmmAlphabet_->getState(s).clone();
    brlenParameters_.push_back(coalHmmAlphabet_->getBranchLengthParametersForState(s));
  }

  applyParameters();
}

void RateAndCoalHmmStateAlphabet::applyParameters()
{
  coalHmmAlphabet_->matchParametersValues(getParameters());
  rDist_          ->matchParametersValues(getParameters());
  
  //Actualize trees:
  for (size_t i = 0; i < nbStates_; i++)
  {
    size_t s = static_cast<size_t>(i % nbCoalStates_);
    size_t r = i / nbCoalStates_;
    matchBranchLengths(trees_[i]->getRootNode(), coalHmmAlphabet_->getState(s).getRootNode(), rDist_->getCategory(r));
    brlenParameters_[i].matchParametersValues(coalHmmAlphabet_->getBranchLengthParametersForState(s));
    for (size_t j = 0; j < brlenParameters_[i].size(); j++)
    {
      brlenParameters_[i][j].setValue(brlenParameters_[i][j].getValue() * rDist_->getCategory(r));
      //cout << i << "\t" << j << "\t" << _brlenParameters[i][j]->getName() << "\t" << _brlenParameters[i][j].getValue() << endl;
    }
  }
}

void RateAndCoalHmmStateAlphabet::fireParameterChanged(const bpp::ParameterList& pl)
{
  applyParameters();
  std::vector<unsigned int> states(coalHmmAlphabet_->getNumberOfStates());
  for (unsigned int i = 0; i < states.size(); i++) states[i] = i;

  //Notify listeners:
  bpp::StateChangedEvent event(states);
  fireStateChangedEvent(event);
}

void RateAndCoalHmmStateAlphabet::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  coalHmmAlphabet_->printUserFriendlyParameters(out);
  bpp::ParameterList pl = rDist_->getIndependentParameters();
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    (out << "coal." << pl[i].getName() << " = " << pl[i].getValue()).endLine();
  }
}

void RateAndCoalHmmStateAlphabet::matchBranchLengths(bpp::Node* node1, const bpp::Node* node2, double scale)
{
  if(node1->getId() != node2->getId()) throw bpp::Exception("RateAndCoalHmmStateAlphabet::matchBranchLengths. Nodes are not identical (Id mismatch).");
  if(node1->getNumberOfSons() != node2->getNumberOfSons()) throw bpp::Exception("RateAndCoalHmmStateAlphabet::matchBranchLengths. Nodes are not identical (# sons mismatch).");
  if(node2->hasDistanceToFather())
    node1->setDistanceToFather(node2->getDistanceToFather() * scale);
  else
    node1->deleteDistanceToFather();
  for(unsigned int i = 0; i < node1->getNumberOfSons(); i++)
  {
    matchBranchLengths(node1->getSon(i), node2->getSon(i), scale);
  }
}

