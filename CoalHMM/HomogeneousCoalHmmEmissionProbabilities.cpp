//
// File: HomogeneousCoalHmmEmissionProbabilities.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:57 2007
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

#include "HomogeneousCoalHmmEmissionProbabilities.h"

// From bpp-core:
#include <Bpp/App/ApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>

using namespace bpp;

// from the STL:
#include <iostream>
#include <algorithm>
using namespace std;

HomogeneousCoalHmmEmissionProbabilities::HomogeneousCoalHmmEmissionProbabilities(
    AverageCoalHmmStateAlphabet* stateAlphabet,
    bpp::SubstitutionModel* model,
    bpp::DiscreteDistribution* rDist,
    bool verbose):
  AbstractCoalHmmEmissionProbabilities(stateAlphabet),
  model_(model), rDist_(rDist),
  treeLikelihoods_(), modelParameters_(), rDistParameters_(), emissionParameterNames_()
{
  init_(verbose);
}

HomogeneousCoalHmmEmissionProbabilities::HomogeneousCoalHmmEmissionProbabilities(
    AverageCoalHmmStateAlphabet* stateAlphabet,
    const bpp::SiteContainer* data,
    bpp::SubstitutionModel* model,
    bpp::DiscreteDistribution* rDist,
    bool verbose):
  AbstractCoalHmmEmissionProbabilities(stateAlphabet),
  model_(model), rDist_(rDist),
  treeLikelihoods_(), modelParameters_(), rDistParameters_(), emissionParameterNames_()
{
  init_(verbose);
  setData(*data);
}

void HomogeneousCoalHmmEmissionProbabilities::init_(bool verbose)
{
  treeLikelihoods_.resize(nbStates_);
  if (verbose)
  {
    ApplicationTools::displayTask("Initializing likelihoods", true);
  }

  for (unsigned int i = 0; i < nbStates_; i++)
  {
    if (verbose) ApplicationTools::displayGauge(i, nbStates_ - 1, '=');
    RHomogeneousTreeLikelihood* tl = new RHomogeneousTreeLikelihood(
        dynamic_cast<const AverageCoalHmmStateAlphabet *>(hiddenAlphabet_)->getState(i),
        model_,
        rDist_,
        false, false, true);
    tl->setMinimumBranchLength(0.);
    treeLikelihoods_[i] = tl;
  }
  if(verbose) ApplicationTools::displayTaskDone();

  //Init parameters:
  ParameterList mpl = model_->getIndependentParameters();
  for(unsigned int j = 0; j < mpl.size(); j++)
  {
    emissionParameterNames_.push_back(mpl[j].getName());
  }
  modelParameters_ = mpl;
  addParameters_(mpl);
  
  ParameterList rpl = rDist_->getIndependentParameters();
  for(unsigned int j = 0; j < rpl.size(); j++)
  {
    emissionParameterNames_.push_back(rpl[j].getName());
  }
  rDistParameters_ = rpl;
  addParameters_(rpl);
}

void HomogeneousCoalHmmEmissionProbabilities::initialize() throw (bpp::Exception)
{
  if (!data_)
    throw Exception("HomogeneousCoalHmmEmissionProbabilities::initialize(). No data associated to this instance.");
  if (isInitialized_)
    throw Exception("HomogeneousCoalHmmEmissionProbabilities::initialize(). Instance already initialized.");

  ApplicationTools::displayTask("Initializing emissions", true);
  for (unsigned int i = 0; i < nbStates_; i++)
  {
    ApplicationTools::displayGauge(i, nbStates_ - 1, '=');
    treeLikelihoods_[i]->setData(*data_);
    treeLikelihoods_[i]->initialize();
  }
  ApplicationTools::displayTaskDone();

  isInitialized_ = true;
  fireParameterChanged(getParameters());
}

/***************************************************************************************************************************/

void HomogeneousCoalHmmEmissionProbabilities::fireParameterChanged(const bpp::ParameterList& pl)
{
  modelParameters_.matchParametersValues(pl);
  rDistParameters_.matchParametersValues(pl);
  for (unsigned int i = 0; i < nbStates_; i++)
  {
    ParameterList p = pl;
    p.addParameters(getHmmStateAlphabet()->getBranchLengthParametersForState(i));
    treeLikelihoods_[i]->matchParametersValues(p);
  }

  for (unsigned int i = 0; i < emissions_.size(); i++)
  {
    for (unsigned int j = 0; j < nbStates_; j++)
    {
      emissions_[i][j] = treeLikelihoods_[j]->getLikelihoodForASite(i);
    }
  }
}

/***************************************************************************************************************************/

void HomogeneousCoalHmmEmissionProbabilities::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  for (unsigned int i = 0; i < modelParameters_.size(); i++)
  {
    (out << modelParameters_[i].getName() << " = " << modelParameters_[i].getValue()).endLine();
  }
  for (unsigned int i = 0; i < rDistParameters_.size(); i++)
  {
    (out << rDistParameters_[i].getName() << " = " << rDistParameters_[i].getValue()).endLine();
  }
}

/***************************************************************************************************************************/

