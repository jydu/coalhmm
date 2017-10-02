//
// File: CoalHmmSimulator.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 03 15:41 2008
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

#include "CoalHmmSimulator.h"

// From bpp-core:
#include <Bpp/App/ApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>

// from the STL:
#include <iostream>
#include <algorithm>
using namespace std;

CoalHmmSimulator::CoalHmmSimulator(
    CoalHmmTransitionMatrix* hiddenModel,
    bpp::SubstitutionModel* model,
    bpp::DiscreteDistribution* rDist,
    bool verbose):
  hiddenAlphabet_(dynamic_cast<const AverageCoalHmmStateAlphabet *>(hiddenModel->getHmmStateAlphabet())),
  hiddenModel_(hiddenModel),
  model_(model), 
  rDist_(rDist),
  simulators_(),
  genealogies_(),
  nbStates_(hiddenModel_->getNumberOfStates())
{
  init_(verbose);
}

/***************************************************************************************************************************/

void CoalHmmSimulator::init_(bool verbose)
{
  simulators_.resize(nbStates_);
  if (verbose)
  {
    ApplicationTools::displayTask("Initializing simulators", true);
  }

  for (size_t i = 0; i < nbStates_; i++)
  {
    if (verbose) ApplicationTools::displayGauge(i, nbStates_ - 1, '=');
      HomogeneousSequenceSimulator* sim = new HomogeneousSequenceSimulator(
        model_,
        rDist_,
        &dynamic_cast<const AverageCoalHmmStateAlphabet*>(hiddenAlphabet_)->getState(i));
      simulators_[i] = sim;
  }
  if (verbose) ApplicationTools::displayTaskDone();

}

/***************************************************************************************************************************/

SiteContainer* CoalHmmSimulator::simulate(size_t numberOfSites) const
{
  genealogies_.resize(numberOfSites);
  double r, cumprob;
  
  //Get the genealogies first:
  vector<double> probs = hiddenModel_->getEquilibriumFrequencies();
  genealogies_[0] = static_cast<size_t>(RandomTools::randMultinomial(1, probs)[0]);
  for( size_t i = 1; i < numberOfSites; i++)
  {
    r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
    cumprob = 0;
    bool test = true;
    for (size_t j = 0; test && j < nbStates_; j++)
    {
      cumprob += hiddenModel_->Pij(genealogies_[i-1],j);
		  if(r <= cumprob)
      {
        genealogies_[i] = j;
        test = false;
      }
    }
    //This should never happen!!!
    if(test) throw Exception("The unexpected happened in CoalHmmSimulator::simulate!");
  }

  //Now simulate sites according to the genealogies.
  //We should take car of the ordering of sequence which varies among genealogies!!!
  //So we parse the genealogies into chunks of a given type:
  vector<size_t> types;
  vector<size_t> sizes;
  size_t current = genealogies_[0];
  size_t count = 1;
  for(size_t i = 1; i < numberOfSites; i++)
  {
    size_t newgen = genealogies_[i];
    if(newgen == current) count++;
    else
    {
      //New segment:
      types.push_back(current);
      sizes.push_back(count);
      current = newgen;
      count = 1;
    }
  }
  //Add the last one:
  types.push_back(current);
  sizes.push_back(count);

  //Now build the containers:
  SiteContainer *sites1 = simulators_[types[0]]->simulate(sizes[0]);

  for (size_t i = 1; i < types.size(); i++)
  {
    SiteContainer *sites2 = simulators_[types[i]]->simulate(sizes[i]);
    SiteContainerTools::merge(*sites1, *sites2);
    delete sites2;
  }

  return sites1;
}

/***************************************************************************************************************************/

