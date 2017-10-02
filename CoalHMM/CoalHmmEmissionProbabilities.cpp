//
// File: CoalHmmEmissionProbabilities.cpp
// Created by: Julien Dutheil
// Created on: Mon Oct 29 16:54 2007
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

#include "CoalHmmEmissionProbabilities.h"

//From bpp-core:
#include <Bpp/App/ApplicationTools.h>

//From bpp-seq:
#include <Bpp/Seq/SiteTools.h>

//From bpp-phyl:
#include <Bpp/Phyl/PatternTools.h>

using namespace bpp;

void AbstractCoalHmmEmissionProbabilities::setData(const bpp::SiteContainer& data) throw (bpp::Exception)
{
  //Initialize index:
  unsigned int nbSites = data.getNumberOfSites();
  siteIndex_.resize(nbSites);

  if (! VectorTools::haveSameElements(hiddenAlphabet_->getSpeciesNames(), data.getSequencesNames()))
    throw Exception("AbstractCoalHmmEmissionProbabilities::setData(). Sequence names must be the same as species names in the alphabet, and there should not be any additional sequence.");
  SiteContainer* zipSites = PatternTools::shrinkSiteSet(data);
  unsigned int nbDistinctSites = zipSites->getNumberOfSites();
  ApplicationTools::displayResult("Number of patterns", nbDistinctSites);

  ApplicationTools::displayTask("Assigning patterns", true);
  for (unsigned int i = 0; i < nbSites; i++)
  {
    ApplicationTools::displayGauge(i, nbSites - 1, '=');
    const Site* site1 = &data.getSite(i);
    for (unsigned int ii = 0; ii < nbDistinctSites; ii++)
    {
      if (SiteTools::areSitesIdentical(zipSites->getSite(ii), *site1))
      {
        siteIndex_[i] = ii;
        break;
      }
    }
  }
  data_ = zipSites;
  ApplicationTools::displayTaskDone();

  //Emission array:
  emissions_.resize(nbDistinctSites);
  for(unsigned int i = 0; i < nbDistinctSites; i++)
    emissions_[i].resize(nbStates_);

  isInitialized_ = false;
}

