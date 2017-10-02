//
// File: HomogeneousCoalHmmEmissionProbabilities.h
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

#ifndef _HOMOGENEOUSCOALHMMPROBABILITIES_H_
#define _HOMOGENEOUSCOALHMMPROBABILITIES_H_

#include "CoalHmmEmissionProbabilities.h"

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

// From bpp-phyl:
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/ClockTreeLikelihood.h>

//From the STL:
#include <vector>
#include <string>

class HomogeneousCoalHmmEmissionProbabilities:
  public AbstractCoalHmmEmissionProbabilities
{
  protected:
    bpp::SubstitutionModel* model_;
    bpp::DiscreteDistribution* rDist_;
    std::vector<bpp::TreeLikelihood *> treeLikelihoods_;
    bpp::ParameterList modelParameters_;
    bpp::ParameterList rDistParameters_;
    std::vector<std::string> emissionParameterNames_;

  public:
    HomogeneousCoalHmmEmissionProbabilities(
        AverageCoalHmmStateAlphabet* stateAlphabet,
        const bpp::SiteContainer* data,
        bpp::SubstitutionModel* model,
        bpp::DiscreteDistribution* rDist,
        bool verbose = true);

    HomogeneousCoalHmmEmissionProbabilities(
        AverageCoalHmmStateAlphabet* stateAlphabet,
        bpp::SubstitutionModel* model,
        bpp::DiscreteDistribution* rDist,
        bool verbose = true);

    HomogeneousCoalHmmEmissionProbabilities(const HomogeneousCoalHmmEmissionProbabilities& ep):
      AbstractCoalHmmEmissionProbabilities(ep),
      model_(ep.model_),
      rDist_(ep.rDist_),
      treeLikelihoods_(ep.treeLikelihoods_.size()),
      modelParameters_(ep.modelParameters_),
      rDistParameters_(ep.rDistParameters_),
      emissionParameterNames_(ep.emissionParameterNames_)
    {
      for(unsigned int i = 0; i < treeLikelihoods_.size(); i++)
      {
        treeLikelihoods_[i]   = dynamic_cast<bpp::TreeLikelihood*>(ep.treeLikelihoods_[i]->clone());
      }
    }

    HomogeneousCoalHmmEmissionProbabilities& operator=(const HomogeneousCoalHmmEmissionProbabilities& ep)
    {
      AbstractCoalHmmEmissionProbabilities::operator=(ep);
      model_                  = ep.model_;
      modelParameters_        = ep.modelParameters_;
      rDist_                  = ep.rDist_;
      rDistParameters_        = ep.rDistParameters_;
      treeLikelihoods_.resize(ep.treeLikelihoods_.size());
      for(unsigned int i = 0; i < treeLikelihoods_.size(); i++)
      {
        treeLikelihoods_[i]   = dynamic_cast<bpp::TreeLikelihood*>(ep.treeLikelihoods_[i]->clone());
      }
      return *this;
    }

    virtual ~HomogeneousCoalHmmEmissionProbabilities()
    {
      for(unsigned int i = 0; i < treeLikelihoods_.size(); i++)
      {
        delete treeLikelihoods_[i];
      }
    }

    HomogeneousCoalHmmEmissionProbabilities* clone() const { return new HomogeneousCoalHmmEmissionProbabilities(*this); }

  public:

    const AverageCoalHmmStateAlphabet* getHmmStateAlphabet() const
    {
      return dynamic_cast<const AverageCoalHmmStateAlphabet*>(AbstractCoalHmmEmissionProbabilities::getHmmStateAlphabet());
    }

    void initialize() throw (bpp::Exception);

    void setParameters(const bpp::ParameterList& pl) throw (bpp::Exception)
    {
      setParametersValues(pl);
    }

    /**
     * @brief actualize substitution parameters and compute emission probabilities.
     */
    void fireParameterChanged(const bpp::ParameterList& pl);

    void printUserFriendlyParameters(bpp::OutputStream& out) const;

  protected:
    void init_(bool verbose);

};
    
#endif //_HOMOGENEOUSCOALHMMPROBABILITIES_H_

