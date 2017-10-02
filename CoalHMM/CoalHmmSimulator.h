//
// File: CoalHmmSimulator.h
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

#ifndef _COALHMMSIMULATOR_H_
#define _COALHMMSIMULATOR_H_

#include "CoalHmmStateAlphabet.h"
#include "CoalHmmTransitionMatrix.h"

//From bpp-core:
#include <Bpp/Clonable.h>

//From bpp-phyl:
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
using namespace bpp;

class CoalHmmSimulator:
public virtual Clonable
{
protected:
    /**
     * @brief The alphabet describing the hidden states.
     */
    const AverageCoalHmmStateAlphabet* hiddenAlphabet_;
    
    const CoalHmmTransitionMatrix* hiddenModel_;
    
    const bpp::SubstitutionModel* model_;
    const bpp::DiscreteDistribution* rDist_;
    std::vector<HomogeneousSequenceSimulator *> simulators_;
    mutable std::vector<size_t> genealogies_;
    size_t nbStates_;
    
public:
    CoalHmmSimulator(
                     CoalHmmTransitionMatrix* hiddenModel,
                     bpp::SubstitutionModel* model,
                     bpp::DiscreteDistribution* rDist,
                     bool verbose = true);
    
    CoalHmmSimulator(const CoalHmmSimulator& sim)
    : hiddenAlphabet_(sim.hiddenAlphabet_),
      hiddenModel_(0),
      model_(sim.model_),
      rDist_(sim.rDist_),
      simulators_(sim.simulators_.size()),
      genealogies_(sim.genealogies_),
      nbStates_()
    {
        for (size_t i = 0; i < simulators_.size(); i++)
        {
            simulators_[i] = dynamic_cast<HomogeneousSequenceSimulator*>(sim.simulators_[i]->clone());
        }
    }
    
    CoalHmmSimulator& operator=(const CoalHmmSimulator& sim)
    {
        hiddenAlphabet_         = sim.hiddenAlphabet_;
        model_                  = sim.model_;
        rDist_                  = sim.rDist_;
        nbStates_               = sim.nbStates_;
        simulators_.resize(sim.simulators_.size());
        for (size_t i = 0; i < simulators_.size(); i++)
        {
            simulators_[i] = dynamic_cast<HomogeneousSequenceSimulator *>(sim.simulators_[i]->clone());
        }
        genealogies_            = sim.genealogies_;
        
        return *this;
    }
    
    virtual ~CoalHmmSimulator()
    {
        for (size_t i = 0; i < simulators_.size(); i++)
        {
            delete simulators_[i];
        }
    }
    
#ifndef NO_VIRTUAL_COV
    CoalHmmSimulator*
#else
    Clonable*
#endif
    clone() const { return new CoalHmmSimulator(*this); }
    
public:
#ifndef NO_VIRTUAL_COV
    const CoalHmmStateAlphabet*
#else
    const HmmStateAlphabet*
#endif
    getHmmStateAlphabet() const { return hiddenAlphabet_; }
    
    bpp::SiteContainer* simulate(size_t numberOfSites) const;
    
    const std::vector<size_t>& getGenealogySequence() const { return genealogies_; }
    
protected:
    void init_(bool verbose);
    
};

#endif //_COALHMMSIMULATOR_H_

