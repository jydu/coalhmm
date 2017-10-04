//
// File: RateAndCoalHmmStateAlphabet.h
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

#ifndef _RATEANDCOALHMMSTATEALPHABET_H_
#define _RATEANDCOALHMMSTATEALPHABET_H_

#include "CoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

//From the STL:
#include <memory>

class RateAndCoalHmmStateAlphabet:
  public AverageCoalHmmStateAlphabet,
  public bpp::AbstractParametrizable
{
  public:
    std::unique_ptr<AverageCoalHmmStateAlphabet> coalHmmAlphabet_;
    std::unique_ptr<bpp::DiscreteDistribution> rDist_;
    std::vector<bpp::ParameterList> brlenParameters_;
    size_t nbCoalStates_;
    size_t nbRates_;
    size_t nbStates_;

  public:
    /**
     * @brief Build a new rate modulated coalhmm alphabet.
     *
     * The alphabet will own the two objects given as input.
     */
    RateAndCoalHmmStateAlphabet(
        AverageCoalHmmStateAlphabet* hAlpha,
        bpp::DiscreteDistribution* rDist) throw (bpp::Exception);

    RateAndCoalHmmStateAlphabet(const RateAndCoalHmmStateAlphabet& alpha):
      AverageCoalHmmStateAlphabet(alpha),
      bpp::AbstractParametrizable(alpha),
      coalHmmAlphabet_(alpha.coalHmmAlphabet_->clone()),
      rDist_(alpha.rDist_->clone()),
      brlenParameters_(alpha.brlenParameters_),
      nbCoalStates_(alpha.nbCoalStates_),
      nbRates_(alpha.nbRates_),
      nbStates_(alpha.nbStates_)
    {}

    RateAndCoalHmmStateAlphabet& operator=(const RateAndCoalHmmStateAlphabet& alpha)
    {
      AverageCoalHmmStateAlphabet::operator=(alpha);
      bpp::AbstractParametrizable::operator=(alpha);
      coalHmmAlphabet_ = std::auto_ptr<AverageCoalHmmStateAlphabet>(alpha.coalHmmAlphabet_->clone());
      rDist_           = std::auto_ptr<bpp::DiscreteDistribution>(alpha.rDist_->clone());
      brlenParameters_ = alpha.brlenParameters_;
      nbCoalStates_    = alpha.nbCoalStates_;
      nbRates_         = alpha.nbRates_;
      nbStates_        = alpha.nbStates_;
      return *this;
    }

    virtual ~RateAndCoalHmmStateAlphabet() {}

    RateAndCoalHmmStateAlphabet* clone() const { return new RateAndCoalHmmStateAlphabet(*this); }

  public:
    const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex) const throw (bpp::HmmBadStateException)
    {
      if (stateIndex >= nbStates_) throw bpp::HmmBadStateException("RateAndCoalHmmStateAlphabet::getBranchLengthParametersForState [index=" + bpp::TextTools::toString(stateIndex) + "].");
      return brlenParameters_[stateIndex];
    }

    size_t getNumberOfRateClasses() const { return nbRates_; }
    
    size_t getNumberOfCoalStates() const { return nbCoalStates_; }

    const bpp::DiscreteDistribution& getRateDistribution() const { return *rDist_.get(); }
    bpp::DiscreteDistribution& getRateDistribution() { return *rDist_.get(); }

    void printUserFriendlyParameters(bpp::OutputStream& out) const;

  protected:
    void fireParameterChanged(const bpp::ParameterList& pl);
    void applyParameters();

    //Recursive method
    static void matchBranchLengths(bpp::Node* node1, const bpp::Node* node2, double scale) throw (bpp::Exception);

};

#endif //_RATEANDCOALHMMSTATEALPHABET_H_

