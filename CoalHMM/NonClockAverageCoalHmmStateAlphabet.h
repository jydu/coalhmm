//
// File: NonClockAverageCoalHmmStateAlphabet.h
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

#ifndef _NONCLOCKAVERAGECOALHMMSTATEALPHABET_H_
#define _NONCLOCKAVERAGECOALHMMSTATEALPHABET_H_

#include "CoalHmmStateAlphabet.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

/**
 * @brief This alphabet is a wrapper class that takes as input another AverageCoalHmmStateAlphabet, and adds mutational parameters allowing for different substitution and/or sequencing error rates.
 * Up to two branches are allowed to have differential rates.
 */
class NonClockAverageCoalHmmStateAlphabet:
  public AverageCoalHmmStateAlphabet,
  public bpp::AbstractParametrizable
{
  public:
    AverageCoalHmmStateAlphabet* alphabet_;
    std::string speciesWithError1_;
    std::string speciesWithError2_;
    double error1_;
    double error2_;
    bool hasSecondError_;
    std::vector<bpp::ParameterList> brlenParameters_;

  public:
    NonClockAverageCoalHmmStateAlphabet(
        AverageCoalHmmStateAlphabet* alphabet,
        const std::string& species,
       	double error = 0);

    NonClockAverageCoalHmmStateAlphabet(
        AverageCoalHmmStateAlphabet* alphabet,
        const std::string& species1,
        const std::string& species2,
       	double error1 = 0,
       	double error2 = 0);

    ~NonClockAverageCoalHmmStateAlphabet() { delete alphabet_; }

    NonClockAverageCoalHmmStateAlphabet(const NonClockAverageCoalHmmStateAlphabet& model) :
      AverageCoalHmmStateAlphabet(model),
      AbstractParametrizable(model),
      alphabet_(model.alphabet_->clone()),
      speciesWithError1_(model.speciesWithError1_),
      speciesWithError2_(model.speciesWithError2_),
      error1_(model.error1_),
      error2_(model.error2_),
      hasSecondError_(model.hasSecondError_),
      brlenParameters_(model.brlenParameters_)
    {}

    NonClockAverageCoalHmmStateAlphabet& operator=(const NonClockAverageCoalHmmStateAlphabet& model)
    {
      AverageCoalHmmStateAlphabet::operator=(model);
      AbstractParametrizable::operator=(model);
      alphabet_          = model.alphabet_->clone();
      speciesWithError1_ = model.speciesWithError1_;
      speciesWithError2_ = model.speciesWithError2_;
      error1_            = model.error1_;
      error2_            = model.error2_;
      hasSecondError_    = model.hasSecondError_;
      brlenParameters_   = model.brlenParameters_;
      return *this;
    }

    NonClockAverageCoalHmmStateAlphabet* clone() const {
      return new NonClockAverageCoalHmmStateAlphabet(*this);
    }

  public:
    bool worksWith(const HmmStateAlphabet* stateAlphabet) const {
      return stateAlphabet == this || stateAlphabet == alphabet_;
    }
    const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex) const {
      if (stateIndex > alphabet_->getNumberOfStates())
        throw bpp::HmmBadStateException("NonClockAverageCoalHmmStateAlphabet::getBranchLengthParametersForState [index=" + bpp::TextTools::toString(stateIndex) + "].");
      return brlenParameters_[stateIndex];
    }
    void printUserFriendlyParameters(bpp::OutputStream& out) const;
    void fireParameterChanged(const bpp::ParameterList& params);

    
};

#endif //_NONCLOCKAVERAGECOALHMMSTATEALPHABET_H_

