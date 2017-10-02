//
// File: TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet.h
// Created by: Julien Dutheil
// Created on: Thu Apr 02 10:30 2009
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

#ifndef _TWOSPECIESWITHOUTGROUPDISCRETIZEDCOALHMMSTATEALPHABET_H_
#define _TWOSPECIESWITHOUTGROUPDISCRETIZEDCOALHMMSTATEALPHABET_H_

#include "TwoSpeciesDiscretizedCoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Numeric/Prob/TruncatedExponentialDiscreteDistribution.h>


class TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet:
  public TwoSpeciesDiscretizedCoalHmmStateAlphabet,
  public bpp::AbstractParametrizable
{
  protected:
    std::vector<bpp::ParameterList> brlenParameters_;
    unsigned int nbClasses_;
    bpp::TruncatedExponentialDiscreteDistribution coalDist_;
    double a_;
    double b_;

  public:
    TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet(
        const std::string& species1, const std::string& species2, const std::string& outgroup,
        double tau1, double tau2, double theta12, unsigned int nbClasses, bool useMedian=false, double minTau1 = 0., double minTau2 = 0.0001, double minTheta = 0.0001);

    virtual ~TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet() {}

    virtual TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet * clone() const
    {
      return new TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet(*this);
    }

  public:
    const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex) const throw (bpp::HmmBadStateException)
    {
      if (stateIndex > nbClasses_) throw bpp::HmmBadStateException("TwoSpeciesWithOutgroupDiscretizedCoalHiddenStateAlphabet::getBranchLengthParametersForState [index=" + bpp::TextTools::toString(stateIndex) + "].");
      return brlenParameters_[stateIndex];
    }

    size_t getNumberOfStates() const { return nbClasses_; }

    double getProbabilityOfClass(size_t c) const { return 1. / static_cast<double>(nbClasses_); }

    void printUserFriendlyParameters(bpp::OutputStream& out) const;

    virtual double getTau1() const { return getParameter_(0).getValue(); }
    virtual double getTau2() const { return getParameter_(1).getValue(); }
    virtual double getTheta12() const { return getParameter_(2).getValue(); }

    /**
     * @return The discrete distribution of coalescent time within the ancestral population.
     */
    virtual const bpp::DiscreteDistribution& getCoalescentDistribution() const { return coalDist_; }

  protected:
    void fireParameterChanged(const bpp::ParameterList& pl);
    void applyParameters();

};

#endif //_TWOSPECIESWITHOUTGROUPDISCRETIZEDCOALHMMSTATEALPHABET_H_

