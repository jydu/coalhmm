//
// File: ThreeSpeciesCoalHmmStateAlphabet08.h
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:07 2007
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

#ifndef _THREESPECIESCOALHMMSTATEALPHABET09_H_
#define _THREESPECIESCOALHMMSTATEALPHABET09_H_

#include "CoalHmmStateAlphabet.h"
#include "ThreeSpeciesCoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/NumTools.h>

class ThreeSpeciesCoalHmmStateAlphabet09:
  public AverageCoalHmmStateAlphabet,
  public ThreeSpeciesCoalHmmStateAlphabet,
  public bpp::AbstractParametrizable
{
  protected:
    std::vector<bpp::ParameterList> brlenParameters_;
    double tau1_, tau2_, tau3_, a_, b_, c_, a2_, b2_, c2_;
    short parametrization_;
    bool useMedian_;

  public:
    ThreeSpeciesCoalHmmStateAlphabet09(
        const std::string& species1, const std::string& species2, const std::string& species3, const std::string& outgroup,
        double tau1, double tau2, double c2, double theta1, double theta2, short parametrization = 1, bool useMedian = false, double minTau = 0.00001, double minTheta = 0.0001, double maxTau=1, double maxTheta=1);

    virtual ~ThreeSpeciesCoalHmmStateAlphabet09() {}

    ThreeSpeciesCoalHmmStateAlphabet09* clone() const
    {
      return new ThreeSpeciesCoalHmmStateAlphabet09(*this);
    }

  public:
    const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex) const
    {
      if (stateIndex > 3)
        throw bpp::HmmBadStateException("ThreeSpeciesHmmStateAlphabet09::getBranchLengthParametersForState [index=" + bpp::TextTools::toString(stateIndex) + "].");
      return brlenParameters_[stateIndex];
    }

    void printUserFriendlyParameters(bpp::OutputStream& out) const;

    double getTau1() const { return tau1_; }
    double getTau2() const { return tau2_; }
    double getTau3() const { return tau3_; }
    double getTheta1() const { return getParameter_(3).getValue(); }
    double getTheta2() const { return getParameter_(4).getValue(); }

  protected:
    void fireParameterChanged(const bpp::ParameterList& pl);
    void applyParameters();

};

#endif //_THREESPECIESCOALHMMSTATEALPHABET09_H_

