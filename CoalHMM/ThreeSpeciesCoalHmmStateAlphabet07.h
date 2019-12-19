//
// File: ThreeSpeciesCoalHmmStateAlphabet07.h
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

#ifndef _THREESPECIESCOALHMMSTATEALPHABET07_H_
#define _THREESPECIESCOALHMMSTATEALPHABET07_H_

#include "CoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Numeric/AbstractParametrizable.h>

class ThreeSpeciesCoalHmmStateAlphabet07:
  public AverageCoalHmmStateAlphabet,
  public bpp::AbstractParametrizable
{
  public:
    std::vector<bpp::ParameterList> brlenParameters_;
    bool reparam_;

  public:
    ThreeSpeciesCoalHmmStateAlphabet07(
        const std::string& species1, const std::string& species2, const std::string& species3, const std::string& outgroup,
        double a, double b, double c, double a2, bool reparam = true);

    virtual ~ThreeSpeciesCoalHmmStateAlphabet07() {}

    ThreeSpeciesCoalHmmStateAlphabet07* clone() const { return new ThreeSpeciesCoalHmmStateAlphabet07(*this); }

  public:
    const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex) const
    {
      if (stateIndex > 3) throw bpp::HmmBadStateException("ThreeSpeciesHiddenStateAlphabet::getBranchLengthParametersForState [index=" + bpp::TextTools::toString(stateIndex) + "].");
      return brlenParameters_[stateIndex];
    }

    void printUserFriendlyParameters(bpp::OutputStream& out) const;

  protected:
    void fireParameterChanged(const bpp::ParameterList& pl);
    void applyParameters();

};

#endif //_THREESPECIESCOALHMMSTATEALPHABET07_H_

