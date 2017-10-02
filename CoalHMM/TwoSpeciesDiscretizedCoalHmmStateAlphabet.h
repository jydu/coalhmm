//
// File: TwoSpeciesDiscretizedCoalHmmStateAlphabet.h
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

#ifndef _TWOSPECIESDISCRETIZEDCOALHMMSTATEALPHABET_H_
#define _TWOSPECIESDISCRETIZEDCOALHMMSTATEALPHABET_H_

#include "CoalHmmStateAlphabet.h"

//From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

/**
 * @brief Common interface for the two-species models.
 */
class TwoSpeciesDiscretizedCoalHmmStateAlphabet:
  public AverageCoalHmmStateAlphabet
{
  public:

    virtual ~TwoSpeciesDiscretizedCoalHmmStateAlphabet() {}

    virtual TwoSpeciesDiscretizedCoalHmmStateAlphabet* clone() const = 0;

  public:
    virtual double getTau1() const = 0;
    virtual double getTheta12() const = 0;

    /**
     * @return The discrete distribution of coalescent time within the ancestral population.
     */
    virtual const bpp::DiscreteDistribution& getCoalescentDistribution() const = 0;

};

#endif //_TWOSPECIESDISCRETIZEDCOALHMMSTATEALPHABET_H_

