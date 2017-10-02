//
// File: ThreeSpeciesCoalHmmStateAlphabet.h
// Created by: Julien Dutheil
// Created on: Thu Jun 12 18:01 2008
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

#ifndef _THREESPECIESCOALHMMSTATEALPHABET_H_
#define _THREESPECIESCOALHMMSTATEALPHABET_H_

#include "CoalHmmStateAlphabet.h"

// From NumCalc:

class ThreeSpeciesCoalHmmStateAlphabet:
  public virtual CoalHmmStateAlphabet
{
  public:
    virtual ~ThreeSpeciesCoalHmmStateAlphabet() {}

#ifndef NO_VIRTUAL_COV
    virtual ThreeSpeciesCoalHmmStateAlphabet* clone() const = 0;
#endif

  public:
    virtual double getTau1() const = 0;
    virtual double getTau2() const = 0;
    virtual double getTau3() const = 0;
    virtual double getTheta1() const = 0;
    virtual double getTheta2() const = 0;

};

#endif //_THREESPECIESCOALHMMSTATEALPHABET_H_

