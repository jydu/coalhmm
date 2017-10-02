//
// File: CtmcTools.h
// Created by: Thomas Mailund and Julien Dutheil
// Created on: Fri Apr 03 14:58 2009
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

#ifndef _CTMCTOOLS_H_
#define _CTMCTOOLS_H_

//From bpp-core:
#include <Bpp/Numeric/Matrix/Matrix.h>

using namespace bpp;

class CtmcTools
{
  public:
    static void singleCtmc(double recRate, double coalrate, double t, Matrix<double> & probs);

    /**
     * @brief Join the output probability vectors of H and C into an input vector for the two-species CTMC.
     */
    static void join(const Matrix<double> & h, const Matrix<double> & c, Matrix<double> & hc);

    /**
     * STATES:
     *
     * 0:  * *
     *     * *
     *
     * 1:  *-*
     *     * *
     *
     * 2:  * *
     *     *-*
     *
     * 3:  *-*
     *     *-*
     *
     * 4:    *
     *     x
     *       *
     *
     * 5:  / *
     *     x
     *       *
     *
     * 6:    *
     *     x
     *     \ *
     *
     * 7:  x-x
     *
     * 8:  *
     *       x
     *     *
     *
     * 9:  * \
     *       x
     *     *
     * 10: *
     *       x
     *     * /
     *
     * 11: x x
     * 
     *
     */
    static void twoRateMatrix(double recRate, double coalRate, Matrix<double> & generator);

    static void twoCtmc(double recrate, double coalRate, double t1, double t2, double coalRate2, Matrix<double> & probs);
};

#endif //_CTMCTOOLS_H_


