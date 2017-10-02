//
// File: CtmcTools.cpp
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

#include "CtmcTools.h"

//From bpp-core:
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;

//From the STL:
#include <cmath>
using namespace std;

void CtmcTools::singleCtmc(double recRate, double coalRate, double t, Matrix<double> & probs)
{
  probs.resize(2, 2);
  double expt = exp(-(recRate + coalRate) * t);
  double relRec = recRate / (recRate + coalRate);
  double relCoal = coalRate / (recRate + coalRate);
  probs(0, 0) = relRec * expt + relCoal;
  probs(0, 1) = relRec - relRec * expt;
  probs(1, 0) = relCoal - relCoal * expt;
  probs(1, 1) = relCoal * expt + relRec;
}

void CtmcTools::join(const Matrix<double> & h, const Matrix<double> & c, Matrix<double> & hc)
{
  hc.resize(1,15);
  MatrixTools::fill(hc, 0);
  hc(0, 0) = h(0, 1) * c(0, 1); // * * x * *
  hc(0, 1) = h(0, 0) * c(0, 1); // *-* x * *
  hc(0, 2) = h(0, 1) * c(0, 0); // * * x *-*
  hc(0, 3) = h(0, 0) * c(0, 0); // *-* x *-*
}

void CtmcTools::twoRateMatrix(double recRate, double coalRate, Matrix<double> & generator)
{
  generator.resize(15, 15);
  MatrixTools::fill(generator, 0);


  // Beginning states
  generator(0, 0) = -6 * coalRate;
  generator(0, 1) = generator(0, 2) = generator(0, 4) = generator(0,5) = generator(0,6) = generator(0, 10) = coalRate;

  generator(1, 1) = -recRate -3 * coalRate;
  generator(1, 0) = recRate;
  generator(1, 3) = generator(1, 8) = generator(1, 11) = coalRate;
  
  generator(2, 2) = -recRate -3 * coalRate;
  generator(2, 0) = recRate;
  generator(2, 3) = generator(2, 9) = generator(2, 12) = coalRate;
  
  generator(3, 3) = -2 * recRate - coalRate;
  generator(3, 1) = generator(3, 2) = recRate;
  generator(3, 13) = coalRate;

  generator(4, 4) = -recRate - 3*coalRate;
  generator(4, 0) = recRate;
  generator(4, 6) = generator(4, 9) = generator(4, 11) = coalRate;

  generator(5, 5) = -recRate - 3*coalRate;
  generator(5, 0) = recRate;
  generator(5, 6) = generator(5, 8) = generator(5, 12) = coalRate;

  generator(6, 6) = -2*recRate - coalRate;
  generator(6, 4) = generator(6, 5) = recRate;
  generator(6, 13) = coalRate;


  // Left states
  generator(7, 7) = -3 * coalRate;
  generator(7, 8) = generator(7, 9) = generator(7, 14) = coalRate;

  generator(8, 8) = -recRate - coalRate;
  generator(8, 7) = recRate;
  generator(8, 13) = coalRate;

  generator(9, 9) = -recRate - coalRate;
  generator(9, 7) = recRate;
  generator(9, 13) = coalRate;


  // Right states
  generator(10, 10) = -3 * coalRate;
  generator(10, 11) = generator(10, 12) = generator(10, 14) = coalRate;

  generator(11, 11) = -recRate - coalRate;
  generator(11, 10) = recRate;
  generator(11, 13) = coalRate;

  generator(12, 12) = -recRate - coalRate;
  generator(12, 10) = recRate;
  generator(12, 13) = coalRate;


  // End states
  generator(13, 13) = -recRate;
  generator(13, 14) = recRate;

  generator(14, 14) = -coalRate;
  generator(14, 13) = coalRate;
}

void CtmcTools::twoCtmc(double recRate, double coalRate, double t1, double t2, double coalRate2, Matrix<double> & probs)
{
  RowMatrix<double> S1, S2, S, Q, D;

  //Single processes:
  singleCtmc(recRate, 1., t1, S1);
  singleCtmc(recRate, coalRate2, t1, S2);
  join(S1, S2, S);

  //Double process:
  twoRateMatrix(recRate, coalRate, Q);
  MatrixTools::scale(Q, t2);
  MatrixTools::exp(Q, D);

  MatrixTools::mult(S, D, probs);
}

