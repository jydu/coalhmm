//
// File: HmmTools.h
// Created by: Julien Dutheil
// Created on: Tue Nov 01 13:03 2007
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

#ifndef _HMMTOOLS_H_
#define _HMMTOOLS_H_

#include <Bpp/Numeric/Hmm/HmmLikelihood.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

// From the STL:
#include <vector>
using namespace std;

class HmmTools
{
  public:
    static void getBestPosteriorHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, int separatorCode, bool append = false);
    static void getBestPosteriorHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, double threshold, int unknownCode, int separatorCode, bool append = false);
    static void getBestPosteriorCoalHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, int separatorCode, bool append = false);
    static void getBestPosteriorCoalHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, double threshold, int unknownCode, int separatorCode, bool append = false);
    static void getBestPosteriorRateHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, int separatorCode, bool append = false);
    static void getBestPosteriorRateHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, double threshold, int unknownCode, int separatorCode, bool append = false);
    static void getPosteriorDivergences(const bpp::HmmLikelihood& hmmLik, bpp::Matrix<double>& divergences, vector<int>& blocks, vector<string>& namePairs);
};

#endif //_HMMTOOLS_H_

