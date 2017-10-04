//
// File: TwoSpeciesDiscretizedCoalHmmTransitionMatrix.h
// Created by: Julien Dutheil
// Created on: Thu Apr 02 13:20 2009
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

#ifndef _TWOSPECIESDISCRETIZEDCOALHMMTRANSITIONMATRIX_H_
#define _TWOSPECIESDISCRETIZEDCOALHMMTRANSITIONMATRIX_H_

#include "CoalHmmTransitionMatrix.h"
#include "TwoSpeciesDiscretizedCoalHmmStateAlphabet.h"

class OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix:
  public AbstractCoalHmmTransitionMatrix
{
  private:
    double tau1_, rho_, theta1_, theta2_, theta12_;
    std::vector<unsigned int> noneCoalescedIndices_, bothCoalescedIndices_, leftCoalescedRightNotIndices_;
    bool extentThetas_;
    
  public:
    OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double theta1, double theta2, double minRho = 0.0001, double minTheta = 0.00001);
    OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double minRho = 0.0001);
    
    OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix* clone() const { return new OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(*this); }
    
  public:
    void fireParameterChanged(const bpp::ParameterList& pl)
    {
        actualizeTransitionMatrix_();
    }
    
    void printUserFriendlyParameters(bpp::OutputStream& out) const;
    
  protected:
    void actualizeTransitionMatrix_();
    static double probSubset_(const bpp::Matrix<double>& pi, const std::vector<unsigned int>& indices);
    static void distSubset_(const bpp::Matrix<double>& pi, const std::vector<unsigned int>& indices, bpp::Matrix<double>& out);
    
  private:
    void initialize_();
};



class TwoSpeciesDiscretizedCoalHmmTransitionMatrix:
  public AbstractCoalHmmTransitionMatrix
{
  private:
    double tau1_, rho_, theta1_, theta2_, theta12_;
    bool extentThetas_;
    
  public:
    TwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double theta1, double theta2, double minRho = 0.0001, double minTheta = 0.00001);
    TwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double minRho = 0.0001);
    
    TwoSpeciesDiscretizedCoalHmmTransitionMatrix* clone() const { return new TwoSpeciesDiscretizedCoalHmmTransitionMatrix(*this); }
    
  public:
    void fireParameterChanged(const bpp::ParameterList& pl)
    {
      actualizeTransitionMatrix_();
    }
    
    void printUserFriendlyParameters(bpp::OutputStream& out) const;
    
  protected:
    void actualizeTransitionMatrix_();
    
  private:
    void initialize_();

    static inline size_t idx(size_t i, size_t j) { return j * (j + 1) / 2 + i; }
};

#endif //_TWOSPECIESDISCRETIZEDCOALHMMTRANSITIONMATRIX_H_

