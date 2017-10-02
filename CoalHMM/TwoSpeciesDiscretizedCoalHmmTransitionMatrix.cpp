//
// File: TwoSpeciesDiscretizedCoalHmmTransitionMatrix.cpp
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

#include <cassert>

#include "TwoSpeciesDiscretizedCoalHmmTransitionMatrix.h"
//This is for the old version:
#include "CtmcTools.h"
//This is for the new one:
#include "generated-two-species/generators.h"
#include "generated-two-species/statesets.cpp"

//From bpp-core:
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;

//From the STL:
#include <deque>
using namespace std;

OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double theta1, double theta2, double minRho, double minTheta):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet),
  tau1_(), rho_(rho), theta1_(theta1), theta2_(theta2), theta12_(),
  noneCoalescedIndices_(), bothCoalescedIndices_(), leftCoalescedRightNotIndices_(),
  extentThetas_(true)
{
  addParameter_(new Parameter("coal.theta1", theta1_, new IntervalConstraint(1, minTheta, true), true));
  addParameter_(new Parameter("coal.theta2", theta2_, new IntervalConstraint(1, minTheta, true), true));
  addParameter_(new Parameter("coal.rho", rho_, new IntervalConstraint(1, 0.0001, true), true));
  initialize_();
}

OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double minRho):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet),
  tau1_(), rho_(rho), theta1_(0), theta2_(0), theta12_(),
  noneCoalescedIndices_(), bothCoalescedIndices_(), leftCoalescedRightNotIndices_(),
  extentThetas_(false)
{
  addParameter_(new Parameter("coal.rho", rho_, new IntervalConstraint(1, minRho, true), true)); 
  initialize_();
}


void OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::initialize_()
{
  unsigned int n = getNumberOfStates();
  transitions_.resize(n, n);
  freqs_.resize(n);
  noneCoalescedIndices_.resize(7); 
  noneCoalescedIndices_[0] = 0; 
  noneCoalescedIndices_[1] = 1; 
  noneCoalescedIndices_[2] = 2; 
  noneCoalescedIndices_[3] = 3; 
  noneCoalescedIndices_[4] = 4; 
  noneCoalescedIndices_[5] = 5; 
  noneCoalescedIndices_[6] = 6; 

  bothCoalescedIndices_.resize(2);
  bothCoalescedIndices_[0] = 13;
  bothCoalescedIndices_[1] = 14;

  leftCoalescedRightNotIndices_.resize(3);  
  leftCoalescedRightNotIndices_[0] = 7;
  leftCoalescedRightNotIndices_[1] = 8;  
  leftCoalescedRightNotIndices_[2] = 9;

  actualizeTransitionMatrix_();
}

double OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::probSubset_(const Matrix<double>& pi, const vector<unsigned int>& indices)
{
  deque<unsigned int> sortedIndices(indices.begin(), indices.end());
  //sort(sortedIndices.begin(), sortedIndices.end()); We assume they are sorted!
  double sum = 0;
  for(unsigned int i = 0; i < pi.getNumberOfColumns(); i++)
  {
    if(i == sortedIndices[0])
    {
      sum += pi((unsigned int)0, i);
      sortedIndices.pop_front();
    }
  }
  return sum;
}

void OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::distSubset_(const Matrix<double>& pi, const vector<unsigned int>& indices, Matrix<double>& out)
{
  out.resize(1,15);
  MatrixTools::fill(out, 0.);
  double s = 0;
  for(unsigned int i = 0; i < indices.size(); i++)
  {
    s+= (out((unsigned int)0, indices[i]) = pi((unsigned int)0, indices[i]));
  }
  for(unsigned int i = 0; i < 15; i++)
    out((unsigned int)0, i) /= s;
}

#define CONSISTENCY_CHECK 0

inline bool checkDouble(double a, double b, double precision=1e-5)
{
  if (abs(a-b) < precision) return true;
  else return false;
}

inline void checkRowSum(const Matrix<double> &M, double expected, int lineno)
{
#if CONSISTENCY_CHECK
  for (unsigned int r = 0; r < M.nRows(); ++r)
    {
      double rowSum = 0.0;
      for (unsigned int c = 0; c < M.nCols(); ++c)
  rowSum += M(r,c);
      if (!checkDouble(rowSum, expected))
  {
    std::cerr << "Row sum consistency violated at line "
        << lineno << std::endl
        << "Expected " << expected << " but got "
        << rowSum << std::endl;
    assert(false);
  }
    }
#endif
}

inline void checkSum(const Matrix<double> &M, double expected, int lineno)
{
#if CONSISTENCY_CHECK
  double sum = 0.0;
  for (unsigned int r = 0; r < M.nRows(); ++r)
    {
      for (unsigned int c = 0; c < M.nCols(); ++c)
  sum += M(r,c);
    }
  if (!checkDouble(sum, expected))
    {
      std::cerr << "Matrix sum consistency violated at line "
    << lineno << std::endl
    << "Expected " << expected << " but got "
    << sum << std::endl;
      assert(false);
    }
#endif
}


void OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::actualizeTransitionMatrix_()
{
  tau1_    = dynamic_cast<const TwoSpeciesDiscretizedCoalHmmStateAlphabet*>(getHmmStateAlphabet())->getTau1(); 
  theta12_ = dynamic_cast<const TwoSpeciesDiscretizedCoalHmmStateAlphabet*>(getHmmStateAlphabet())->getTheta12();
  if(extentThetas_)
  {
    theta1_  = getParameter_(0).getValue();
    theta2_  = getParameter_(1).getValue();
    rho_     = getParameter_(2).getValue();
  }
  else
  {
    theta1_  = theta12_;
    theta2_  = theta12_;
    rho_     = getParameter_(0).getValue();
  }

  RowMatrix<double> H, C, Q, P, Pi, joint, tmp;
  //These ones are actually vectors:
  RowMatrix<double> pi_0, pi_i, pi_j, pi_i_minus_one, pi_j_minus_one, pi_i_minus_one_none, pi_j_plus_one, pi_i_left, pi_n; 

  //Single species output:
  CtmcTools::singleCtmc(rho_, 1./theta1_, tau1_, H);
  CtmcTools::singleCtmc(rho_, 1./theta2_, tau1_, C);
  checkRowSum(H, 1.0, __LINE__);
  checkRowSum(C, 1.0, __LINE__);

  //Two species rate matrix:
  CtmcTools::twoRateMatrix(rho_, 1./theta12_, Q);
  checkRowSum(Q, 0.0, __LINE__);

  //Distribution at the beginning
  //of the two species process:
  CtmcTools::join(H, C, pi_0);
  checkRowSum(pi_0, 1.0, __LINE__);

  vector<double> breakPoints = dynamic_cast<const TwoSpeciesDiscretizedCoalHmmStateAlphabet*>(getHmmStateAlphabet())->getCoalescentDistribution().getBounds();
  double tau2 = *breakPoints.rbegin();

  bool isTruncated = tau2 < 200; //otherwise it's not significantly different from a "full" exponential distribution.
  //// first get the joint probabilities

  // all but the last...
  unsigned int n = getNumberOfStates();
  joint.resize(n, n);
  MatrixTools::fill(joint, 0);

  //double norm = 1.;
  if (isTruncated)
  {
    //norm = 1. - exp(-tau2 / theta12_);
    for (unsigned int i = 0; i < n; i++)
    {
      tmp = Q;
      MatrixTools::scale(tmp, breakPoints[i]);
      MatrixTools::exp(tmp, P);
      checkRowSum(P, 1.0, __LINE__);

      tmp = Q;
      MatrixTools::scale(tmp, breakPoints[i+1] - breakPoints[i]);
      MatrixTools::exp(tmp, Pi);
      checkRowSum(Pi, 1.0, __LINE__);

      //distribution when entering state i...
      MatrixTools::mult(pi_0, P, pi_i_minus_one);
      checkRowSum(pi_i_minus_one, 1.0, __LINE__);

      //distribution assuming none were coalesced before state i...
      double probNone = probSubset_(pi_i_minus_one, noneCoalescedIndices_);
      distSubset_(pi_i_minus_one, noneCoalescedIndices_, pi_i_minus_one_none);
      checkRowSum(pi_i_minus_one_none, 1.0, __LINE__);

      MatrixTools::mult(pi_i_minus_one_none, Pi, pi_i);
      checkRowSum(pi_i, 1.0, __LINE__);

      //both coalescing in state i...
      joint(i, i) = probSubset_(pi_i, bothCoalescedIndices_) * probNone;

      double probLeft = probSubset_(pi_i, leftCoalescedRightNotIndices_);
      distSubset_(pi_i, leftCoalescedRightNotIndices_, pi_i_left); 
      checkRowSum(pi_i_left, 1.0, __LINE__);

      for (unsigned int j = i + 1; j < n; j++)
      {
        //right coalescing in state j that is *not* the final state
        tmp = Q;
        MatrixTools::scale(tmp, breakPoints[j] - breakPoints[i+1]);
        MatrixTools::exp(tmp, P);
  checkRowSum(P, 1.0, __LINE__);

        MatrixTools::mult(pi_i_left, P, pi_j_minus_one);
  checkRowSum(pi_j_minus_one, 1.0, __LINE__);

        double probNotRight = probSubset_(pi_j_minus_one, leftCoalescedRightNotIndices_);
        distSubset_(pi_j_minus_one, leftCoalescedRightNotIndices_, pi_j);
  checkRowSum(pi_j, 1.0, __LINE__);

        tmp = Q;
        MatrixTools::scale(tmp, breakPoints[j+1] - breakPoints[j]);
        MatrixTools::exp(tmp, P);
  checkRowSum(P, 1.0, __LINE__);

        MatrixTools::mult(pi_j, P, pi_j_plus_one);
  checkRowSum(pi_j_plus_one, 1.0, __LINE__);

        double prob_j = probSubset_(pi_j_plus_one, bothCoalescedIndices_);

        joint(i, j) = prob_j * probNotRight * probLeft * probNone;
        joint(j, i) = joint(i, j); //symmetries...
      }         
    }
  }
  else
  {
    for(unsigned int i = 0; i < n - 1; i++)
    {
      tmp = Q;
      MatrixTools::scale(tmp, breakPoints[i]);
      MatrixTools::exp(tmp, P);
      checkRowSum(P, 1.0, __LINE__);

      tmp = Q;
      MatrixTools::scale(tmp, breakPoints[i+1] - breakPoints[i]);
      MatrixTools::exp(tmp, Pi);

#if 0
      MatrixTools::print(Q);
      std::cout << (breakPoints[i+1] - breakPoints[i]) << std::endl;
      MatrixTools::print(Pi);
#endif
      checkRowSum(Pi, 1.0, __LINE__);

      //distribution when entering state i...
      MatrixTools::mult(pi_0, P, pi_i_minus_one);
      checkRowSum(pi_i_minus_one, 1.0, __LINE__);

      //distribution assuming none were coalesced before state i...
      double probNone = probSubset_(pi_i_minus_one, noneCoalescedIndices_);
      distSubset_(pi_i_minus_one, noneCoalescedIndices_, pi_i_minus_one_none);
      checkRowSum(pi_i_minus_one_none, 1.0, __LINE__);

      MatrixTools::mult(pi_i_minus_one_none, Pi, pi_i);
      checkRowSum(pi_i, 1.0, __LINE__);

      //both coalescing in state i...
      joint(i, i) = probSubset_(pi_i, bothCoalescedIndices_) * probNone;

      double probLeft = probSubset_(pi_i, leftCoalescedRightNotIndices_);
      distSubset_(pi_i, leftCoalescedRightNotIndices_, pi_i_left); 
      checkRowSum(pi_i_left, 1.0, __LINE__);

      for(unsigned int j = i + 1; j < n - 1; j++)
      {
        //right coalescing in state j that is *not* the final state
        tmp = Q;
        MatrixTools::scale(tmp, breakPoints[j] - breakPoints[i+1]);
        MatrixTools::exp(tmp, P);
  checkRowSum(P, 1.0, __LINE__);

        MatrixTools::mult(pi_i_left, P, pi_j_minus_one);
  checkRowSum(pi_j_minus_one, 1.0, __LINE__);

        double probNotRight = probSubset_(pi_j_minus_one, leftCoalescedRightNotIndices_);
        distSubset_(pi_j_minus_one, leftCoalescedRightNotIndices_, pi_j);
  checkRowSum(pi_j, 1.0, __LINE__);

        tmp = Q;
        MatrixTools::scale(tmp, breakPoints[j+1] - breakPoints[j]);
        MatrixTools::exp(tmp, P);
  checkRowSum(P, 1.0, __LINE__);

        MatrixTools::mult(pi_j, P, pi_j_plus_one);
  checkRowSum(pi_j_plus_one, 1.0, __LINE__);

        double prob_j = probSubset_(pi_j_plus_one, bothCoalescedIndices_);

        joint(i, j) = prob_j * probNotRight * probLeft * probNone;
        joint(j, i) = joint(i, j); //symmetries...
      }
            
      // right coalescing in the final state
      tmp = Q;
      MatrixTools::scale(tmp, breakPoints[n-1] - breakPoints[i+1]);
      MatrixTools::exp(tmp, P);
      checkRowSum(P, 1.0, __LINE__);

      MatrixTools::mult(pi_i_left, P, pi_n);
      checkRowSum(pi_n, 1.0, __LINE__);
        
      joint(i, n-1) = probSubset_(pi_n, leftCoalescedRightNotIndices_) * probLeft * probNone;
      joint(n-1, i) = joint(i, n-1); // symmetries...
    }


    // joint for the last state
    tmp = Q;
    MatrixTools::scale(tmp, breakPoints[n-1]);
    MatrixTools::exp(tmp, P);
    checkRowSum(P, 1.0, __LINE__);
    
    MatrixTools::mult(pi_0, P, pi_n);
    checkRowSum(pi_n, 1.0, __LINE__);

    joint(n-1, n-1) = probSubset_(pi_n, noneCoalescedIndices_);
  }

  //Now get the transition probabilities
  // using that the breakpoints give us
  // equal probability for the left
  // nucleotide to coalesce in each
  // state...
  transitions_ = joint;
  //MatrixTools::scale(transitions_, static_cast<double>(nbClasses_) / norm);
  for (unsigned int i = 0; i < n; i++)
  {
    double sum = 0.;
    for (unsigned int j = 0; j < n; j++)
      sum += transitions_(i, j);
    for (unsigned int j = 0; j < n; j++)
      transitions_(i, j) /= sum;
  }
  checkRowSum(transitions_, 1.0, __LINE__);

 
  //Perform some testing:
  for(unsigned int i = 0; i < n; i++)
  {
    double s = 0;
    for(unsigned int j = 0; j < n; j++)
    {
      s += transitions_(i,j);
    }
    if(abs(s-1.) > 0.01)
    {
      cout << "ERROR!!! Probabilities do not sum to one in transition matrix:" << endl;
      MatrixTools::print(transitions_);
      cout << "tau1 = " << tau1_ << endl;
      cout << "theta1 = " << theta1_ << endl;
      cout << "theta2 = " << theta2_ << endl;
      cout << "theta12 = " << theta12_ << endl;
      cout << "rho = " << rho_ << endl;
    }
  }

  //CHECKME
  for (unsigned int i = 0; i < n; i++)
    freqs_[i] = 1. / static_cast<double>(n);
}

void OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix::printUserFriendlyParameters(OutputStream& out) const
{
  (out << "theta1  = " << theta1_).endLine();
  (out << "theta2  = " << theta2_).endLine();
  (out << "theta12 = " << theta12_).endLine();
  (out << "rho     = " << rho_).endLine();
}




TwoSpeciesDiscretizedCoalHmmTransitionMatrix::TwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double theta1, double theta2, double minRho, double minTheta):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet), 
  tau1_(), rho_(rho), theta1_(theta1), theta2_(theta2), theta12_(), extentThetas_(true)
{
  addParameter_(new Parameter("coal.theta1", theta1_, new IntervalConstraint(1, minTheta, true), true));
  addParameter_(new Parameter("coal.theta2", theta2_, new IntervalConstraint(1, minTheta, true), true));
  addParameter_(new Parameter("coal.rho", rho_, new IntervalConstraint(1, 0.0001, true), true));
  initialize_();
}

TwoSpeciesDiscretizedCoalHmmTransitionMatrix::TwoSpeciesDiscretizedCoalHmmTransitionMatrix(const TwoSpeciesDiscretizedCoalHmmStateAlphabet* hiddenAlphabet, double rho, double minRho):
  AbstractCoalHmmTransitionMatrix(hiddenAlphabet),
  tau1_(), rho_(rho), theta1_(0), theta2_(0), theta12_(), extentThetas_(false)
{
  addParameter_(new Parameter("coal.rho", rho_, new IntervalConstraint(1, minRho, true), true)); 
  initialize_();
}

void TwoSpeciesDiscretizedCoalHmmTransitionMatrix::initialize_()
{
  unsigned int n = getNumberOfStates();
  transitions_.resize(n, n);
  freqs_.resize(n);
  actualizeTransitionMatrix_();
}

void TwoSpeciesDiscretizedCoalHmmTransitionMatrix::actualizeTransitionMatrix_()
{
  tau1_    = dynamic_cast<const TwoSpeciesDiscretizedCoalHmmStateAlphabet*>(getHmmStateAlphabet())->getTau1(); 
  theta12_ = dynamic_cast<const TwoSpeciesDiscretizedCoalHmmStateAlphabet*>(getHmmStateAlphabet())->getTheta12();
  if (extentThetas_)
  {
    theta1_  = getParameter_(0).getValue();
    theta2_  = getParameter_(1).getValue();
    rho_     = getParameter_(2).getValue();
  }
  else
  {
    theta1_  = theta12_;
    theta2_  = theta12_;
    rho_     = getParameter_(0).getValue();
  }

  double R    = rho_;
  double C_H  = 1. / theta1_;
  double C_C  = 1. / theta2_;
  double C_HC = 1. / theta12_;
  
  RowMatrix<double> Q1 = bcs::getRateMatrix_C_H(C_C, C_H, R);
  RowMatrix<double> Q2 = bcs::getRateMatrix_HC(C_HC, R);
  
  vector<double> breakpoints = dynamic_cast<const TwoSpeciesDiscretizedCoalHmmStateAlphabet*>(getHmmStateAlphabet())->getCoalescentDistribution().getBounds();
  breakpoints.erase(--breakpoints.end()); //Removes most upper bound, assumed to be infty
  unsigned int n = breakpoints.size();
  
  vector<RowMatrix<double> > P(idx(0, n + 1));
  unsigned int CTMCsize = Q2.getNumberOfRows();
  
  // Fill the diagonal with identities and one above the diagonal with exponentiated matrices
  for (unsigned int i = 0; i <= n; ++i) {
    MatrixTools::getId(CTMCsize, P[idx(i, i)]);
  }
  for (unsigned int i = 0; i < n - 1; ++i) {
    RowMatrix<double> tmp(Q2);
    MatrixTools::scale(tmp, breakpoints[i + 1] - breakpoints[i]);
    MatrixTools::exp(tmp, P[idx(i, i + 1)]);
  }
  
  // Handle non-diagonal matrices that are not single intervals by multiplying together
  // transitions from individual intervals.  This is simply basic CTMC theory: P[i,i+2]=P[i,i+1]P[i+1,i+2].
  for (unsigned int i = 0; i < n - 2; ++i) {
    for (unsigned int j = i + 2; j < n; ++j) {
      MatrixTools::mult(P[idx(i, j - 1)], P[idx(j - 1, j)], P[idx(i, j)]);
    }
  }
    
  // Make the probability of going to an end-state 1 in the last interval
  // It does not matter which end-state we end up in
  for (unsigned int i = 0; i < n; ++i) {
    P[idx(i, n)].resize(CTMCsize, CTMCsize);
    MatrixTools::fill(P[idx(i, n)], 0.0);
    for (unsigned int s = 0; s < CTMCsize; ++s) {
      P[idx(i, n)](s, HC__Ee[0]) = 1.0;
    }
  }
    
  
  // Get probability vector at speciation time...
  RowMatrix<double> tmp(Q1);
  MatrixTools::scale(tmp, breakpoints[0]);
  RowMatrix<double> exptmp;
  MatrixTools::exp(tmp, exptmp);
  vector<double> pi0(CTMCsize);
  for (unsigned int s = 0; s < Q1.getNumberOfColumns(); ++s) {
    pi0[s] = exptmp(::start, s);
  }
    
  // Make joint probability matrix
  RowMatrix<double> joint; joint.resize(n, n);
    
  // First the diagonal
  for (unsigned int i = 0; i < n; ++i) {
    Matrix<double>& P0i = P[idx(0,i)];
    Matrix<double>& Piip1 = P[idx(i,i+1)];
    double p = 0.0;
    for (unsigned int b1 = 0; b1 < no_HC__Bb; ++b1) {
      double part1 = pi0[HC__Bb[b1]];
      for (unsigned int b2 = 0; b2 < no_HC__Bb; ++b2) {
        double part2 = part1 * P0i(HC__Bb[b1],HC__Bb[b2]);
        for (unsigned int e = 0; e < no_HC__Ee; ++e) {
          p += part2 * Piip1(HC__Bb[b2],HC__Ee[e]);
        }
      }
    }
    joint(i,i) = p;
  }
    
  // Then the upper triangle
  for (unsigned int i = 0; i < n; ++i) {
    Matrix<double>& P0i = P[idx(0,i)];
    Matrix<double>& Piip1 = P[idx(i,i+1)];
    for (unsigned int j = i+1; j < n; ++j) {
      Matrix<double>& Pip1j = P[idx(i+1,j)];
      Matrix<double>& Pjjp1 = P[idx(j,j+1)];
      double p = 0.0;
      for (unsigned int b1 = 0; b1 < no_HC__Bb; ++b1) {
        double part1 = pi0[HC__Bb[b1]];
        for (unsigned int b2 = 0; b2 < no_HC__Bb; ++b2) {
          double part2 = part1 * P0i(HC__Bb[b1],HC__Bb[b2]);
          for (unsigned int l1 = 0; l1 < no_HC__Eb; ++l1) {
            double part3 = part2 * Piip1(HC__Bb[b2],HC__Eb[l1]);
            for (unsigned int l2 = 0; l2 < no_HC__Eb; ++l2) {
              double part4 = part3 * Pip1j(HC__Eb[l1],HC__Eb[l2]);
              for (unsigned int e = 0; e < no_HC__Ee; ++e) {
                p += part4 * Pjjp1(HC__Eb[l2],HC__Ee[e]);
              }
            }
          }
        }
      }
      joint(i,j) = p;
    }
  }
  // ... and then the lower triangle by symmetry
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      joint(i,j) = joint(j,i);
    }
  }
      
  // Finally, get transition probability matrix from joint probability matrix
  for (unsigned int i = 0; i < n; ++i) {
    double sum = 0.0;
    for (unsigned int j = 0; j < n; ++j) {
      sum += joint(i,j);
    }
    for (unsigned int j = 0; j < n; ++j) {
      transitions_(i,j) = joint(i,j)/sum;
    }
  }

  //CHECKME
  for (unsigned int i = 0; i < n; i++)
    freqs_[i] = 1. / (double)n;
}

void TwoSpeciesDiscretizedCoalHmmTransitionMatrix::printUserFriendlyParameters(OutputStream& out) const
{
  (out << "theta1  = " << theta1_).endLine();
  (out << "theta2  = " << theta2_).endLine();
  (out << "theta12 = " << theta12_).endLine();
  (out << "rho     = " << rho_).endLine();
}


