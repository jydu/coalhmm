//
// File: CoalHmmTransitionMatrix.h
// Created by: Julien Dutheil
// Created on: Mon Oct 29 09:39 2007
//

// This file is part of the CoalHMM program and library.
//
// CoalHMM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.
//
// CoalHMM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CoalHMM.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _COALHMMTRANSITIONMATRIX_H_
#define _COALHMMTRANSITIONMATRIX_H_

#include "CoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Hmm/HmmTransitionMatrix.h>
#include <Bpp/Numeric/AbstractParametrizable.h>

// From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>

class CoalHmmTransitionMatrix :
  public virtual bpp::HmmTransitionMatrix
{
public:
#ifndef NO_VIRTUAL_COV
  const CoalHmmStateAlphabet* getHmmStateAlphabet() const = 0;
#endif

  /**
   * @brief Print parameter values.
   */
  virtual void printUserFriendlyParameters(bpp::OutputStream& out) const = 0;
};

class AbstractCoalHmmTransitionMatrix :
  public virtual CoalHmmTransitionMatrix,
  public bpp::AbstractParametrizable
{
private:
  const CoalHmmStateAlphabet* hiddenAlphabet_;
  size_t nbStates_;

protected:
  bpp::RowMatrix<double> transitions_;
  std::vector<double> freqs_;

public:
  AbstractCoalHmmTransitionMatrix(const CoalHmmStateAlphabet* hiddenAlphabet) :
    bpp::AbstractParametrizable("coal."),
    hiddenAlphabet_(hiddenAlphabet),
    nbStates_(hiddenAlphabet->getNumberOfStates()),
    transitions_(),
    freqs_()
  {}

  AbstractCoalHmmTransitionMatrix(const AbstractCoalHmmTransitionMatrix& tm):
    bpp::AbstractParametrizable(tm),
    hiddenAlphabet_(tm.hiddenAlphabet_),
    nbStates_(tm.nbStates_),
    transitions_(tm.transitions_),
    freqs_(tm.freqs_)
  {}

  AbstractCoalHmmTransitionMatrix& operator=(const AbstractCoalHmmTransitionMatrix& tm)
  {
    AbstractParametrizable::operator=(tm);
    hiddenAlphabet_ = tm.hiddenAlphabet_;
    nbStates_ = tm.nbStates_;
    transitions_ = tm.transitions_;
    freqs_ = tm.freqs_;
    return *this;
  }

public:
#ifndef NO_VIRTUAL_COV
  const CoalHmmStateAlphabet*
#else
  const HmmStateAlphabet*
#endif
  getHmmStateAlphabet() const { return hiddenAlphabet_; }

  void setHmmStateAlphabet(const bpp::HmmStateAlphabet* hiddenAlphabet)
  {
    const CoalHmmStateAlphabet* cHiddenAlphabet = dynamic_cast<const CoalHmmStateAlphabet*>(hiddenAlphabet);
    if (!cHiddenAlphabet)
      throw bpp::HmmUnvalidAlphabetException("AbstractCoalHmmTransitionMatrix::setHmmStateAlphabet. HmmStateAphabet should be a CoalHmmStateAlphabet.");
    hiddenAlphabet_ = cHiddenAlphabet;
  }

  size_t getNumberOfStates() const { return nbStates_; }

  double Pij(size_t i, size_t j) const { return transitions_(i, j); }

  const bpp::Matrix<double>& getPij() const { return transitions_; }

  const std::vector<double>& getEquilibriumFrequencies() const { return freqs_; }
};

#endif // _COALHMMTRANSITIONMATRIX_H_

