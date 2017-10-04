//
// File: CoalHmmStateAlphabet.h
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

#ifndef _COALHMMSTATEALPHABET_H_
#define _COALHMMSTATEALPHABET_H_


// From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>

// From bpp-core:
#include <Bpp/BppVector.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Hmm/HmmStateAlphabet.h>

// From the STL:
#include <vector>

class CoalHmmStateAlphabet:
  public virtual bpp::HmmStateAlphabet
{
  protected:
    std::vector<bpp::StateListener*> listeners_;
    std::vector<std::string> species_;

  public:
    CoalHmmStateAlphabet() : listeners_(), species_() {}
    virtual ~CoalHmmStateAlphabet() {}

  public:
    bool worksWith(const HmmStateAlphabet* stateAlphabet) const { return stateAlphabet == this; }

    void addStateListener(bpp::StateListener* listener)
    {
      listeners_.push_back(listener);
    }

    virtual const std::vector<std::string>& getSpeciesNames() const { return species_; }

    /**
     * @brief Print parameter values.
     */
    virtual void printUserFriendlyParameters(bpp::OutputStream& out) const = 0;

  protected:
    void fireStateChangedEvent(bpp::StateChangedEvent& event)
    {
      for (size_t i = 0; i < listeners_.size(); i++)
        listeners_[i]->stateChanged(event);
    }
};

/**
 * @brief This abstract class designs models for which each hidden state is an average genealogy.
 */
class AverageCoalHmmStateAlphabet:
  public virtual CoalHmmStateAlphabet
{
  protected:
    std::vector<bpp::TreeTemplate<bpp::Node>*> trees_;

  public:
    AverageCoalHmmStateAlphabet() : trees_() {}
    virtual ~AverageCoalHmmStateAlphabet() {}

    AverageCoalHmmStateAlphabet* clone() const = 0;

  public:
    virtual const bpp::TreeTemplate<bpp::Node>& getState(size_t stateIndex) const throw (bpp::HmmBadStateException)
    {
      if (stateIndex >= trees_.size()) throw bpp::HmmBadStateException("CoalHmmStateAlphabet::getState (const) [index=." + bpp::TextTools::toString(stateIndex) + "].");
      return *trees_[stateIndex];
    }
    virtual bpp::TreeTemplate<bpp::Node>& getState(size_t stateIndex) throw (bpp::HmmBadStateException)
    {
      if (stateIndex >= trees_.size()) throw bpp::HmmBadStateException("CoalHmmStateAlphabet::getState [index=." + bpp::TextTools::toString(stateIndex) + "].");
      return *trees_[stateIndex];
    }

    size_t getNumberOfStates() const { return trees_.size(); }

    /**
     * @brief Get 'standard' parameter set for a state.
     *
     * This method return parameters of the kind "BrLen+NodeID", with values corresponding
     * to the current value of this alphabet parameters.
     *
     * @return A set of branch length parameters for the specified state.
     * @param stateIndex The state of interest.
     */
    virtual const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex) const = 0;

    /**
     * @brief Print parameter values.
     */
    virtual void printUserFriendlyParameters(bpp::OutputStream& out) const = 0;

};

/**
 * @brief This interface designs models for which hidden states is actually a set of genealogies with distinct banch lengths.
 * These models are currently deprecated, and were to be used with the 'corrected' HmmLikelihood objects.
 *
 * @warning: these models are not directly related to the other type of discretized models (called 'Timer' models),
 * for which each discretized genealogy is a full hidden states, and is hence an 'Average' coalescent model,
 * each of these genealogies begin actually an average genealogy for a class!
 *
 * @deprecated
 */
class DiscretizedCoalHmmStateAlphabet:
  public virtual CoalHmmStateAlphabet
{
  protected:
    std::vector< bpp::BppVector<bpp::TreeTemplate<bpp::Node>*> > trees_;

  public:
    DiscretizedCoalHmmStateAlphabet() : trees_() {}
    virtual ~DiscretizedCoalHmmStateAlphabet() {}

  public:
    virtual const bpp::BppVector<bpp::TreeTemplate<bpp::Node>*>& getState(size_t stateIndex) const throw (bpp::HmmBadStateException)
    {
      if (stateIndex >= trees_.size()) throw bpp::HmmBadStateException("DiscretizedCoalHmmStateAlphabet::getState (const) [index=." + bpp::TextTools::toString(stateIndex) + "].");
      return trees_[stateIndex];
    }
    virtual bpp::BppVector<bpp::TreeTemplate<bpp::Node>*>& getState(size_t stateIndex) throw (bpp::HmmBadStateException)
    {
      if (stateIndex >= trees_.size()) throw bpp::HmmBadStateException("DiscretizedCoalHmmStateAlphabet::getState [index=." + bpp::TextTools::toString(stateIndex) + "].");
      return trees_[stateIndex];
    }

    size_t getNumberOfStates() const { return trees_.size(); }

    virtual size_t getNumberOfClasses() const { return trees_[0].size(); }

    virtual double getProbabilityOfClass(size_t c) const = 0;

    /**
     * @brief Get 'standard' parameter set for a state.
     *
     * This method return parameters of the kind "BrLen+NodeID", with values corresponding
     * to the current value of this alphabet parameters.
     *
     * @return A set of branch length parameters for the specified state.
     * @param stateIndex The state of interest.
     * @param c          The class index.
     */
    virtual const bpp::ParameterList& getBranchLengthParametersForState(size_t stateIndex, size_t c) const = 0;

    /**
     * @brief Print parameter values.
     */
    virtual void printUserFriendlyParameters(bpp::OutputStream& out) const = 0;

};

#endif //_COALHMMSTATEALPHABET_H_

