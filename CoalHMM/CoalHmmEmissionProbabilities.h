//
// File: CoalHmmEmissionProbabilities.h
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:20 2007
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

#ifndef _COALHMMEMISSIONPROBABILITIES_H_
#define _COALHMMEMISSIONPROBABILITIES_H_

#include "CoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Hmm/HmmEmissionProbabilities.h>

// From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>

// From the STL:
#include <vector>

class CoalHmmEmissionProbabilities:
  public virtual bpp::HmmEmissionProbabilities
{
  public:
#ifndef NO_VIRTUAL_COV
    virtual CoalHmmEmissionProbabilities* clone() const = 0;
    virtual const CoalHmmStateAlphabet* getHmmStateAlphabet() const = 0;
#endif

};

class AbstractCoalHmmEmissionProbabilities:
  public virtual CoalHmmEmissionProbabilities,
  public bpp::AbstractParametrizable
{
  protected:
    /**
     * @brief the hidden states alphabet.
     */
    const CoalHmmStateAlphabet* hiddenAlphabet_;

    /**
     * @brief A pointer toward the compressed sequence data.
     */
    const bpp::SiteContainer* data_;

    /**
     * @brief A index vector containing for each position in the original alignment the actual position in the compressed data.
     */
    std::vector<size_t> siteIndex_;

    /**
     * @brief Store the emission probabilities for each pattern.
     */
    std::vector< std::vector<double> > emissions_;

    size_t nbStates_;

    bool isInitialized_;

  public:
    AbstractCoalHmmEmissionProbabilities(const CoalHmmStateAlphabet* hiddenAlphabet):
      AbstractParametrizable(""),
      hiddenAlphabet_(hiddenAlphabet),
      data_(0),
      siteIndex_(0),
      emissions_(0),
      nbStates_(hiddenAlphabet_->getNumberOfStates()),
      isInitialized_(false)
    {}
   
    AbstractCoalHmmEmissionProbabilities(const AbstractCoalHmmEmissionProbabilities& ep):
      AbstractParametrizable(ep),
      hiddenAlphabet_(ep.hiddenAlphabet_),
      data_(ep.data_),
      siteIndex_(ep.siteIndex_),
      emissions_(ep.emissions_),
      nbStates_(ep.nbStates_),
      isInitialized_(ep.isInitialized_)
    {}

    AbstractCoalHmmEmissionProbabilities& operator=(const AbstractCoalHmmEmissionProbabilities& ep)
    {
      AbstractParametrizable::operator=(ep);
      hiddenAlphabet_ = ep.hiddenAlphabet_;
      data_           = ep.data_;
      siteIndex_      = ep.siteIndex_;
      emissions_      = ep.emissions_;
      nbStates_       = ep.nbStates_;
      isInitialized_  = ep.isInitialized_;
      return *this;
    }

  public:
#ifndef NO_VIRTUAL_COV
    const CoalHmmStateAlphabet* getHmmStateAlphabet() const { return hiddenAlphabet_; }
#else
    const bpp::HmmStateAlphabet* getHmmStateAlphabet() const { return hiddenAlphabet_; }
#endif
    
    void setHmmStateAlphabet(const bpp::HmmStateAlphabet* hiddenAlphabet)
    {
      const CoalHmmStateAlphabet* cHiddenAlphabet = dynamic_cast<const CoalHmmStateAlphabet*>(hiddenAlphabet);
      if (!cHiddenAlphabet)
        throw bpp::HmmUnvalidAlphabetException("AbstractCoalHmmEmissionProbabilities::setHmmStateAlphabet. HmmStateAphabet should be a CoalHmmStateAlphabet.");
      hiddenAlphabet_ = cHiddenAlphabet;
    }

    void setData(const bpp::SiteContainer& data);
    
    bool isInitialized() const { return isInitialized_; }

    double getEmissionProbability(size_t pos, size_t stateIndex) const
    {
     if (stateIndex >= nbStates_) throw bpp::HmmBadStateException("AbstractCoalHmmEmissionProbabilities::getEmissionProbability. Unvalid state [index=" + bpp::TextTools::toString(stateIndex) + "].");
     if(pos >= siteIndex_.size()) throw bpp::Exception("AbstractCoalHmmEmissionProbabilities::getEmissionProbability. Bad site position: " + bpp::TextTools::toString(pos));
     return emissions_[siteIndex_[pos]][stateIndex];
    }

    double getLogEmissionProbability(size_t pos, size_t stateIndex) const
    {
     if (stateIndex >= nbStates_) throw bpp::HmmBadStateException("AbstractCoalHmmEmissionProbabilities::getEmissionProbability. Unvalid state [index=" + bpp::TextTools::toString(stateIndex) + "].");
     if(pos >= siteIndex_.size()) throw bpp::Exception("AbstractCoalHmmEmissionProbabilities::getEmissionProbability. Bad site position: " + bpp::TextTools::toString(pos));
     return log(emissions_[siteIndex_[pos]][stateIndex]);
    }

    void getEmissionProbabilities(size_t pos, std::vector<double>& probs) const
    {
      probs.resize(nbStates_);
      for(size_t i = 0; i < probs.size(); ++i)
      {
        probs[i] = getEmissionProbability(pos, i);
      }
    }

    void getEmissionProbabilitiesForEachPosition(std::vector< std::vector<double> >& probs, bool append = false) const
    {
      size_t offset = append ? probs.size() : 0;
      probs.resize(offset + siteIndex_.size());
      for(size_t i = 0; i < siteIndex_.size(); ++i)
      {
        getEmissionProbabilities(i, probs[offset + i]);
      }
    }

    void getLogEmissionProbabilities(size_t pos, std::vector<double>& probs) const
    {
      probs.resize(nbStates_);
      for (size_t i = 0; i < probs.size(); ++i)
      {
        probs[i] = getLogEmissionProbability(pos, i);
      }
    }

    void getLogEmissionProbabilitiesForEachPosition(std::vector< std::vector<double> >& probs, bool append = false) const
    {
      size_t offset = append ? probs.size() : 0;
      probs.resize(offset + siteIndex_.size());
      for (size_t i = 0; i < siteIndex_.size(); ++i)
      {
        getLogEmissionProbabilities(i, probs[offset + i]);
      }
    }

    size_t getNumberOfPositions() const
    {
      if (!data_)
        throw bpp::Exception("AbstractCoalHmmEmissionProbabilities::getNumberOfPositions. No data associated to this object.");
      return siteIndex_.size();
    }

    double operator()(size_t pos, size_t state) const { return emissions_[siteIndex_[pos]][state]; }

    const std::vector<double>& operator()(size_t pos) const { return emissions_[siteIndex_[pos]]; }

    size_t getObservedStateIndex(size_t pos) const { return siteIndex_[pos]; }
    
    void getEmissionProbabilities(std::vector< std::vector<double> >& probs) const { probs = emissions_; }

};

#endif //_COALHMMEMISSIONPROBABILITIES_H_

