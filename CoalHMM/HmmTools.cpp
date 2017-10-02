//
// File: HmmTools.cpp
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

#include "HmmTools.h"
#include "RateAndCoalHmmStateAlphabet.h"

// From bpp-core:
#include <Bpp/Numeric/VectorTools.h>

using namespace std;

/*************************************************************************************/

void HmmTools::getBestPosteriorHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, int separatorCode, bool append)
{
  if (!append) states.clear();
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);
  vector<size_t>::const_iterator brkpts = hmmLik.getBreakPoints().begin();
  size_t nextBrkPt = probs.size();
  if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;

  for (size_t i = 0; i < probs.size(); i++)
  {
    if (i == nextBrkPt)
    {
      states.push_back(separatorCode);
      nextBrkPt++;
      if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;
      else nextBrkPt = probs.size();
    }
    states.push_back(bpp::VectorTools::whichMax(probs[i]));
    //cout << i << "\t" << VectorTools::sum(probs[i]) << endl;
  }
}

void HmmTools::getBestPosteriorCoalHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, int separatorCode, bool append) throw (bpp::Exception)
{
  if (!append) states.clear();
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);
  vector<size_t>::const_iterator brkpts = hmmLik.getBreakPoints().begin();
  size_t nextBrkPt = probs.size();
  try
  {
    // Is it a Markov-Modulated model?
    const RateAndCoalHmmStateAlphabet* hAlpha = dynamic_cast<const RateAndCoalHmmStateAlphabet*>(&hmmLik.getHmmStateAlphabet());
    size_t n = hAlpha->getNumberOfStates();
    size_t m = hAlpha->getNumberOfCoalStates();

    if (!hAlpha) throw bpp::Exception("HmmTools::getBestPosteriorCoalHiddenStates. Unvalid alphabet.");
    vector<double> p(m);
    for (size_t i = 0; i < probs.size(); i++)
    {
      if (i == nextBrkPt)
      {
        states.push_back(separatorCode);
        nextBrkPt++;
        if(brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;
        else nextBrkPt = probs.size();
      }
      for (size_t j = 0; j < m; j++) p[j] = 0;
      for (size_t j = 0; j < n; j++) p[j % m] += probs[i][j];
      states.push_back(bpp::VectorTools::whichMax(p));
    }
  }
  catch (exception& e)
  {
    throw bpp::Exception("HmmTools::getBestPosteriorCoalHiddenStates(). Not a Markov-Modulated model.");
  }
}

void HmmTools::getBestPosteriorRateHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, int separatorCode, bool append) throw (bpp::Exception)
{
  if (!append) states.clear();
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);
  vector<size_t>::const_iterator brkpts = hmmLik.getBreakPoints().begin();
  size_t nextBrkPt = probs.size();
  try
  {
    // Is it a Markov-Modulated model?
    const RateAndCoalHmmStateAlphabet* hAlpha = dynamic_cast<const RateAndCoalHmmStateAlphabet*>(&hmmLik.getHmmStateAlphabet());
    size_t n = hAlpha->getNumberOfStates();
    size_t m = hAlpha->getNumberOfRateClasses();

    if (!hAlpha) throw bpp::Exception("HmmTools::getBestPosteriorRateHiddenStates. Unvalid alphabet.");
    vector<double> p(m);
    for (size_t i = 0; i < probs.size(); i++)
    {
      if (i == nextBrkPt)
      {
        states.push_back(separatorCode);
        nextBrkPt++;
        if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;
        else nextBrkPt = probs.size();
      }
      for (size_t j = 0; j < m; j++) p[j] = 0;
      for (size_t j = 0; j < n; j++) p[j / m] += probs[i][j];
      states.push_back(bpp::VectorTools::whichMax(p));
    }
  }
  catch (exception& e)
  {
    throw bpp::Exception("HmmTools::getBestPosteriorRateHiddenStates(). Not a Markov-Modulated model.");
  }
}

/*************************************************************************************/

void HmmTools::getBestPosteriorHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, double threshold, int unknownCode, int separatorCode, bool append)
{
  size_t offset = append ? states.size() : 0;
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);
  states.resize(offset + probs.size());
  vector<size_t>::const_iterator brkpts = hmmLik.getBreakPoints().begin();
  size_t nextBrkPt = probs.size();
  if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;

  for (size_t i = 0; i < probs.size(); i++)
  {
    if (i == nextBrkPt)
    {
      states.push_back(separatorCode);
      nextBrkPt++;
      if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;
      else nextBrkPt = probs.size();
    }
    double prob = bpp::VectorTools::max(probs[i]);
    if (prob >= threshold)
      states[offset + i] = bpp::VectorTools::whichMax(probs[i]);
    else
      states[offset + i] = unknownCode;
  }
}

void HmmTools::getBestPosteriorCoalHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, double threshold, int unknownCode, int separatorCode, bool append) throw (bpp::Exception)
{
  size_t offset = append ? states.size() : 0;
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);
  states.resize(offset + probs.size());
  try
  {
    // Is it a Markov-Modulated model?
    const RateAndCoalHmmStateAlphabet* hAlpha = dynamic_cast<const RateAndCoalHmmStateAlphabet*>(&hmmLik.getHmmStateAlphabet());
    size_t n = hAlpha->getNumberOfStates();
    size_t m = hAlpha->getNumberOfCoalStates();

    if (!hAlpha) throw bpp::Exception("HmmTools::getBestPosteriorCoalHiddenStates. Unvalid alphabet.");
    vector<double> p(m);
    for (size_t i = 0; i < probs.size(); i++)
    {
      for (size_t j = 0; j < m; j++) p[j] = 0;
      for (size_t j = 0; j < n; j++) p[j % m] += probs[i][j];
      double prob = bpp::VectorTools::max(p);
      if (prob >= threshold)
        states[offset + i] = bpp::VectorTools::whichMax(p);
      else
        states[offset + i] = unknownCode;
    }
  }
  catch (exception& e)
  {
    throw bpp::Exception("HmmTools::getBestPosteriorCoalHiddenStates(). Not a Markov-Modulated model.");
  }
}

void HmmTools::getBestPosteriorRateHiddenStates(const bpp::HmmLikelihood& hmmLik, std::vector<int>& states, double threshold, int unknownCode, int separatorCode, bool append) throw (bpp::Exception)
{
  size_t offset = append ? states.size() : 0;
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);
  states.resize(offset + probs.size());
  try
  {
    // Is it a Markov-Modulated model?
    const RateAndCoalHmmStateAlphabet* hAlpha = dynamic_cast<const RateAndCoalHmmStateAlphabet*>(&hmmLik.getHmmStateAlphabet());
    size_t n = hAlpha->getNumberOfStates();
    size_t m = hAlpha->getNumberOfRateClasses();

    if (!hAlpha) throw bpp::Exception("HmmTools::getBestPosteriorRateHiddenStates. Unvalid alphabet.");
    vector<double> p(m);
    for (size_t i = 0; i < probs.size(); i++)
    {
      for (size_t j = 0; j < m; j++) p[j] = 0;
      for (size_t j = 0; j < n; j++) p[j / m] += probs[i][j];
      double prob = bpp::VectorTools::max(p);
      if(prob >= threshold)
        states[offset + i] = bpp::VectorTools::whichMax(p);
      else
        states[offset + i] = unknownCode;
    }
  }
  catch (exception& e)
  {
    throw bpp::Exception("HmmTools::getBestPosteriorRateHiddenStates(). Not a Markov-Modulated model.");
  }
}

/*************************************************************************************/

void HmmTools::getPosteriorDivergences(const bpp::HmmLikelihood& hmmLik, bpp::Matrix<double>& divergences, vector<int>& blocks, vector<string>& namePairs) throw (bpp::Exception)
{
  const CoalHmmStateAlphabet& hAlpha = dynamic_cast<const CoalHmmStateAlphabet&>(hmmLik.getHmmStateAlphabet());
  size_t nbStates = hAlpha.getNumberOfStates();

  //first get posterior probabilities for all positions:
  vector< vector<double> > probs;
  hmmLik.getHiddenStatesPosteriorProbabilities(probs, false);

  //Then get divergence values for all hidden states:
  vector<string> names = hAlpha.getSpeciesNames();
  vector< vector<double> > statesDiv(nbStates);
  for (size_t i = 0; i < nbStates; ++i)
  {
    const bpp::Tree& tree = dynamic_cast<const bpp::Tree&>(hAlpha.getState(i));
    auto_ptr<bpp::DistanceMatrix> mat(bpp::TreeTools::getDistanceMatrix(tree));
    for (size_t j = 0; j < (names.size() - 1); ++j) {
      for (size_t k = j + 1; k < names.size(); ++k) {
        statesDiv[i].push_back((*mat)(names[j], names[k])); // do not use j,k directly, as species might not be in the same order for all states...
      }
    }
  }
 
  //Get column heads:
  namePairs.resize(0);
  for (size_t j = 0; j < (names.size() - 1); ++j) {
    for (size_t k = j + 1; k < names.size(); ++k) {
      namePairs.push_back(names[j] + "-" + names[k]);
    }
  }

  //Now make weighted average per position:
  divergences.resize(probs.size(), statesDiv[0].size());
  blocks.resize(probs.size());

  vector<size_t>::const_iterator brkpts = hmmLik.getBreakPoints().begin();
  size_t nextBrkPt = probs.size();
  if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;

  int currentBlock = 1;
  for (size_t i = 0; i < probs.size(); i++)
  {
    if (i == nextBrkPt)
    {
      currentBlock++;
      nextBrkPt++;
      if (brkpts != hmmLik.getBreakPoints().end()) nextBrkPt = *brkpts;
      else nextBrkPt = probs.size();
    }

    blocks[i] = currentBlock; 
    size_t p = 0;
    for (size_t j = 0; j < (names.size() - 1); ++j) {
      for (size_t k = j + 1; k < names.size(); ++k) {
        double d = 0;
        for (size_t l = 0; l < nbStates; ++l) {
          d += probs[i][l] * statesDiv[l][p];
        }
        divergences(i, p) = d;
        p++;
      }
    }
  }
}



/*************************************************************************************/

