//
// File: TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet.cpp
// Created by: Julien Dutheil
// Created on: Thu Apr 02 10:30 2009
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

#include "TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet.h"

TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet::TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet(
    const std::string& species1, const std::string& species2,
    double tau1, double theta12, size_t nbClasses, bool useMedian, double minTau, double minTheta) :
  AbstractParametrizable("coal."),
  brlenParameters_(),
  nbClasses_(nbClasses),
  coalDist_(nbClasses, 1. / theta12),
  a_()
{
  coalDist_.setMedian(useMedian);

  //These constraints are to avoid numerical issues. They comprise any reasonable biological values...
  minTau = std::max(0., minTau);
  minTheta = std::max(0., minTheta);
  addParameter_(new bpp::Parameter("coal.tau1"   , tau1   , new bpp::IntervalConstraint(1, minTau, true), true));
  addParameter_(new bpp::Parameter("coal.theta12", theta12, new bpp::IntervalConstraint(1, minTheta, true), true));

  species_.push_back(species1);
  species_.push_back(species2);
  
  //Now build trees:
  trees_.resize(1);
  bpp::Node* root, * n1, * n2;
  bpp::TreeTemplate<bpp::Node> * tree;
  ////species topology:
  trees_.resize(nbClasses_);
  brlenParameters_.resize(nbClasses_);
  root = new bpp::Node(0);
  n1 = new bpp::Node(1, species1);
  n2 = new bpp::Node(2, species2);
  root->addSon(n1);
  root->addSon(n2);
  tree = new bpp::TreeTemplate<bpp::Node>(root);
  tree->resetNodesId();

  for(unsigned int i = 0; i < nbClasses_; i++)
  {
    trees_[i] = tree->clone();
    brlenParameters_[i].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n1->getId()), 0));
    brlenParameters_[i].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n2->getId()), 0));
  }
  delete tree;

  applyParameters();
}

void TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet::applyParameters()
{
  double tau1    = getParameter_(0).getValue();
  double theta12 = getParameter_(1).getValue();
  bpp::ParameterList pl = coalDist_.getParameters();
  pl.setParameterValue("Exponential.lambda", 1. / theta12);
  coalDist_.setParametersValues(pl);
  
  for(unsigned int i = 0; i < nbClasses_; i++)
  {
    a_  = tau1 + coalDist_.getCategory(i);
    a_  = std::max(0.000001, a_);

    //Actualize trees:
    bpp::TreeTemplate<bpp::Node> * tree = trees_[i];
    tree->getRootNode()->getSon(0)->setDistanceToFather(a_);
    tree->getRootNode()->getSon(1)->setDistanceToFather(a_);
     
    //Actualize branch lengths parameters:
    bpp::ParameterList * brlenParameters = &brlenParameters_[i];
    (*brlenParameters)[0].setValue(a_);
    (*brlenParameters)[1].setValue(a_);
  }
}

void TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet::fireParameterChanged(const bpp::ParameterList& pl)
{
  applyParameters();
  std::vector<unsigned int> states(nbClasses_);
  for (unsigned int i = 0; i < nbClasses_; i++)
    states[i] = 0;

  //Notify listeners:
  bpp::StateChangedEvent event(states);
  fireStateChangedEvent(event);
}

void TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  double tau1    = getTau1();
  double theta12 = getTheta12();
  (out << "tau1    = " << tau1).endLine();
  (out << "theta12 = " << theta12).endLine();
}

