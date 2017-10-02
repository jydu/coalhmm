//
// File: ThreeSpeciesCoalHmmStateAlphabet09.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 26 12:10 2007
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

#include "ThreeSpeciesCoalHmmStateAlphabet09.h"

ThreeSpeciesCoalHmmStateAlphabet09::ThreeSpeciesCoalHmmStateAlphabet09(
    const std::string& species1, const std::string& species2, const std::string& species3, const std::string& outgroup,
    double tau1, double tau2, double c2, double theta1, double theta2, short parametrization, bool useMedian, double minTau, double minTheta, double maxTau, double maxTheta) :
  AbstractParametrizable("coal."), 
  brlenParameters_(),
  tau1_(tau1), tau2_(tau2), tau3_(-1), a_(-1), b_(-1), c_(-1), a2_(-1), b2_(-1), c2_(-1),
  parametrization_(parametrization),
  useMedian_(useMedian)
{
  //These constraints are to avoid numerical issues. They comprise any reasonable biological values...
  minTau = std::max(0., minTau);
  minTheta = std::max(0., minTheta);
  if (parametrization == 1) {
    addParameter_(new bpp::Parameter("coal.tau1"  , tau1  , new bpp::IntervalConstraint(minTau, maxTau, true, true)    , true));
    addParameter_(new bpp::Parameter("coal.tau2"  , tau2  , new bpp::IntervalConstraint(minTau, maxTau, true, true)    , true));
  } else {
    addParameter_(new bpp::Parameter("coal.sigma1", tau1 / (tau1 + tau2), new bpp::IntervalConstraint(minTau / (tau1 + tau2), 1. - minTau / (tau1 + tau2), true, true), true));
    addParameter_(new bpp::Parameter("coal.tau12" , tau1 + tau2         , new bpp::IntervalConstraint(minTau, maxTau, true, true), true));
  }
  addParameter_(new bpp::Parameter("coal.c2"    , c2    , new bpp::IntervalConstraint(minTau, maxTau, true, true)    , true));
  addParameter_(new bpp::Parameter("coal.theta1", theta1, new bpp::IntervalConstraint(minTheta, maxTheta, true, true), true));
  addParameter_(new bpp::Parameter("coal.theta2", theta2, new bpp::IntervalConstraint(minTheta, maxTheta, true, true), true));

  species_.push_back(species1);
  species_.push_back(species2);
  species_.push_back(species3);
  species_.push_back(outgroup);

  //Now build trees:
  trees_.resize(4);
  brlenParameters_.resize(4);
  bpp::Node* root, * nOut, * n1, * n2, * n3, * n4, * n5;
  ////species topology:
  root = new bpp::Node(0);
  nOut = new bpp::Node(6, outgroup);
  n1 = new bpp::Node(1);
  n2 = new bpp::Node(2);
  n3 = new bpp::Node(3, species1);
  n4 = new bpp::Node(4, species2);
  n5 = new bpp::Node(5, species3);
  n2->addSon(n3);
  n2->addSon(n4);
  n1->addSon(n2);
  n1->addSon(n5);
  root->addSon(n1);
  root->addSon(nOut);
  trees_[0] = new bpp::TreeTemplate<bpp::Node>(root);
  trees_[0]->resetNodesId();
  brlenParameters_[0].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n1->getId()), 0));
  brlenParameters_[0].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n2->getId()), 0));
  brlenParameters_[0].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n3->getId()), 0));
  brlenParameters_[0].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n4->getId()), 0));
  brlenParameters_[0].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n5->getId()), 0));
  brlenParameters_[0].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(nOut->getId()), 0));
  
  ////Alternative topology 1:
  root = new bpp::Node(0);
  nOut = new bpp::Node(6, outgroup);
  n1 = new bpp::Node(1);
  n2 = new bpp::Node(2);
  n3 = new bpp::Node(3, species1);
  n4 = new bpp::Node(4, species2);
  n5 = new bpp::Node(5, species3);
  n2->addSon(n3);
  n2->addSon(n4);
  n1->addSon(n2);
  n1->addSon(n5);
  root->addSon(n1);
  root->addSon(nOut);
  trees_[1] = new bpp::TreeTemplate<bpp::Node>(root);
  trees_[1]->resetNodesId();
  brlenParameters_[1].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n1->getId()), 0));
  brlenParameters_[1].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n2->getId()), 0));
  brlenParameters_[1].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n3->getId()), 0));
  brlenParameters_[1].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n4->getId()), 0));
  brlenParameters_[1].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n5->getId()), 0));
  brlenParameters_[1].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(nOut->getId()), 0));
  
  ////Alternative topology 2:
  root = new bpp::Node(0);
  nOut = new bpp::Node(6, outgroup);
  n1 = new bpp::Node(1);
  n2 = new bpp::Node(2);
  n3 = new bpp::Node(3, species1);
  n4 = new bpp::Node(4, species3);
  n5 = new bpp::Node(5, species2);
  n2->addSon(n3);
  n2->addSon(n4);
  n1->addSon(n2);
  n1->addSon(n5);
  root->addSon(n1);
  root->addSon(nOut);
  trees_[2] = new bpp::TreeTemplate<bpp::Node>(root);
  trees_[2]->resetNodesId();
  brlenParameters_[2].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n1->getId()), 0));
  brlenParameters_[2].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n2->getId()), 0));
  brlenParameters_[2].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n3->getId()), 0));
  brlenParameters_[2].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n4->getId()), 0));
  brlenParameters_[2].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n5->getId()), 0));
  brlenParameters_[2].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(nOut->getId()), 0));

  ////Alternative topology 3:
  root = new bpp::Node(0);
  nOut = new bpp::Node(6, outgroup);
  n1 = new bpp::Node(1);
  n2 = new bpp::Node(2);
  n3 = new bpp::Node(3, species2);
  n4 = new bpp::Node(4, species3);
  n5 = new bpp::Node(5, species1);
  n2->addSon(n3);
  n2->addSon(n4);
  n1->addSon(n2);
  n1->addSon(n5);
  root->addSon(n1);
  root->addSon(nOut);
  trees_[3] = new bpp::TreeTemplate<bpp::Node>(root);
  trees_[3]->resetNodesId();
  brlenParameters_[3].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n1->getId()), 0));
  brlenParameters_[3].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n2->getId()), 0));
  brlenParameters_[3].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n3->getId()), 0));
  brlenParameters_[3].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n4->getId()), 0));
  brlenParameters_[3].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(n5->getId()), 0));
  brlenParameters_[3].addParameter(bpp::Parameter("BrLen" + bpp::TextTools::toString(nOut->getId()), 0));

  applyParameters();
}

void ThreeSpeciesCoalHmmStateAlphabet09::applyParameters()
{
  if (parametrization_ == 1) {
    tau1_   = getParameter_(0).getValue();
    tau2_   = getParameter_(1).getValue();
  } else {
    double sigma1 = getParameter_(0).getValue();
    double tau12  = getParameter_(1).getValue();
    tau1_ = sigma1 * tau12;
    tau2_ = tau12 - tau1_;
  }
  double theta1 = getParameter_(3).getValue();
  double theta2 = getParameter_(4).getValue();
  if (useMedian_)
  {
    a_  = tau1_ - theta1 * log(1. - (1. - exp(-tau2_ / theta1)) / 2.);
    b_  = tau1_ + tau2_ + 0.6931472 * theta2 - a_;
    a2_ = tau1_ + tau2_ + (0.6931472/3.) * theta2;
    b2_ = 0.6931472 * theta2;
  }
  else
  {
    a_  = tau1_ + theta1 - tau2_ * exp(-tau2_ / theta1) / (1. - exp(-tau2_ / theta1));
    b_  = tau1_ + tau2_ + theta2 - a_;
    a2_ = tau1_ + tau2_ + (1./3.) * theta2;
    b2_ = theta2;
  }
  c2_ = getParameter_(2).getValue();
  c_ = c2_ + theta2 / 3.;
  tau3_ = a2_ + b2_ + c2_ - tau1_ - tau2_;
  
  a_  = std::max(0.000001, a_);
  b_  = std::max(0.000001, b_);
  c_  = std::max(0.000001, c_);
  a2_ = std::max(0.000001, a2_);
  b2_ = std::max(0.000001, b2_);
  c2_ = std::max(0.000001, c2_);

  //Actualize trees:
  trees_[0]->getRootNode()->getSon(0)->setDistanceToFather(c_);
  trees_[0]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b_);
  trees_[0]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a_);
  trees_[0]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a_);
  trees_[0]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a_ + b_);
  trees_[0]->getRootNode()->getSon(1)->setDistanceToFather(a_ + b_ + c_);

  trees_[1]->getRootNode()->getSon(0)->setDistanceToFather(c2_);
  trees_[1]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b2_);
  trees_[1]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a2_);
  trees_[1]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a2_);
  trees_[1]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a2_ + b2_);
  trees_[1]->getRootNode()->getSon(1)->setDistanceToFather(a2_ + b2_ + c2_);

  trees_[2]->getRootNode()->getSon(0)->setDistanceToFather(c2_);
  trees_[2]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b2_);
  trees_[2]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a2_);
  trees_[2]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a2_);
  trees_[2]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a2_ + b2_);
  trees_[2]->getRootNode()->getSon(1)->setDistanceToFather(a2_ + b2_ + c2_);

  trees_[3]->getRootNode()->getSon(0)->setDistanceToFather(c2_);
  trees_[3]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b2_);
  trees_[3]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a2_);
  trees_[3]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a2_);
  trees_[3]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a2_ + b2_);
  trees_[3]->getRootNode()->getSon(1)->setDistanceToFather(a2_ + b2_ + c2_);

  //Actualize branch lengths parameters:
  brlenParameters_[0][0].setValue(c_);
  brlenParameters_[0][1].setValue(b_);
  brlenParameters_[0][2].setValue(a_);
  brlenParameters_[0][3].setValue(a_);
  brlenParameters_[0][4].setValue(a_ + b_);
  brlenParameters_[0][5].setValue(a_ + b_ + c_);

  brlenParameters_[1][0].setValue(c2_);
  brlenParameters_[1][1].setValue(b2_);
  brlenParameters_[1][2].setValue(a2_);
  brlenParameters_[1][3].setValue(a2_);
  brlenParameters_[1][4].setValue(a2_ + b2_);
  brlenParameters_[1][5].setValue(a2_ + b2_ + c2_);

  brlenParameters_[2][0].setValue(c2_);
  brlenParameters_[2][1].setValue(b2_);
  brlenParameters_[2][2].setValue(a2_);
  brlenParameters_[2][3].setValue(a2_);
  brlenParameters_[2][4].setValue(a2_ + b2_);
  brlenParameters_[2][5].setValue(a2_ + b2_ + c2_);
  
  brlenParameters_[3][0].setValue(c2_);
  brlenParameters_[3][1].setValue(b2_);
  brlenParameters_[3][2].setValue(a2_);
  brlenParameters_[3][3].setValue(a2_);
  brlenParameters_[3][4].setValue(a2_ + b2_);
  brlenParameters_[3][5].setValue(a2_ + b2_ + c2_);
}

void ThreeSpeciesCoalHmmStateAlphabet09::fireParameterChanged(const bpp::ParameterList& pl)
{
  applyParameters();
  std::vector<unsigned int> states(4);
  states[0] = 0;
  states[1] = 1;
  states[2] = 2;
  states[3] = 3;

  //Notify listeners:
  bpp::StateChangedEvent event(states);
  fireStateChangedEvent(event);
}

void ThreeSpeciesCoalHmmStateAlphabet09::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  double tau1   = getTau1();
  double tau2   = getTau2();
  double tau3   = getTau3();
  double theta1 = getTheta1();
  double theta2 = getTheta2();
  (out << "tau1 = " << tau1).endLine();
  (out << "tau2 = " << tau2).endLine();
  (out << "tau3 = " << tau3).endLine();
  (out << "theta1 = " << theta1).endLine();
  (out << "theta2 = " << theta2).endLine();
  (out << "a = " << a_).endLine();
  (out << "b = " << b_).endLine();
  (out << "c = " << c_).endLine();
  (out << "a2 = " << a2_).endLine();
  (out << "b2 = " << b2_).endLine();
  (out << "c2 = " << c2_).endLine();
}

