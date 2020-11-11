//
// File: ThreeSpeciesCoalHmmStateAlphabet07.cpp
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

#include "ThreeSpeciesCoalHmmStateAlphabet07.h"

#include <memory>
using namespace std;

ThreeSpeciesCoalHmmStateAlphabet07::ThreeSpeciesCoalHmmStateAlphabet07(
    const std::string& species1, const std::string& species2, const std::string& species3, const std::string& outgroup,
    double a, double b, double c, double a2, bool reparam) :
  AbstractParametrizable(""), brlenParameters_(), reparam_(reparam)
{
  double p = (a2 - a) / b;
  double c2 = c - 0.5 * b * (1. - p);
  addParameter_(new bpp::Parameter("coal.a", a, bpp::Parameter::R_PLUS_STAR));
  addParameter_(new bpp::Parameter("coal.b", b, bpp::Parameter::R_PLUS_STAR));
  if (reparam)
  {
    addParameter_(new bpp::Parameter("coal.c2", c2, bpp::Parameter::R_PLUS_STAR));
    addParameter_(new bpp::Parameter("coal.p", p, make_shared<bpp::IntervalConstraint>(0, 1, false, false))); //For estimation purposes, minimum 0.3??.
  }
  else
  {
    addParameter_(new bpp::Parameter("coal.c", c, bpp::Parameter::R_PLUS_STAR));
    addParameter_(new bpp::Parameter("coal.a2", p, bpp::Parameter::R_PLUS_STAR));
  }

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

void ThreeSpeciesCoalHmmStateAlphabet07::applyParameters()
{
  double a  = getParameterValue("coal.a");
  double b  = getParameterValue("coal.b");
  double c;
  double a2, b2, c2;
  if(reparam_)
  {
    c2 = getParameterValue("coal.c2");
    double p  = getParameterValue("coal.p");
    a2 = p * b + a;
    b2 = (3. / 2.) * (1. - p) * b;
    c  = c2 + 0.5 * b * (1. - p);
  }
  else
  {
    c  = getParameterValue("coal.c");
    a2 = getParameterValue("coal.a2");
    b2 = (3./2.) * (a + b -a2);
    c2 = a + b + c - a2 - b2;
  }

  a  = std::max(0.000001, a);
  b  = std::max(0.000001, b);
  c  = std::max(0.000001, c);
  a2  = std::max(0.000001, a2);
  b2  = std::max(0.000001, b2);
  c2  = std::max(0.000001, c2);

  //cout << (a + b + c - a2 - b2 - c2) << endl;
  //cout << "a=" << a << endl;
  //cout << "b=" << b << endl;
  //cout << "c=" << c << endl;
  //cout << "a2=" << a2 << endl;
  //cout << "b2=" << b2 << endl;
  //cout << "c2=" << c2 << endl;

  //Actualize trees:
  trees_[0]->getRootNode()->getSon(0)->setDistanceToFather(c);
  trees_[0]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b);
  trees_[0]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a);
  trees_[0]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a);
  trees_[0]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a + b);
  trees_[0]->getRootNode()->getSon(1)->setDistanceToFather(a + b + c);

  trees_[1]->getRootNode()->getSon(0)->setDistanceToFather(c2);
  trees_[1]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b2);
  trees_[1]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a2);
  trees_[1]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a2);
  trees_[1]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a2 + b2);
  trees_[1]->getRootNode()->getSon(1)->setDistanceToFather(a2 + b2 + c2);

  trees_[2]->getRootNode()->getSon(0)->setDistanceToFather(c2);
  trees_[2]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b2);
  trees_[2]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a2);
  trees_[2]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a2);
  trees_[2]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a2 + b2);
  trees_[2]->getRootNode()->getSon(1)->setDistanceToFather(a2 + b2 + c2);

  trees_[3]->getRootNode()->getSon(0)->setDistanceToFather(c2);
  trees_[3]->getRootNode()->getSon(0)->getSon(0)->setDistanceToFather(b2);
  trees_[3]->getRootNode()->getSon(0)->getSon(0)->getSon(0)->setDistanceToFather(a2);
  trees_[3]->getRootNode()->getSon(0)->getSon(0)->getSon(1)->setDistanceToFather(a2);
  trees_[3]->getRootNode()->getSon(0)->getSon(1)->setDistanceToFather(a2 + b2);
  trees_[3]->getRootNode()->getSon(1)->setDistanceToFather(a2 + b2 + c2);

  //Actualize branch lengths parameters:
  brlenParameters_[0][0].setValue(c);
  brlenParameters_[0][1].setValue(b);
  brlenParameters_[0][2].setValue(a);
  brlenParameters_[0][3].setValue(a);
  brlenParameters_[0][4].setValue(a + b);
  brlenParameters_[0][5].setValue(a + b + c);

  brlenParameters_[1][0].setValue(c2);
  brlenParameters_[1][1].setValue(b2);
  brlenParameters_[1][2].setValue(a2);
  brlenParameters_[1][3].setValue(a2);
  brlenParameters_[1][4].setValue(a2 + b2);
  brlenParameters_[1][5].setValue(a2 + b2 + c2);

  brlenParameters_[2][0].setValue(c2);
  brlenParameters_[2][1].setValue(b2);
  brlenParameters_[2][2].setValue(a2);
  brlenParameters_[2][3].setValue(a2);
  brlenParameters_[2][4].setValue(a2 + b2);
  brlenParameters_[2][5].setValue(a2 + b2 + c2);
  
  brlenParameters_[3][0].setValue(c2);
  brlenParameters_[3][1].setValue(b2);
  brlenParameters_[3][2].setValue(a2);
  brlenParameters_[3][3].setValue(a2);
  brlenParameters_[3][4].setValue(a2 + b2);
  brlenParameters_[3][5].setValue(a2 + b2 + c2);
}

void ThreeSpeciesCoalHmmStateAlphabet07::fireParameterChanged(const bpp::ParameterList& pl)
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

void ThreeSpeciesCoalHmmStateAlphabet07::printUserFriendlyParameters(bpp::OutputStream& out) const
{
  double a  = getParameterValue("coal.a");
  double b  = getParameterValue("coal.b");
  double c;
  double a2, b2, c2;
  if(reparam_)
  {
    c2 = getParameterValue("coal.c2");
    double p  = getParameterValue("coal.p");
    a2 = p * b + a;
    b2 = (3. / 2.) * (1. - p) * b;
    c  = c2 + 0.5 * b * (1. - p);
  }
  else
  {
    c  = getParameterValue("coal.c");
    a2 = getParameterValue("coal.a2");
    b2 = (3./2.) * (a + b -a2);
    c2 = a + b + c - a2 - b2;
  }
  (out << "coal.a = " << a).endLine();
  (out << "coal.b = " << b).endLine();
  (out << "coal.c = " << c).endLine();
  (out << "coal.a2 = " << a2).endLine();
  (out << "coal.b2 = " << b2).endLine();
  (out << "coal.c2 = " << c2).endLine();
}

