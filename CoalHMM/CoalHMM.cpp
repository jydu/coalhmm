//
// File: CoalHMM.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 26 09:44 2007
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
#include "ThreeSpeciesCoalHmmStateAlphabet09.h"
#include "TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet.h"
#include "TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet.h"
#include "NonClockAverageCoalHmmStateAlphabet.h"

#include "IlsCoalHmmTransitionMatrix07.h"
#include "IlsCoalHmmTransitionMatrix09.h"
#include "RateAndCoalHmmTransitionMatrix.h"
#include "TwoSpeciesDiscretizedCoalHmmTransitionMatrix.h"

#include "HomogeneousCoalHmmEmissionProbabilities.h"
#include "HmmTools.h"
#include "CoalHmmSimulator.h"

// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/NumCalcApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Hmm/RescaledHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/LowMemoryRescaledHmmLikelihood.h>
#include <Bpp/Numeric/Hmm/LogsumHmmLikelihood.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/CompressedVectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Io/Fasta.h>

// From bpp-phyl:
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/OptimizationTools.h>

using namespace bpp;

// From the STL:
#include <string>

using namespace std;

void help()
{
  cout << "coalhmm param=optionfile p1=v1 p2=v2 ..." << endl;
}

void printStates(const vector<int>& states, ostream& out)
{
  string codes = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  for(size_t i = 0; i < states.size(); i++)
  {
    switch(states[i])
    {
      case(-3):
        out << "?"; break;
      case(-10):
        out << "|"; break;
      default:
        out << codes[static_cast<size_t>(states[i])];
    }
  }
  out << endl;
}



AverageCoalHmmStateAlphabet* getThreeSpeciesAverageCoalHiddenAlphabet(map<string, string>& params)
{
  string species1 = ApplicationTools::getStringParameter("species1" , params, "Ingroup1");
  string species2 = ApplicationTools::getStringParameter("species2" , params, "Ingroup2");
  string species3 = ApplicationTools::getStringParameter("species3" , params, "Ingroup3");
  string outgroup = ApplicationTools::getStringParameter("outgroup" , params, "Outgroup");
  
  string impl = ApplicationTools::getStringParameter("implementation", params, "09", "", false, false);
  if (impl == "07")
  {
    double a  = ApplicationTools::getDoubleParameter("a"  , params, 2);
    double b  = ApplicationTools::getDoubleParameter("b"  , params, 2);
    double c  = ApplicationTools::getDoubleParameter("c"  , params, 2);
    double a2 = ApplicationTools::getDoubleParameter("a2", params, 1);
    ApplicationTools::displayResult("coal.a" , a);
    ApplicationTools::displayResult("coal.b" , b);
    ApplicationTools::displayResult("coal.c" , c);
    ApplicationTools::displayResult("coal.a2", a2);
    ApplicationTools::displayResult("Genealogies", string("3 species \'07"));
    return new ThreeSpeciesCoalHmmStateAlphabet07(species1, species2, species3, outgroup, a, b, c, a2, true);
  }
  else if (impl == "07nr")
  {
    double a  = ApplicationTools::getDoubleParameter("a"  , params, 2);
    double b  = ApplicationTools::getDoubleParameter("b"  , params, 2);
    double c  = ApplicationTools::getDoubleParameter("c"  , params, 2);
    double a2 = ApplicationTools::getDoubleParameter("a2", params, 1);
    ApplicationTools::displayResult("coal.a" , a);
    ApplicationTools::displayResult("coal.b" , b);
    ApplicationTools::displayResult("coal.c" , c);
    ApplicationTools::displayResult("coal.a2", a2);
    return new ThreeSpeciesCoalHmmStateAlphabet07(species1, species2, species3, outgroup, a, b, c, a2, false);
  }
  else if (impl == "09")
  {
    ApplicationTools::displayResult("Genealogies", string("3 species \'09"));
    double tau1   = ApplicationTools::getDoubleParameter("tau1"  , params, 2);
    double tau2   = ApplicationTools::getDoubleParameter("tau2"  , params, 2);
    double c2     = ApplicationTools::getDoubleParameter("c2"  , params, 2);
    double theta1 = ApplicationTools::getDoubleParameter("theta1", params, 1);
    double theta2 = ApplicationTools::getDoubleParameter("theta2", params, 1);
    bool useMedian = ApplicationTools::getBooleanParameter("median", params, false);
    double tau_min   = ApplicationTools::getDoubleParameter("tau.min", params, 0.0, "", true, false);
    double tau_max   = ApplicationTools::getDoubleParameter("tau.max", params, 1, "", true, false);
    double theta_min = ApplicationTools::getDoubleParameter("theta.min", params, 0.0001, "", true, false);
    double theta_max = ApplicationTools::getDoubleParameter("theta.max", params, 1, "", true, false);
    short parametrization = ApplicationTools::getParameter<short>("parametrization", params, 1);
    if (parametrization == 1)
      ApplicationTools::displayResult("Speciation times parametrized as", string("tau1 + tau2"));
    else
      ApplicationTools::displayResult("Speciation times parametrized as", string("sigma1 + tau12"));
    AverageCoalHmmStateAlphabet* hAlpha = new ThreeSpeciesCoalHmmStateAlphabet09(species1, species2, species3, outgroup, tau1, tau2, c2, theta1, theta2, parametrization, useMedian, tau_min, theta_min, tau_max, theta_max);
    if (parametrization == 1) {
      ApplicationTools::displayResult("coal.tau1", TextTools::toString(tau1) + " " + hAlpha->getParameter("tau1").getConstraint()->getDescription());
      ApplicationTools::displayResult("coal.tau2", TextTools::toString(tau2) + " " + hAlpha->getParameter("tau2").getConstraint()->getDescription());
    } else {
      ApplicationTools::displayResult("coal.sigma1", TextTools::toString(hAlpha->getParameter("sigma1").getValue()) + " " + hAlpha->getParameter("sigma1").getConstraint()->getDescription());
      ApplicationTools::displayResult("coal.tau12" , TextTools::toString(hAlpha->getParameter("tau12").getValue())  + " " + hAlpha->getParameter("tau12").getConstraint()->getDescription());
    }
    ApplicationTools::displayResult("coal.c2"  , TextTools::toString(c2)   + " " + hAlpha->getParameter("c2").getConstraint()->getDescription());
    ApplicationTools::displayResult("coal.theta1", TextTools::toString(theta1) + " " + hAlpha->getParameter("theta1").getConstraint()->getDescription());
    ApplicationTools::displayResult("coal.theta2", TextTools::toString(theta2) + " " + hAlpha->getParameter("theta2").getConstraint()->getDescription());
    return hAlpha;
  }
  else throw Exception("Bad implementation specified: " + impl);
}



AverageCoalHmmStateAlphabet* getTwoSpeciesAverageCoalHiddenAlphabet(map<string, string>& params)
{
  string species1 = ApplicationTools::getStringParameter("species1" , params, "Ingroup1");
  string species2 = ApplicationTools::getStringParameter("species2" , params, "Ingroup2");

  bool withOutgroup = params.find("outgroup") != params.end() && !TextTools::isEmpty(params["outgroup"]);
  size_t nbClasses = ApplicationTools::getParameter<size_t>("numberOfClasses", params, 3, "");
  ApplicationTools::displayResult("Number of hidden states", nbClasses);
  if (withOutgroup)
  {
    string outgroup = ApplicationTools::getStringParameter("outgroup" , params, "Outgroup");
    double tau1     = ApplicationTools::getDoubleParameter("tau1"   , params, 0.001);
    double tau2     = ApplicationTools::getDoubleParameter("tau2"   , params, 0.01);
    double theta12  = ApplicationTools::getDoubleParameter("theta12", params, 0.001);
    bool useMedian  = ApplicationTools::getBooleanParameter("median", params, false);
    double tau1_min  = ApplicationTools::getDoubleParameter("tau1.min", params, 0.0, "", true, false);
    double tau2_min  = ApplicationTools::getDoubleParameter("tau2.min", params, 0.0001, "", true, false);
    double theta_min = ApplicationTools::getDoubleParameter("theta.min", params, 0.0001, "", true, false);
    TwoSpeciesDiscretizedCoalHmmStateAlphabet* hAlpha = new TwoSpeciesWithOutgroupDiscretizedCoalHmmStateAlphabet(species1, species2, outgroup, tau1, tau2, theta12, nbClasses, useMedian, tau1_min, tau2_min, theta_min);
    ApplicationTools::displayResult("coal.tau1", TextTools::toString(tau1) + " " + hAlpha->getParameter("tau1").getConstraint()->getDescription());
    ApplicationTools::displayResult("coal.tau2", TextTools::toString(tau2) + " " + hAlpha->getParameter("tau2").getConstraint()->getDescription());
    ApplicationTools::displayResult("coal.theta12", TextTools::toString(theta12) + " " + hAlpha->getParameter("theta12").getConstraint()->getDescription());
    return hAlpha;
  }
  else
  {
    double tau1      = ApplicationTools::getDoubleParameter("tau1"   , params, 0.001);
    double theta12   = ApplicationTools::getDoubleParameter("theta12", params, 0.001);
    bool useMedian   = ApplicationTools::getBooleanParameter("median", params, false);
    double tau_min   = ApplicationTools::getDoubleParameter("tau.min", params, 0.0, "", true, false);
    double theta_min = ApplicationTools::getDoubleParameter("theta.min", params, 0.0001, "", true, false);
    TwoSpeciesDiscretizedCoalHmmStateAlphabet* hAlpha = new TwoSpeciesWithoutOutgroupDiscretizedCoalHmmStateAlphabet(species1, species2, tau1, theta12, nbClasses, useMedian, tau_min, theta_min);
    ApplicationTools::displayResult("coal.tau1", TextTools::toString(tau1) + " " + hAlpha->getParameter("tau1").getConstraint()->getDescription());
    ApplicationTools::displayResult("coal.theta12", TextTools::toString(theta12) + " " + hAlpha->getParameter("theta12").getConstraint()->getDescription());
    return hAlpha;
  }
}

CoalHmmTransitionMatrix* getHiddenModel(const string& model, map<string, string>& params, CoalHmmStateAlphabet* hAlpha)
{
  CoalHmmTransitionMatrix* hModel = 0;
  if (model == "ILS")
  {
    string impl = ApplicationTools::getStringParameter("implementation", params, "09", "", false, false);
    if (impl == "07")
    {
      double s = ApplicationTools::getDoubleParameter("s", params, 0.1);
      double u = ApplicationTools::getDoubleParameter("u", params, 0.1);
      if (params.find("v") != params.end())
      {
        double v = ApplicationTools::getDoubleParameter("v", params, 0.1);
        ApplicationTools::displayResult("coal.s", s);
        ApplicationTools::displayResult("coal.u", u);
        ApplicationTools::displayResult("coal.v", v);
        hModel = new IlsCoalHmmTransitionMatrix07(hAlpha, s, u, v, v, true);
      }
      else
      {
        double v1 = ApplicationTools::getDoubleParameter("v1", params, 0.1);
        double v2 = ApplicationTools::getDoubleParameter("v2", params, 0.1);
        ApplicationTools::displayResult("coal.s", s);
        ApplicationTools::displayResult("coal.u", u);
        ApplicationTools::displayResult("coal.v1", v1);
        ApplicationTools::displayResult("coal.v2", v2);
        hModel = new IlsCoalHmmTransitionMatrix07(hAlpha, s, u, v1, v2, false);
      }
    }
    else if (impl == "09" || impl == "1rho")
    {
      double rho = ApplicationTools::getDoubleParameter("rho", params, 0.01);
      double rho_min = ApplicationTools::getDoubleParameter("rho.min", params, 0.0, "", true, false);
      double rho_max = ApplicationTools::getDoubleParameter("rho.max", params, 1000, "", true, false);
      hModel = new IlsCoalHmmTransitionMatrix09(dynamic_cast<ThreeSpeciesCoalHmmStateAlphabet*>(hAlpha), rho, rho_min, rho_max);
      ApplicationTools::displayResult("coal.rho", TextTools::toString(rho) + " " + hModel->getParameter("rho").getConstraint()->getDescription());
    }
    else if (impl == "2rhos12")
    {
      double rho12 = ApplicationTools::getDoubleParameter("rho12", params, 0.01);
      double rho3  = ApplicationTools::getDoubleParameter("rho3", params, 0.01);
      double rho_min = ApplicationTools::getDoubleParameter("rho.min", params, 0.0, "", true, false);
      double rho_max = ApplicationTools::getDoubleParameter("rho.max", params, 1000, "", true, false);
      hModel = new IlsCoalHmmTransitionMatrix09(dynamic_cast<ThreeSpeciesCoalHmmStateAlphabet*>(hAlpha), rho12, rho3, IlsCoalHmmTransitionMatrix09::TWO_RHOS_12, rho_min, rho_max);
      ApplicationTools::displayResult("coal.rho12", TextTools::toString(rho12) + " " + hModel->getParameter("rho12").getConstraint()->getDescription());
      ApplicationTools::displayResult("coal.rho3", TextTools::toString(rho3) + " " + hModel->getParameter("rho3").getConstraint()->getDescription());
    }
    else if (impl == "3rhos")
    {
      double rho1 = ApplicationTools::getDoubleParameter("rho1", params, 0.01);
      double rho2 = ApplicationTools::getDoubleParameter("rho2", params, 0.01);
      double rho3 = ApplicationTools::getDoubleParameter("rho3", params, 0.01);
      double rho_min = ApplicationTools::getDoubleParameter("rho.min", params, 0.0, "", true, false);
      double rho_max = ApplicationTools::getDoubleParameter("rho.max", params, 1000, "", true, false);
      hModel = new IlsCoalHmmTransitionMatrix09(dynamic_cast<ThreeSpeciesCoalHmmStateAlphabet*>(hAlpha), rho1, rho2, rho3, rho_min, rho_max);
      ApplicationTools::displayResult("coal.rho1", TextTools::toString(rho1) + " " + hModel->getParameter("rho1").getConstraint()->getDescription());
      ApplicationTools::displayResult("coal.rho2", TextTools::toString(rho2) + " " + hModel->getParameter("rho2").getConstraint()->getDescription());
      ApplicationTools::displayResult("coal.rho3", TextTools::toString(rho3) + " " + hModel->getParameter("rho3").getConstraint()->getDescription());
    }
    else throw Exception("Bad implementation specified: " + impl);
  }
  else if (model == "Divergence")
  {
    size_t nbsp = ApplicationTools::getParameter<size_t>("numberOfSpecies", params, 2);
    if (nbsp == 2)
    {
      double rho = ApplicationTools::getDoubleParameter("rho", params, 0.1);
      double rho_min = ApplicationTools::getDoubleParameter("rho.min", params, 0.0, "", true, false);
      bool old = ApplicationTools::getBooleanParameter("old", params, false, "", true, false);
      if( (params.find("theta1") != params.end() && params["theta1"] == "theta12")
       || (params.find("theta2") != params.end() && params["theta2"] == "theta12"))
      {
        if (old)
          hModel = new OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(dynamic_cast<TwoSpeciesDiscretizedCoalHmmStateAlphabet *>(hAlpha), rho, rho_min);
        else  
          hModel = new TwoSpeciesDiscretizedCoalHmmTransitionMatrix(dynamic_cast<TwoSpeciesDiscretizedCoalHmmStateAlphabet *>(hAlpha), rho, rho_min);
      }
      else
      {
        double theta1 = ApplicationTools::getDoubleParameter("theta1", params, 0.002);
        double theta2 = ApplicationTools::getDoubleParameter("theta2", params, 0.002);
        if (old)
          hModel = new OldTwoSpeciesDiscretizedCoalHmmTransitionMatrix(dynamic_cast<TwoSpeciesDiscretizedCoalHmmStateAlphabet *>(hAlpha), rho, theta1, theta2);
        else  
          hModel = new TwoSpeciesDiscretizedCoalHmmTransitionMatrix(dynamic_cast<TwoSpeciesDiscretizedCoalHmmStateAlphabet *>(hAlpha), rho, theta1, theta2);
        ApplicationTools::displayResult("coal.theta1", TextTools::toString(theta1) + " " + hModel->getParameter("theta1").getConstraint()->getDescription());
        ApplicationTools::displayResult("coal.theta2", TextTools::toString(theta2) + " " + hModel->getParameter("theta2").getConstraint()->getDescription());
      }
      ApplicationTools::displayResult("coal.rho", TextTools::toString(rho) + " " + TextTools::toString(hModel->getParameter("rho").getConstraint()->getDescription()));
    }
    else throw Exception("Sorry, the divergence model is only available for 2 species right now!");
  }
  else throw Exception("Unknown hidden substitution model: " + model);

  return hModel;
}



string getModel(const string& method, map<string, string>& params, AverageCoalHmmStateAlphabet*& hAlpha, CoalHmmTransitionMatrix*& hModel, DiscreteDistribution*& rDist)
{
  map<string, string> args;
  string model;
  KeyvalTools::parseProcedure(method, model, args);
	
  if (model == "ILS")
  {
    size_t nbsp = ApplicationTools::getParameter<size_t>("numberOfSpecies", args, 3);
    if (nbsp == 3)
    {
      hAlpha = getThreeSpeciesAverageCoalHiddenAlphabet(args);
    }
    else throw Exception("Only model with 3 species + outgroup is supported for now.");
    hModel = getHiddenModel(model, args, hAlpha);

    rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    return model;
  }
  else if(model == "Divergence")
  {
    size_t nbsp = ApplicationTools::getParameter<size_t>("numberOfSpecies", args, 2);
    if (nbsp == 2)
    {
      hAlpha = getTwoSpeciesAverageCoalHiddenAlphabet(args);
    }
    else throw Exception("Only model with 2 species +/- outgroup is supported for now.");
    hModel = getHiddenModel(model, args, hAlpha);

    rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    return model;
  }
  else if (model == "MM")
  {
    string subproc = ApplicationTools::getStringParameter("model", args, "none");
    AverageCoalHmmStateAlphabet* mhAlpha;
    CoalHmmTransitionMatrix* mhModel;
    DiscreteDistribution* mmDist;
    string subtype = getModel(subproc, params, mhAlpha, mhModel, mmDist);

    double w = ApplicationTools::getDoubleParameter("omega", args, 0.0001);
    ApplicationTools::displayResult("omega", TextTools::toString(w));
    hAlpha = new RateAndCoalHmmStateAlphabet(mhAlpha, mmDist);
    hModel = new RateAndCoalHmmTransitionMatrix(dynamic_cast<RateAndCoalHmmStateAlphabet*>(hAlpha), mhModel, w);
    rDist  = new ConstantDistribution(1.);
    return model + "+" + subtype;
  }
  else if (model == "Unclock")
  {
    string subproc = ApplicationTools::getStringParameter("model", args, "none");
    AverageCoalHmmStateAlphabet* chAlpha;
    string subtype = getModel(subproc, params, chAlpha, hModel, rDist);
 
    double error = ApplicationTools::getDoubleParameter("error", args, 0.);
    string species = ApplicationTools::getStringParameter("species", args, "");
    ApplicationTools::displayResult(species + "_error", TextTools::toString(error));
    hAlpha = new NonClockAverageCoalHmmStateAlphabet(chAlpha, species, error);
    return model + "+" + subtype;
  }
  else
  {
    throw Exception("Unknown coal-HMM method: " + method);
  }
}


class IterationCounter: public OptimizationListener
{
  protected:
    size_t _counter;
    size_t _maxIt;

  public:
    IterationCounter(size_t maxit): _counter(0), _maxIt(maxit) {}
  
  public:
    void optimizationInitializationPerformed(const OptimizationEvent &event) {}
    void optimizationStepPerformed(const OptimizationEvent &event)
    {
      _counter++;
      if(_counter >= _maxIt) throw Exception("Maximum number of optimization reached in optimization.");
    }
    bool listenerModifiesParameters () const { return false; }
    size_t getCount() const { return _counter; }
};



int main(int argc, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                    CoalHMM, version 1.0.3                      *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. 04/10/17 *" << endl;
  cout << "*          T. Mailund                                            *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (argc == 1)
  {
    help();
    return(0);
  }
  
  try {
 
  BppApplication coalhmm(argc, argv, "CoalHMM");
  coalhmm.startTimer();

  ApplicationTools::startTimer();

  map<string, string> params = coalhmm.getParams();


  /////////////////////////////////////////////////////////
  // Emission model:
  /////////////////////////////////////////////////////////
  
  Alphabet* alpha = SequenceApplicationTools::getAlphabet(params, "", false);
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "        ..:: Sequence model ::..").endLine();
  (*ApplicationTools::message).endLine();
  SubstitutionModel* sModel = PhylogeneticsApplicationTools::getSubstitutionModel(alpha, 0, 0, params);
  DiscreteDistribution* rDist = 0;
  AverageCoalHmmStateAlphabet* hAlpha = 0;
  CoalHmmTransitionMatrix* hModel = 0;
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "       ..:: Genealogy model ::..").endLine();
  (*ApplicationTools::message).endLine();
  string method = ApplicationTools::getStringParameter("coalmethod", params, "none");
  string modelType = getModel(method, params, hAlpha, hModel, rDist);
  ApplicationTools::displayResult("CoalHMM method:", modelType);
  
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "           ..:: Analysis ::..").endLine();
  (*ApplicationTools::message).endLine();
  string analysis = ApplicationTools::getStringParameter("analysis", params, "estimate");
  ApplicationTools::displayResult("Running mode", analysis);
  if (analysis == "estimate")
  {
    HomogeneousCoalHmmEmissionProbabilities* hmmEmissions = new HomogeneousCoalHmmEmissionProbabilities(hAlpha, sModel, rDist, true);

    bool multiParts = ApplicationTools::getBooleanParameter("input.sequence.multiparts", params, false);
    string prefix = ApplicationTools::getStringParameter("input.sequence.multiparts.prefix", params, "", "", true, false);
    SiteContainer* alignment = 0;
    vector<size_t> breakPoints;
    bool resetChain = false;
    if (multiParts)
    {
      resetChain = ApplicationTools::getBooleanParameter("input.sequence.multiparts.reset", params, true);
      ApplicationTools::displayResult("Reset chain for each chunk", (resetChain ? "yes" : "no"));
      string listPath = ApplicationTools::getAFilePath("input.sequence.list", params, true, true);
      size_t nbAli = 0;
      size_t totNbSites = 0;
      ifstream listFile(listPath.c_str(), ios::in);
      while (!listFile.eof())
      {
        string fileName = FileTools::getNextLine(listFile);
        if (TextTools::isEmpty(fileName)) continue;
        fileName = prefix + fileName;
        params["input.sequence.file"] = fileName;
        unique_ptr<VectorSiteContainer> tmp(SequenceApplicationTools::getSiteContainer(alpha, params, "", false, false));
        unique_ptr<SiteContainer> tmp2(SequenceApplicationTools::getSitesToAnalyse(*tmp, params));
        unique_ptr<SiteContainer> sites(PatternTools::getSequenceSubset(*tmp2, hAlpha->getSpeciesNames()));
        size_t nbSites = sites->getNumberOfSites();
        if (nbSites == 0)
        {
          ApplicationTools::displayMessage("Part " + fileName + " is empty and was ignored.");
        }
        else
        {
          string dname = fileName;
          if (dname.size() > 28)
          {
            dname = "[...]" + TextTools::resizeLeft(dname, 22);
          }
          ApplicationTools::displayResult(string("Adding part ") + dname, TextTools::toString(nbSites) + string(" sites found."));
          if (!alignment) alignment = new CompressedVectorSiteContainer(*sites);
          else
          {
            SiteContainerTools::merge(*alignment, *sites);
            breakPoints.push_back(totNbSites);
          }
          nbAli++;
          totNbSites += nbSites;
        }
      }
      ApplicationTools::displayResult("Number of parts", nbAli);
      ApplicationTools::displayResult("Total number of sites", totNbSites);
    }
    else
    {
      unique_ptr<VectorSiteContainer> tmp(SequenceApplicationTools::getSiteContainer(alpha, params));
      unique_ptr<SiteContainer> tmp2(PatternTools::getSequenceSubset(*tmp, hAlpha->getSpeciesNames()));
      alignment = SequenceApplicationTools::getSitesToAnalyse(*tmp2, params);
      size_t nbSites = alignment->getNumberOfSites();
      ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
    }

    hmmEmissions->setData(*alignment);

    hmmEmissions->initialize();
    
    string hmmalgo = ApplicationTools::getStringParameter("analysis.method", params, "rescaled");
    ApplicationTools::displayResult("HMM numerical method", hmmalgo);
    HmmLikelihood* hmmLik = 0;
    if (hmmalgo == "rescaled")
      hmmLik = new RescaledHmmLikelihood(hAlpha, hModel, hmmEmissions, "");
    else if(hmmalgo == "logsum")
      hmmLik = new LogsumHmmLikelihood(hAlpha, hModel, hmmEmissions, true, "");
    else if(hmmalgo == "lowmem")
      hmmLik = new LowMemoryRescaledHmmLikelihood(hAlpha, hModel, hmmEmissions, "");
    else throw Exception("Unknow HMM method, should be 'rescaled', 'logsum' or 'lowmem'.");

    Sequence* annotation = 0;
    string annotationPath = ApplicationTools::getAFilePath("input.sequence.annotation.file", params, false, true);
    ApplicationTools::displayResult("Annotation file", annotationPath);
    if (annotationPath != "none")
    {
      Fasta reader;
      OrderedSequenceContainer* cont = reader.readSequences(annotationPath, &AlphabetTools::DEFAULT_ALPHABET);
      annotation = new BasicSequence(cont->getSequence(0));
      delete cont;
    }

    /////////////////////////////////////////////////////////
    // Optimize parameters:
    /////////////////////////////////////////////////////////
 
    bool optimize = ApplicationTools::getBooleanParameter("optimize", params, false);
    IterationCounter* itC = 0;
    if (optimize)
    {
      ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-hmmLik->getValue(), 20));
      unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2);
  
      string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false);
      shared_ptr<OutputStream> messageHandler = 
        (mhPath == "none") ? nullptr :
          (mhPath == "std") ? ApplicationTools::message :
            shared_ptr<StlOutputStream>(new StlOutputStream(new ofstream(mhPath.c_str(), ios::out)));
      ApplicationTools::displayResult("Message handler", mhPath);

      string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false);
      shared_ptr<OutputStream> profiler = 
        (prPath == "none") ? nullptr :
          (prPath == "std") ? ApplicationTools::message :
             shared_ptr<StlOutputStream>(new StlOutputStream(new ofstream(prPath.c_str(), ios::out)));
      if (profiler) profiler->setPrecision(20);
      ApplicationTools::displayResult("Profiler", prPath);

      vector<string> hAlphabetParameterNames = hmmLik->getHmmStateAlphabet().getParameters().getParameterNames();
      vector<string> hModelParameterNames    = hmmLik->getHmmTransitionMatrix().getParameters().getParameterNames();
      vector<string> emissionsParameterNames = hmmLik->getHmmEmissionProbabilities().getParameters().getParameterNames();
      ParameterList pl = hmmLik->getParameters();
  
      // Should I ignore some parameters?
      string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", "", true, false);
      StringTokenizer st(paramListDesc, ",");
      while (st.hasMoreToken())
      {
        string param = st.nextToken();
        try
        {
          pl.deleteParameter(param);
        }
        catch (ParameterNotFoundException& pnfe)
        {
          ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
        }
      }
  
      size_t maxIt = ApplicationTools::getParameter<size_t>("optimization.max_number_iterations", params, 1000);
      itC = new IterationCounter(maxIt); 
      ApplicationTools::displayResult("Maximum number of iterations", maxIt);
    
      unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000);
      ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
  
      double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001);
      ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

      bool preOpt = ApplicationTools::getBooleanParameter("optimization.pre", params, false);
      ApplicationTools::displayResult("Pre-optimization of mutation parameters", (preOpt ? "yes" : "no"));
      if (preOpt)
      {
        TreeTemplate<Node> tree(dynamic_cast<AverageCoalHmmStateAlphabet&>(hmmLik->getHmmStateAlphabet()).getState(0));
        tree.unroot();
        RHomogeneousTreeLikelihood tl(tree, *alignment, sModel, rDist, true, false, true);
        tl.initialize();
        // Should I ignore some parameters?
        string paramListDescPre = ApplicationTools::getStringParameter("optimization.pre.ignore_parameter", params, "", "", true, false);
        StringTokenizer stPre(paramListDescPre, ",");
        ParameterList parameters = tl.getParameters();
        while (stPre.hasMoreToken())
        {
          string param = stPre.nextToken();
          try
          {
            parameters.deleteParameter(param);
          }
          catch (ParameterNotFoundException& pnfe)
          {
            ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
          }
        }
        OptimizationTools::optimizeNumericalParameters2(&tl, parameters, 0, tolerance, nbEvalMax, messageHandler.get(), profiler.get(), false, false, 1);
        ParameterList pli = tl.getSubstitutionModelParameters();
        pli.addParameters(tl.getRateDistributionParameters());
        hmmLik->setParameters(pli);
        ApplicationTools::displayResult("New log likelihood", TextTools::toString(-hmmLik->getValue(), 20));
      }
      delete alignment;
      
      string optMethod = ApplicationTools::getStringParameter("optimization.method", params, "fullD");
      ApplicationTools::displayResult("Optimization method", optMethod);
      Optimizer* optimizer = 0;
    
      bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false);
      ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));
      Function* f = hmmLik;
      if (reparam)
      {
        f = new ReparametrizationFunctionWrapper(hmmLik);

        //Reset parameters to remove constraints:
        pl = f->getParameters().createSubList(pl.getParameterNames());
      }
      if (optMethod == "fullD")
      {
        ThreePointsNumericalDerivative* fun = new ThreePointsNumericalDerivative(f);
        vector<string> names = pl.getParameterNames();
        fun->setParametersToDerivate(names);
        optimizer = new PseudoNewtonOptimizer(fun);
        dynamic_cast<PseudoNewtonOptimizer*>(optimizer)->disableCG(); //Conjugate gradient does not perform very well with coalhmm likelihood function as it seems...
        NaNListener nanL(optimizer, hmmLik); 
        optimizer->addOptimizationListener(itC);
        optimizer->addOptimizationListener(&nanL);
        optimizer->setVerbose(optVerbose);
        optimizer->setProfiler(profiler.get());
        optimizer->setMessageHandler(messageHandler.get());
        optimizer->setMaximumNumberOfEvaluations(nbEvalMax);
        optimizer->getStopCondition()->setTolerance(tolerance);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        optimizer->init(pl);
        try { optimizer->optimize(); } catch (Exception & ex) { cerr << ex.what() << endl; }
        delete optimizer;
        delete fun;
        ApplicationTools::displayTaskDone();
      }
      else if (optMethod == "splitD")
      {
        unsigned int nbSteps = ApplicationTools::getParameter<unsigned int>("optimization.number_of_steps", params, 10);
        ApplicationTools::displayResult("Number of steps", TextTools::toString(nbSteps));
 
        ThreePointsNumericalDerivative* fun = new ThreePointsNumericalDerivative(f);
        fun->setParametersToDerivate(pl.getParameterNames());

        MetaOptimizerInfos* metaInf = new MetaOptimizerInfos();
      
        Optimizer *opt1 = new PseudoNewtonOptimizer(fun);
        metaInf->addOptimizer("Emission parameters", opt1, emissionsParameterNames, 2, MetaOptimizerInfos::IT_TYPE_FULL);
    
        Optimizer *opt2 = new PseudoNewtonOptimizer(fun);
        metaInf->addOptimizer("Hidden alphabet parameters", opt2, hAlphabetParameterNames, 2, MetaOptimizerInfos::IT_TYPE_FULL);
    
        Optimizer *opt3 = new PseudoNewtonOptimizer(fun);
        metaInf->addOptimizer("Hidden substitution parameters", opt3, hModelParameterNames, 2, MetaOptimizerInfos::IT_TYPE_FULL);

        optimizer = new MetaOptimizer(fun, metaInf, nbSteps);
        optimizer->addOptimizationListener(itC);
        optimizer->setVerbose(optVerbose);
        optimizer->setProfiler(profiler.get());
        optimizer->setMessageHandler(messageHandler.get());
        optimizer->setMaximumNumberOfEvaluations(nbEvalMax);
        optimizer->getStopCondition()->setTolerance(tolerance);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        optimizer->init(pl);
        optimizer->optimize();
        delete optimizer;
        delete fun;
        ApplicationTools::displayTaskDone();
      }
      else if (optMethod == "DB")
      {
        unsigned int nbSteps = ApplicationTools::getParameter<unsigned int>("optimization.number_of_steps", params, 10);
        ApplicationTools::displayResult("Number of steps", TextTools::toString(nbSteps));
 
        ThreePointsNumericalDerivative* fun = new ThreePointsNumericalDerivative(f);
        fun->setParametersToDerivate(hAlphabetParameterNames);

        MetaOptimizerInfos* metaInf = new MetaOptimizerInfos();
      
        Optimizer *opt1 = new SimpleMultiDimensions(f);
        metaInf->addOptimizer("Emission parameters", opt1, emissionsParameterNames, 0, MetaOptimizerInfos::IT_TYPE_STEP);
    
        Optimizer *opt2 = new PseudoNewtonOptimizer(fun);
        metaInf->addOptimizer("Hidden alphabet parameters", opt2, hAlphabetParameterNames, 2, MetaOptimizerInfos::IT_TYPE_FULL);
    
        Optimizer *opt3 = new SimpleMultiDimensions(f);
        metaInf->addOptimizer("Hidden substitution parameters", opt3, hModelParameterNames, 0, MetaOptimizerInfos::IT_TYPE_STEP);

        optimizer = new MetaOptimizer(fun, metaInf, nbSteps);
        optimizer->addOptimizationListener(itC);
        optimizer->setVerbose(optVerbose);
        optimizer->setProfiler(profiler.get());
        optimizer->setMessageHandler(messageHandler.get());
        optimizer->setMaximumNumberOfEvaluations(nbEvalMax);
        optimizer->getStopCondition()->setTolerance(tolerance);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        optimizer->init(pl);
        optimizer->optimize();
        delete optimizer;
        delete fun;
        ApplicationTools::displayTaskDone();
      }
      else throw Exception("Unknown optimization method: " + optMethod);

      //Post optimization (in case of numerical derivative lack of precision):
      bool final = ApplicationTools::getBooleanParameter("optimization.final", params, true);
      if (final)
      {
        ApplicationTools::displayTask("\n\nFinal optimization step:");
        pl.matchParametersValues(f->getParameters());
        optimizer = new PowellMultiDimensions(f);
        optimizer->setVerbose(optVerbose);
        optimizer->setProfiler(profiler.get());
        optimizer->setMessageHandler(messageHandler.get());
        optimizer->setMaximumNumberOfEvaluations(nbEvalMax);
        optimizer->getStopCondition()->setTolerance(tolerance);
        optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
        try //In case sthg happens...
        {
          optimizer->init(pl);
        } 
        catch (ConstraintException & e)
        {
          pl.printParameters(cout);
          throw(e);
        }
        optimizer->optimize();
        delete optimizer;
        ApplicationTools::displayTaskDone();
      }
      if (f != hmmLik) delete f;

      //Print estimated parameters to file:
      string estPath = ApplicationTools::getAFilePath("output.estimated.parameters", params, false, false);
      ApplicationTools::displayResult("Estimated parameters file", estPath);
      if (estPath != "none")
      {
        StlOutputStream out(new ofstream(estPath.c_str(), ios::out));
        //hmmLik->printUserFriendlyParameters(out);
        ParameterList plEst;
        plEst = hAlpha->getParameters();
        for (size_t i = 0; i < plEst.size(); ++i)
          (out << plEst[i].getName() << " = " << plEst[i].getValue()).endLine();
        plEst = hModel->getParameters();
        for (size_t i = 0; i < plEst.size(); ++i)
          (out << plEst[i].getName() << " = " << plEst[i].getValue()).endLine();
        PhylogeneticsApplicationTools::printParameters(sModel, out);
        PhylogeneticsApplicationTools::printParameters(rDist, out);
      } 
    }
    else
    {
      delete alignment;
    }
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-hmmLik->getValue(), 20));

    /////////////////////////////////////////////////////////
    //Print user-friendly parameters:
    /////////////////////////////////////////////////////////
    string ufpPath = ApplicationTools::getAFilePath("output.userfriendly.parameters", params, false, false);
    ApplicationTools::displayResult("User-friendly parameters file", ufpPath);
    if (ufpPath != "none")
    {
      shared_ptr<OutputStream> out;
      if (ufpPath == "std")
      {
        out = ApplicationTools::message;
      }
      else
      {
        out = shared_ptr<StlOutputStream>(new StlOutputStream(new ofstream(ufpPath.c_str(), ios::out)));
      }

      hAlpha->printUserFriendlyParameters(*out);
      hModel->printUserFriendlyParameters(*out);
      hmmEmissions->printUserFriendlyParameters(*out);
      (*out << "logL = ").setPrecision(20);
      (*out << - hmmLik->getValue()).endLine();
      (*out << "time = " << ApplicationTools::getTime()).endLine();
      if (itC)
      {
        (*out << "nbIterations = " << itC->getCount()).endLine();
      }
    }
  
    /////////////////////////////////////////////////////////
    //Print posterior states:
    /////////////////////////////////////////////////////////
  
    string psPath;
    if (method == "MMILS")
    {
      //Print rate states:
      psPath = ApplicationTools::getAFilePath("output.posterior.states", params, false, false);
      ApplicationTools::displayResult("Posterior states file", psPath);
      if (psPath != "none")
      {
        vector<int> states;
        ofstream out(psPath.c_str(), ios::out);
  
        HmmTools::getBestPosteriorCoalHiddenStates(*hmmLik, states, -10);
        out << ">Posterior states" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorCoalHiddenStates(*hmmLik, states, 0.5, -3, -10);
        out << ">Posterior states 0.5" << endl;
        printStates(states, out);
   
        HmmTools::getBestPosteriorCoalHiddenStates(*hmmLik, states, 0.9, -3, -10);
        out << ">Posterior states 0.9" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorCoalHiddenStates(*hmmLik, states, 0.95, -3, -10);
        out << ">Posterior states 0.95" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorCoalHiddenStates(*hmmLik, states, 0.99, -3, -10);
        out << ">Posterior states 0.99" << endl;
        printStates(states, out);

        out.close();
      }

      //Print coal rates:
      psPath = ApplicationTools::getAFilePath("output.posterior.rates", params, false, false);
      ApplicationTools::displayResult("Posterior rates file", psPath);
      if (psPath != "none")
      {
        vector<int> states;
        ofstream out(psPath.c_str(), ios::out);
  
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, -10);
        out << ">Posterior states" << endl;
        printStates(states, out);
  
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.5, -3, -10);
        out << ">Posterior states 0.5" << endl;
        printStates(states, out);
   
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.9, -3, -10);
        out << ">Posterior states 0.9" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.95, -3, -10);
        out << ">Posterior states 0.95" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.99, -3, -10);
        out << ">Posterior states 0.99" << endl;
        printStates(states, out);
      
        out.close();
      }

      //Print states and rates:
      psPath = ApplicationTools::getAFilePath("output.posterior.states_and_rates", params, false, false);
      ApplicationTools::displayResult("Posterior states & rates file", psPath);
      if (psPath != "none")
      {
        vector<int> states;
        ofstream out(psPath.c_str(), ios::out);
  
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, -10);
        out << ">Posterior states" << endl;
        printStates(states, out);
  
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.5, -3, -10);
        out << ">Posterior states 0.5" << endl;
        printStates(states, out);
   
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.9, -3, -10);
        out << ">Posterior states 0.9" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.95, -3, -10);
        out << ">Posterior states 0.95" << endl;
        printStates(states, out);
    
        HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.99, -3, -10);
        out << ">Posterior states 0.99" << endl;
        printStates(states, out);
 
        out.close();
      }
    }
    else
    {
      psPath = ApplicationTools::getAFilePath("output.posterior.states", params, false, false);
      ApplicationTools::displayResult("Posterior states file", psPath);
      if (psPath != "none")
      {
        if (hAlpha->getNumberOfStates() > 36)
        {
          ApplicationTools::displayError("Can't output states in fasta format with more than 36 different states.");
        }
        else
        {
          vector<int> states;
          ofstream out(psPath.c_str(), ios::out);
    
          HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, -10);
          out << ">Posterior states" << endl;
          printStates(states, out);
    
          HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.5, -3, -10);
          out << ">Posterior states 0.5" << endl;
          printStates(states, out);
    
          HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.9, -3, -10);
          out << ">Posterior states 0.9" << endl;
          printStates(states, out);
    
          HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.95, -3, -10);
          out << ">Posterior states 0.95" << endl;
          printStates(states, out);
    
          HmmTools::getBestPosteriorHiddenStates(*hmmLik, states, 0.99, -3, -10);
          out << ">Posterior states 0.99" << endl;
          printStates(states, out);
    
          out.close();
        }
      }
    }
  
    /////////////////////////////////////////////////////////
    //Print posterior values:
    //For the drosophila analysis... we should improve the format later!
    /////////////////////////////////////////////////////////
  
    string pvPath = ApplicationTools::getAFilePath("output.posterior.values", params, false, false);
    ApplicationTools::displayResult("Posterior values file", pvPath);
    if (pvPath != "none")
    {
      vector< vector<double> > probs;
      hmmLik->getHiddenStatesPosteriorProbabilities(probs, false);
    
      ofstream out(pvPath.c_str(), ios::out);
    
      for (size_t j = 0; j < probs[0].size(); j++)
      {
        out << "\"V" << j << "\" ";
      } 
      out << "\"Chunk\" ";
      if(annotation)
        out << "\"Annotation\"" << endl;
      else
        out << endl;

      vector<size_t>::const_iterator brkptsIt = hmmLik->getBreakPoints().begin();
      size_t current = 1;
      for (size_t i = 0; i < probs.size(); i++)
      {
        if (brkptsIt != hmmLik->getBreakPoints().end() && i == *brkptsIt)
        {
          current++;
          brkptsIt++;
        }
        out << "\"" << (i+1) << "\"";
        for (size_t j = 0; j < probs[0].size(); j++)
        {
          out << " " << probs[i][j];
        }
        out << " " << current;
        if (annotation)
          out << " " << annotation->getChar(i) << endl;
        else
          out << endl;
      }
    
      out.close();
    }

    //Print hidden states:
    string hsPath = ApplicationTools::getAFilePath("output.hidden_states", params, false, false);
    ApplicationTools::displayResult("Hidden states file", hsPath);
    if (hsPath != "none")
    {
      ofstream out(hsPath.c_str(), ios::out);
      for (size_t i = 0; i < hAlpha->getNumberOfStates(); i++)
      {
        out << TreeTools::treeToParenthesis(dynamic_cast<Tree&>(hAlpha->getState(i))) << endl;
      }
      out.close();
    }

    /////////////////////////////////////////////////////////
    //Print divergences for hidden states:
    /////////////////////////////////////////////////////////
    
    string hsdPath = ApplicationTools::getAFilePath("output.hidden_states.divergences", params, false, false);
    ApplicationTools::displayResult("Hidden states divergence file", hsdPath);
    if (hsdPath != "none")
    {
      ofstream out(hsdPath.c_str(), ios::out);
      out << "# The first column is the index of the state, then come the divergence times..." << endl;
      for (size_t i = 0; i < hAlpha->getNumberOfStates(); i++)
      {
        Tree& tree = dynamic_cast<Tree&>(hAlpha->getState(i));
        vector<int> ids = tree.getNodesId();
        out << i;
        for (size_t j = 0; j < ids.size(); j++)
        {
          if (!tree.isLeaf(ids[j]))
          {
            out << "\t" << TreeTools::getHeight(tree, ids[j]);
          }
        }
        out << endl;
      }
      out.close();
    }

    /////////////////////////////////////////////////////////
    // Compute local divergences:
    /////////////////////////////////////////////////////////
 
    psPath = ApplicationTools::getAFilePath("output.posterior.divergences", params, false, false);
    ApplicationTools::displayResult("Posterior divergences file", psPath);
    if (psPath != "none")
    {
      vector<int> states;
      ofstream out(psPath.c_str(), ios::out);
      RowMatrix<double> divergences;
      vector<int> blocks;
      vector<string> pairNames;
      HmmTools::getPosteriorDivergences(*hmmLik, divergences, blocks, pairNames);
      for (size_t i = 0; i < pairNames.size(); ++i)
        out << pairNames[i] << "\t";
      out << "Block" << endl;
      for (size_t i = 0; i < divergences.getNumberOfRows(); ++i) {
        for (size_t j = 0; j < divergences.getNumberOfColumns(); ++j) {
          out << divergences(i, j) << "\t";
        }
        out << blocks[i] << endl;
      }
      out.close();
    }


    /////////////////////////////////////////////////////////
    // Compute variance covariance matrix of parameters:
    /////////////////////////////////////////////////////////
 
    string vcvmPath = ApplicationTools::getAFilePath("output.estimated.vcv", params, false, false);
    ApplicationTools::displayResult("Parameters Var-Covar matrix output file", vcvmPath);
    //Restore parameters in case a transformation was used during optimization.
    //Otherwise, just update values.
    ParameterList pl = hmmLik->getParameters();
    if (vcvmPath != "none")
    {
      ParameterList pl2;
      vector<string> vcvPNames = ApplicationTools::getVectorParameter<string>("output.estimated.vcv.select", params, ',', "all", "", false, false);
      if(vcvPNames.size() == 1 && vcvPNames[0] == "all")
      {
        pl2 = pl;
      }
      else
      {
        for (size_t i = 0; i < vcvPNames.size(); i++)
        {
          vcvPNames[i] = TextTools::removeWhiteSpaces(vcvPNames[i]);
          if (!pl.hasParameter(vcvPNames[i]))
            ApplicationTools::displayWarning("Parameter " + vcvPNames[i] + " not known and was hence ignored.");
          else pl2.addParameter(pl.getParameter(vcvPNames[i]));
        }
      }
      ThreePointsNumericalDerivative* fun = new ThreePointsNumericalDerivative(hmmLik);
      fun->enableSecondOrderCrossDerivatives(true);
      vector<string> names = pl2.getParameterNames();
      fun->setParametersToDerivate(names);
      RowMatrix<double>* hessian = NumTools::computeHessianMatrix(*fun, pl2);
      string hessianPath = ApplicationTools::getAFilePath("output.estimated.hessian", params, false, false);
      ApplicationTools::displayResult("Parameters Hessian matrix output file", hessianPath);
     
      ofstream outh(hessianPath.c_str(), ios::out);
      outh << names[0];
      for(size_t j = 1; j < names.size(); j++)
        outh << "," << names[j];
      outh << endl;
      for (size_t i = 0; i < names.size(); i++)
      {
        outh << names[i];
        for (size_t j = 0; j < names.size(); j++)
        {
          outh << "," << (*hessian)(i,j);
        }
        outh << endl;
      }
      outh.close();

      RowMatrix<double> vcvMat;
      MatrixTools::inv(*hessian, vcvMat);
      ofstream out(vcvmPath.c_str(), ios::out);
      out << names[0];
      for (size_t j = 1; j < names.size(); j++)
        out << "," << names[j];
      out << endl;
      for (size_t i = 0; i < names.size(); i++)
      {
        out << names[i];
        for (size_t j = 0; j < names.size(); j++)
          out << "," << vcvMat(i,j);
        out << endl;
      }
      out.close();
    }



    /////////////////////////////////////////////////////////
    // Likelihood profiling analysis
    /////////////////////////////////////////////////////////
    
    string likProfPath = ApplicationTools::getAFilePath("output.likelihood_profiling.file", params, false, false);
    if (likProfPath != "none")
    {
      ParameterGrid* grid = NumCalcApplicationTools::getParameterGrid(params, "", true, false);
      VVdouble* values = FunctionTools::computeGrid(*hmmLik, *grid);
      ofstream out(likProfPath.c_str(), ios::out);
      for (unsigned int i = 0; i < grid->getNumberOfDimensions(); i++)
        out << grid->getDimensionName(i) << "\t";
      out << "lnL" << endl;
      for (size_t i = 0; i < values->size(); i++)
      {
        for (size_t j = 0; j < (*values)[i].size(); j++)
        {
          if (j > 0) out << "\t";
          if (j == (*values)[i].size() - 1) out << setprecision(20);
          else out << setprecision(8);
          out << (*values)[i][j];
        }
        out << endl;
      }
      out.close();
      delete grid;
      delete values;
    }

    //Cleaning...
    delete hmmLik;
    if (itC) delete itC;

  }
  else if (analysis == "simulate")
  {
    CoalHmmSimulator* hmmsim = new CoalHmmSimulator(hModel, sModel, rDist, false);

    size_t nbSites = ApplicationTools::getParameter<size_t>("number_of_sites", params, 100);
    ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
    ApplicationTools::displayTask("Perform simulations");
    SequenceContainer * sites = hmmsim->simulate(nbSites);
    ApplicationTools::displayTaskDone();
    
    // Write to file:
    SequenceApplicationTools::writeSequenceFile(*sites, params);

    // Write topology sequence:
    string tsPath = ApplicationTools::getAFilePath("output.true_states", params, false, false);
    ApplicationTools::displayResult("True states file", tsPath);
    if (tsPath != "none")
    {
      ofstream tsFile(tsPath.c_str(), ios::out);
      tsFile << ">genealogies" << endl;
      vector<size_t> ts = hmmsim->getGenealogySequence();
      for(size_t i = 0; i < ts.size(); i++) tsFile << ts[i];
      tsFile << endl;
      tsFile.close();
    }
 
    delete hmmsim;
  }
  else throw Exception("Unknown running mode.");

  delete alpha;
  delete sModel;
  delete rDist;
  coalhmm.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return(1);
  }

  return (0);
}

