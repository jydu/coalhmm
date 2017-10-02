//
// File: IlsCoalHmmTransitionMatrix07.h
// Created by: Julien Dutheil
// Created on: Mon Oct 29 09:08 2007
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

#ifndef _BASICCOALHMMTRANSITIONMATRIX07_H_
#define _BASICCOALHMMTRANSITIONMATRIX07_H_

#include "CoalHmmTransitionMatrix.h"

class IlsCoalHmmTransitionMatrix07:
public AbstractCoalHmmTransitionMatrix
{
protected:
    double sigma_, mu_, nu_, lambda_;
    bool onev_;
    
public:
    IlsCoalHmmTransitionMatrix07(const CoalHmmStateAlphabet* hiddenAlphabet, double s, double u, double v1, double v2, bool onev = true);
    
    IlsCoalHmmTransitionMatrix07* clone() const { return new IlsCoalHmmTransitionMatrix07(*this); }
    
public:
    void fireParameterChanged(const bpp::ParameterList& pl)
    {
        sigma_ = getParameter_(0).getValue();
        nu_    = getParameter_(1).getValue();
        mu_    = getParameter_(2).getValue();
        if(!onev_)
            lambda_ = getParameter_(3).getValue();
        actualizeTransitionMatrix();
    }
    
    void printUserFriendlyParameters(bpp::OutputStream& out) const;
    
protected:
    void actualizeTransitionMatrix();
    
};

#endif //_ILSCOALHMMTRANSITIONMATRIX07_H_

