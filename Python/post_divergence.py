#!/usr/bin/python

#
# File: post_divergence.py
# Created by: Julien Dutheil
# Created on: Oct 23 2009
#

#This file is part of the CoalHMM program and library.
#
#CoalHMM is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#any later version.
#
#CoalHMM is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with CoalHMM.  If not, see <http://www.gnu.org/licenses/>.

import os
from math import *
import sys
import getopt

# sep is a separator character used in the output file.
def process (inf, outf, sep, cols, theta, tau, median) :
    # Compute the [truncated] exponential distribution
    bounds = range(cols +1)
    values = range(cols)
    if tau > 0 :
        #Truncated exp
        total_prob = 1 - exp(-tau/theta);
        bounds[0] = 0
        for i in range(1, cols+1) :
            a = bounds[i-1]
            b = theta * log(float(cols) / (float(cols) - i * total_prob));
            bounds[i] = b;
            if median :
                values[i-1] = theta * log(2 * cols / (2 * cols - (2 * i - 1) * total_prob)); 
            else :
                values[i-1] = theta + (a*exp(-a/theta) - b*exp(-b/theta)) / (exp(-a/theta) - exp(-b/theta));
        bounds[cols] = tau;     
    else :
        #Full exp
        bounds[0] = 0
        for i in range(1, cols+1) :
            a = bounds[i-1]
            if i == cols :
                b = exp(100) #infinity (sort of...)
            else :
                b = theta * log(float(cols) / float(cols - i));
            bounds[i] = b
            if median :
                values[i-1] = theta * log(2*cols / (2*cols - i) + 1); 
            else :
                values[i-1] = cols * ((a + theta) * exp(-a / theta) - (b + theta) * exp(-b / theta));

    outf.write("Divergence%sShannon\n" % (sep))
    # Skip the first line.
    inf.next()
    for line in inf:
        shannon_index = 0
        post = 0.0
        i = 0
        for s in line.split()[1:1 + cols]:
            prob = float(s)
            #if prob > max_prob:
            #    max_prob = prob
            #    max_i    = i
            if prob > 0 :
                shannon_index -= prob * log(prob)
            post += prob * values[i]
            i += 1
        outf.write("%f%s%f\n" %
                   (post, sep, shannon_index))

# sep is a separator character used in the output file.
def process_fnames (in_name, out_name, sep, cols, theta, tau, median):
    inf  = open (in_name, 'r')
    outf = open (out_name, 'w')
    process (inf, outf, sep, cols, theta, tau, median)
    inf.close()
    outf.close()

def help():
    print "This script takes as input the posterior probabilities file output"
    print "from the CoalHMM program (support only the Divergence model for now),"
    print "and the estimated theta value. It computes for each site the"
    print "posterior estimate of divergence and the Shannon information index."
    print "The input format is:"
    print
    print "\t<ignored> <state0> <state1> ... <state<columns>> <ignored <ignored>"
    print
    print "The first row is ignored and is supposed to contain a header line."
    print
    print "The output format is:"
    print
    print "\t<posterior divergence> <Shannon index>"
    print
    print "The fields are separated by a separator character."
    print

def usage():
    print "Usage: ./post_divergence.py <flags> infile outfile"
    print "Flags:"
    print "\t-c<columns> / --columns=<columns>:"
    print "\tSets the number of columns (defaults to 4.)"
    print "\t-h / --help: Prints extended help message."
    print "\t-s<sep> / --separator=<sep>:"
    print "\tSets the separator string to be used in the output file (tab is the default.)"
    print "\t-t<unsigned double> / --theta=<unsigned double>:"
    print "\tSets the estimated value of theta."
    print "\t-T<double> / --tau=<double>:"
    print "\tSets the estimated value of tau, if any. A negative value [the default]"
    print "\tset tau to infinity, which corresponds to a model without outgroup."
    print "\t-m / --median:"
    print "\tTell if the median value was used instead of the mean in the model."


try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:c:t:T:m", ["help", "separator=",
                                                       "columns=", "theta=", "tau=", "median"])
    sep    = '\t'
    cols   = 4
    theta  = -1
    tau    = -1
    median = False
    for o, a in opts:
        if o in ("-s", "--separator"):
            sep = a
        elif o in ("-c", "--columns"):
            cols = int(a)
        elif o in ("-t", "--theta"):
            theta = float(a)
        elif o in ("-T", "--tau"):
            tau = float(a)
        elif o in ("-m", "--median"):
            median = True
        elif o in ("-h", "--help"):
            help()
            sys.exit()
    if theta <= 0 :
        print "Argument 'theta' is missing or unvalid."
        sys.exit(2)

    in_name, out_name = args
except:
    usage()
    sys.exit(2)

process_fnames (in_name, out_name, sep, cols, theta, tau, median)
