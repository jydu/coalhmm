#!/usr/bin/python

#
# File: post_mmhmm2mm.py
# Created by: Julien Dutheil
# Created on: Aug 05 2009
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
import math
import sys
import getopt

# sep is a separator character used in the output file.
def process (inf, outf, sep, nRates, nGens):
    for i in range(nGens) :
        outf.write("%sV%i" % (sep, i))
    outf.write("\n")

    # Skip the first line.
    inf.next()

    # Now loops over all lines.
    count = 0
    for line in inf:
        count += 1
        probs = line.split()[1:1 + nRates * nGens]
        outf.write("%d" % count)
        for i in range(nGens) :
            prob = 0
            for j in range(nRates) :
                prob += float(probs[i * nRates + j])
            outf.write("%s%f" % (sep, prob))
        #TODO: here we should copy the annotation from the original file, if any.
        outf.write("\n")

# sep is a separator character used in the output file.
def process_fnames (in_name, out_name, sep, nRates, nGens):
    inf  = open (in_name, 'r')
    outf = open (out_name, 'w')
    process (inf, outf, sep, nRates, nGens)
    inf.close()
    outf.close()

def help():
    print "This script takes as input the posterior probabilities file output"
    print "from the CoalHMM program and generated from a Markov-Modulated model"
    print "and convert it to a 'standard', non-modulated posterior probability file,"
    print "with sum probabilities for gneealogies states only, summing over all"
    print "possible rate classes."
    print
    print "The first row is ignored and is supposed to contain a header line."
    print

def usage():
    print "Usage: ./post_mmhmm2hmm.py <flags> infile outfile"
    print "Flags:"
    print "\t-r<integer> / --nb_rates=<integer>:"
    print "\tSets the number of rate classes (defaults to 4)."
    print "\t-g<integer> / --nb_genealogies=<integer>:"
    print "\tSets the number of genealogies classes (defaults to 4)."
    print "\t-h / --help: Prints extended help message."
    print "\t-s<sep> / --separator=<sep>:"
    print "\tSets the separator string to be used in the output file (tab is the default.)"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:r:g:", ["help", "separator=", "nb_rates=", "nb_genealogies="])
    sep = '\t'
    nRates = 4
    nGens = 4
    for o, a in opts:
        if o in ("-s", "--separator"):
            sep = a
        elif o in ("-r", "--nb_rates"):
            nRates = int(a)
        elif o in ("-g", "--nb_genealogies"):
            nGens = int(a)
        elif o in ("-h", "--help"):
            help()
            sys.exit()
    in_name, out_name = args
except:
    usage()
    sys.exit(2)

process_fnames (in_name, out_name, sep, nRates, nGens)
