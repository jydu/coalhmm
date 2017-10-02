#!/usr/bin/python

#
# File: post_decode.py
# Created by: Joergen Fogh, Julien Dutheil
# Created on: Jul 08 2009
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
def process (inf, outf, sep, cols):
    outf.write("State%sProbmax%sShannon\n" % (sep, sep))
    # Skip the first line.
    inf.next()
    for line in inf:
        shannon_index = 0
        max_prob = 0.0
        sum_prob = 0.0
        i = 0
        for s in line.split()[1:1 + cols]:
            prob = float(s)
            sum_prob += prob
            if prob > max_prob:
                max_prob = prob
                max_i    = i
            if prob > 0 :
                shannon_index -= prob * math.log(prob)
            i += 1
        if sum_prob > 1.0001 :
            outf.write("NA%sNA%sNA\n" %
                       (sep, sep))
        else :
            outf.write("%d%s%f%s%f\n" %
                       (max_i, sep, max_prob, sep, shannon_index))

# sep is a separator character used in the output file.
def process_fnames (in_name, out_name, sep, cols):
    inf  = open (in_name, 'r')
    outf = open (out_name, 'w')
    process (inf, outf, sep, cols)
    inf.close()
    outf.close()

def help():
    print "This script takes as input the posterior probabilities file output"
    print "from the CoalHMM program and outputs the index of the state with"
    print "the maximum probability, the corresponding probability and the"
    print "Shannon information index for each base pair. The input format is:"
    print
    print "\t<ignored> <state0> <state1> ... <state<columns>> <ignored <ignored>"
    print
    print "The first row is ignored and is supposed to contian a header line."
    print
    print "The output format is:"
    print
    print "\t<most probable state> <probability of state> <Shannon index>"
    print
    print "The fields are separated by a separator character."
    print

def usage():
    print "Usage: ./post_decode.py <flags> infile outfile"
    print "Flags:"
    print "\t-c<columns> / --columns=<columns>:"
    print "\tSets the number of columns (defaults to 4.)"
    print "\t-h / --help: Prints extended help message."
    print "\t-s<sep> / --separator=<sep>:"
    print "\tSets the separator string to be used in the output file (tab is the default.)"

try:
    opts, args = getopt.getopt(sys.argv[1:], "hs:c:", ["help", "separator=",
                                                       "columns="])
    sep = '\t'
    cols = 4
    for o, a in opts:
        if o in ("-s", "--separator"):
            sep = a
        elif o in ("-c", "--columns"):
            cols = int(a)
        elif o in ("-h", "--help"):
            help()
            sys.exit()
    in_name, out_name = args
except:
    usage()
    sys.exit(2)

process_fnames (in_name, out_name, sep, cols)
