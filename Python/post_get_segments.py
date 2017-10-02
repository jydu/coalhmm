#!/usr/bin/python

#
# File: post_get_segments.py
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
import sys
from optparse import OptionParser

def stretches(inf, outf, sep, threshold):
    print >> outf, "State", "Start", "End", "AverageShannon"
    prev = -1
    cnt = 0
    start = 0
    end = 0
    # Skip the first line.
    inf.next()
    for line in inf:
        frags = line.split(sep)
        if frags[0] == "NA" :
            # Pb, probabilities were not summing to one. We assume that state is the same as previous position!
            pass
        else :
            state = int(frags[0])
            index = float(frags[2])
        if index < threshold:
            if 0 <= prev :
                print >> outf, prev, start, end, (float(isum) / cnt)
            prev = -1
            end += 1
        elif 0 > prev:
            prev = state
            cnt = 1
            start = end = end + 1
            isum = index
        elif state != prev:
            print >> outf, prev, start, end, (float(isum) / cnt)
            prev = state
            start = end = end + 1
            cnt = 1
            isum = index
        else:
            cnt += 1
            end += 1
            isum += index
    if 0 <= prev:
        print >> outf, prev, start, end, (float(isum) / cnt)

parser = OptionParser("usage: %prog [options] infile outfile")
parser.add_option("-s", "--separator", dest="sep", default='\t',
                  help="set the separator character. (tab is the default)")
parser.add_option("-t", "--threshold",
                  type="float", dest="threshold", default=0.0,
                  help="ignore all base pairs with Shannon index below this value.")

try:
    options, (in_name, out_name) = parser.parse_args()
except:
    parser.print_help()
    sys.exit(2)

inf = open (in_name, 'r')
outf = open (out_name, 'w')
stretches(inf, outf, options.sep, options.threshold)
outf.close()
inf.close()
