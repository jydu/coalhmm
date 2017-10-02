
#!/usr/local/bin/python

# This python scripts allows you to perform simulations as describe for instance in Dutheil et al, Genetics, 2009.
# It depends on:
# CoaSim
# The Bio++ Program Suite
# The coalhmm executable

# Authors: Troels Toftebjerg Hansen, Kasper Munch, Julien Dutheil

import CoaSim
from CoaSim.popStructure import Population as Po, Sample as Sa, Merge as Me, Migration as Mi

import re
import os
from random import random

coalhmm_exe = "~/.local/bin/coalhmm --noninteractive=yes"
postdecode_exe = "~/.local/bin/post_decode.py"
coalhmmsim_exe = "~/.local/bin/coalhmm --noninteractive=yes param=options.coalhmmsim"
bppseqgen_exe = "~/.local/bin/bppseqgen --noninteractive=yes param=options.bppseqgen"

def runcoasim(outputdir, populationSpec, NeRef, L, r, outgroup, migrationSpec = False):
    out = open(outputdir + "/trees.dnd", "w")
    
    if migrationSpec :
        arg = CoaSim.simulate([], populationSpec, migrationSpec=migrationSpec, rho = 4 * NeRef * L * r, keepEmptyIntervals=True)
    else :
        arg = CoaSim.simulate([], populationSpec, rho = 4 * NeRef * L * r, keepEmptyIntervals=True)

    for i in arg.intervals:
      s = str(i.tree)
      if outgroup > 0:
        s.replace(')', ' ')
        s = s.replace(')', ' ')
        s = s.replace('(', ' ')
        s = s.replace("'0'", ' ')
        s = s.replace("'1'", ' ')
        s = s.replace("'2'", ' ')
        s = s.replace(':', ' ')
        s = s.replace(';', ' ')
        s = s.replace(',', ' ')
        s = s.replace('  ', ' ')
        s = s.replace('  ', ' ')
        s = s.replace('  ', ' ')
        s = s.strip()
        ss = s.split()
        max = 0
        for e in ss:
        	if float(e) > max:
        		max = float(e)
      		
        s = str(i.tree)
        s = s.replace(';', '')
        val = outgroup - max
        s = '(' + s + ' : ' + str(val) + ", 'outgroup' : " + str(outgroup) + ');'
      
      out.writelines( str(i.start) + " " + str(i.end) + " " + s  + "\n")

    out.close()
    
def runbppseqgen(outputdir, seqfile, logfile, L, NeRef, u, g, extra_options=""):
    i = 0
    fin = open(outputdir + "/trees.dnd", "r")
    scaledtree = open(outputdir + "/scaledtrees.dnd", "w")
    for tree in fin :
        scaledtree.writelines(transformtree(tree, NeRef, u, g))
        i = i + 1
        if i % 1000 == 0 :
            print "Transformed " + str(i) + " trees"

    scaledtree.close()
    print "running bppseqgen"
    os.system(bppseqgen_exe + " number_of_sites=" + \
        str(L) + " input.tree.file=" + outputdir + "/scaledtrees.dnd output.sequence.file=" + seqfile + ".fasta " + extra_options + " > " + logfile)


treelengthmatcher = re.compile(': -?\d*\.?\d*e?-?\d*[,)]')
	
def transformtree(tree, NeRef, u, g):
    lens = treelengthmatcher.findall(tree);
    for l in lens :
        #print l
        last = l[-1]
        lc = l.replace(': ', '')
        lc = lc.replace(last, '')
        fl = float(lc) * 2 * NeRef * u * g
        #print fl
        tree = tree.replace(l, ': ' + str(fl) + last)

    #strip beginning
    #tree = tree.lstrip('0123456789.e- ')
    return tree
		
#def runcoalhmmsim(outputdir, seqfile, logfile, method="ILS09", outgroup="", extra_options="", suffix="") :
#    if method == 'ILS09' :
#        param = "tau1=" + str(tau1 * u * 2 * NeRef * g) + " tau2=" + str((tau12-tau1) * u * 2 * NeRef * g) + " c2=" +\
#            str(((tau123- tau12) * NeRef - 4/3 * N123) * u * 2 * g) + " theta1=" + str(N12 * u * 2 * g) + " theta2=" + \
#            str(N123 * u * 2 * g) + " rho=" + str(r / (u * g)) + " DATA=" + seqfile + " SUFFIX=" + suffix
#    elif method == 'Div' :
#        if outgroup == "":
#            param = "tau1=" + str(tau1 * u * 2 * NeRef * g) + " tau2=" + " theta12=" + str(N12 * u * 2 * g) + " theta1=theta12 theta2=theta12" + \
#                " rho=" + str(r / (u * g)) + " outgroup=" + " DATA=" + seqfile + " SUFFIX=" + suffix
#        else :
#            param = "tau1=" + str(tau1 * u * 2 * NeRef * g) + " tau2=" + str((tau12-tau1) * u * 2 * NeRef * g) + \
#                " theta12=" + str(N12 * u * 2 * g) + " theta1=theta12 theta2=theta12" + \
#                str(N123 * u * 2 * g) + " rho=" + str(r / (u * g)) + " outgroup=" + outgroup + " DATA=" + seqfile + " SUFFIX=" + suffix     
#    else :
#        print "Unknown coalhmm method : %s" % method
#
#    print "running coalhmmsim"
#    os.system(coalhmmsim_exe + " " + param + " " + extra_options + " > " + logfile)

def jitter(x, sd = 1.) :
    r = (random() * 2 - 1) * sd
    return x * (1. + r)

def runcoalhmmILS09(outputdir, seqfile, logfile, tau1, tau2, tau3, N12, N123, NeRef, u, g, r, output_post, extra_options="", suffix="", noise=0.2) :
    param = "tau1=" + str(tau1 * u * 2 * NeRef * g) + " tau2=" + str(tau2 * u * 2 * NeRef * g) + " c2=" +\
        str((tau3 * NeRef - 4/3 * N123) * u * 2 * g) + " theta1=" + str(N12 * u * 2 * g) + " theta2=" + \
        str(N123 * u * 2 * g) + " rho=" + str(r / (u * g)) + " DATA=" + seqfile + " SUFFIX=" + suffix
    if not output_post :
      param += " output.posterior.values=None" 

    print "running coalhmm"
    os.system(coalhmm_exe + " param=options.coalhmm_ILS09.bpp " + param + " " + extra_options + " > " + logfile)

def runpostdecode(seqfile, nb_states) :
    print(postdecode_exe + " -c" + str(nb_states) + " " + seqfile + ".hmm.values.csv " + seqfile + ".hmm.postdecode.csv")
    os.system(postdecode_exe + " -c" + str(nb_states) + " " + seqfile + ".hmm.values.csv " + seqfile + ".hmm.postdecode.csv")
 
def runcoalhmmDiv(outputdir, seqfile, logfile, tau1, tau2, N12, NeRef, u, g, r, output_post, extra_options="", suffix="", outgroup="", nbClasses=10, noise=0.2) :
    if outgroup == "":
        init_tau1 = jitter(tau1 * u * 2 * NeRef * g, noise)
        init_theta12 = jitter(N12 * u * 2 * g, noise)
        init_rho = jitter(r / (u * g), noise)
        param = "tau1=" + str(init_tau1) + "tau2= theta12=" + str(init_theta12) + " theta1=theta12 theta2=theta12" + \
            " rho=" + str(init_rho) + " outgroup=" + " nbClasses=" + str(nbClasses) + " DATA=" + seqfile + " SUFFIX=" + suffix + "." + str(nbClasses) + "classes"
    else :
        init_tau1 = jitter(tau1 * u * 2 * NeRef * g, noise)
        init_tau2 = jitter(tau2 * u * 2 * NeRef * g, noise)
        init_theta12 = jitter(N12 * u * 2 * g, noise)
        init_rho = jitter(r / (u * g), noise)
        param = "tau1=" + str(init_tau1) + " tau2=" + str(init_tau2) + \
            " theta12=" + str(init_theta12) + " theta1=theta12 theta2=theta12" + \
            " rho=" + str(init_rho) + " \"outgroup=" + outgroup + "\" nbClasses=" + str(nbClasses) + " DATA=" + seqfile + " SUFFIX=" + suffix + "." + str(nbClasses) + "classes"     
    if not output_post :
      param += " output.posterior.values=None" 

    print "running coalhmm"
    print(coalhmm_exe  + " param=options.coalhmm_Div.bpp " + param + " " + extra_options + " > " + logfile)
    os.system(coalhmm_exe  + " param=options.coalhmm_Div.bpp " + param + " " + extra_options + " > " + logfile)



   
