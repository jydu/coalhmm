from sim_pipeline import *
import sys
import os
import os.path

outputdir_root = "Results_ILS09/Sim"
NeRef = N1 = N2 = N3 = 30000.0
N12 = 90000.0
N123 = 65000.0
g = 20
tau1 = 4.2e6 / g / (2 * NeRef)
tau12 = 6.6e6 / g / (2 * NeRef)
L = 1e6 # 1Mb
r = 1.6e-8
u = 1e-9
tau123 = 18e6 / g / (2 * NeRef) #Fixed divergence with the outgroup
tau3 = tau123 - tau12
tau2 = tau12 - tau1


for i in range(int(sys.argv[1]), int(sys.argv[2])+1) :
    outputdir = outputdir_root + str(i)
    if not os.path.exists(outputdir) :
        os.makedirs(outputdir)

    populationSpec = Po(N123 / NeRef,
                Me(tau12, [Po(N3 / NeRef, Sa(1)),
                          Po(N12 / NeRef, Me(tau1,
                                           [Po(N1 / NeRef, Sa(1)), Po(N2 / NeRef, Sa(1))]))]))

    runcoasim(outputdir, populationSpec, NeRef, L, r, outgroup=tau123)	
    
    runbppseqgen(outputdir, seqfile=outputdir + "/seq.ras" , logfile=outputdir + "/bppseqgen.ras.out", L=L, NeRef=NeRef, u=u, g=g)
    
    runcoalhmmILS09(outputdir, seqfile=outputdir + "/seq.ras", logfile=outputdir + "/coalhmm.out",
        tau1=tau1, tau2=tau2, tau3=tau3, N12=N12, N123=N123, NeRef=NeRef, u=u, g=g ,r=r, output_post=True)
     
    runpostdecode(seqfile=outputdir + "/seq.ras", nb_states=4)

    print  "Done with", i
 
