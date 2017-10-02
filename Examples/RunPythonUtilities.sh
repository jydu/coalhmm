#! /bin/sh

#First we should run the ILS example:
coalhmm param=ILS.bpp

#Get the posterior decoding out of the posterior probabilites:
python ../Python/post_decode.py -c4 sim-ILS09_GTR_Uniform.posterior.values.csv sim-ILS09_GTR_Uniform.posterior.states.csv

#Now we take the posterior decoding and retrieve all "segments", that is continuous parts on the genome in a given genealogy.
#It is also possible to set a minimum information threshold to remove sites with too few information.
python ../Python/post_get_segments.py sim-ILS09_GTR_Uniform.posterior.states.csv sim-ILS09_GTR_Uniform.posterior.segments.csv

