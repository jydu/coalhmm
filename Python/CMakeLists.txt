# CMake script for Bio++ CoalHMM
# Author: Julien Dutheil
# Created: 03/09/2009

set(coalhmm_SCRIPTS
 post_decode.py
 post_get_segments.py
 post_mmhmm2hmm.py
 post_divergence.py
)

# Install python scripts:
install(PROGRAMS ${coalhmm_SCRIPTS} DESTINATION bin)

