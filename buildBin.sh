#! /bin/sh
arch=x86_64 #i686
version=1.0.0-1

strip CoalHMM/coalhmm
tar cvzf CoalHMM-${arch}-bin-static-${version}.tar.gz CoalHMM/coalhmm

