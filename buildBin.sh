#! /bin/sh
arch=`uname -m`
version=1.0.4-1

strip CoalHMM/coalhmm
tar cvzf CoalHMM-${arch}-bin-static-${version}.tar.gz CoalHMM/coalhmm

