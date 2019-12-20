#! /bin/sh
arch=x86_64 #i686
version=1.0.0-2

cd build
strip CoalHMM/coalhmm
tar cvzf CoalHMM-${arch}-bin-static-${version}.tar.gz CoalHMM/coalhmm

