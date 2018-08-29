#!/bin/bash

pacman --needed --noconfirm -S blas
#pacman --noconfirm -R blas
#sudo -H -u archsci bash -c 'yaourt -S --needed --noconfirm openblas'

# clone Sagecal repo
mkdir -p /opt/dirac
git clone https://github.com/nlesc-dirac/sagecal-on-spark.git /opt/dirac

## SAGECAL
cd /opt/dirac/excon/JAVA && make -f Makefile.nocuda
cd /opt/dirac/scripts && sh download_sagecal_dataset.sh

rm -rf /var/cache/pacman/pkg/* /home/archsci/temp
