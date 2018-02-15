#!/bin/bash

#cp /etc/pacman.d/mirrorlist /etc/pacman.d/mirrorlist.backup
#rankmirrors -n 6 /etc/pacman.d/mirrorlist.backup > /etc/pacman.d/mirrorlist

# DEPENDENCIES
pacman -Syyuu --needed --noconfirm --force

pacman -S --needed --noconfirm \
    cmake extra-cmake-modules \
    gcc-fortran \
    mpfr hdf5 cfitsio wcslib \
    python2-numpy boost fftw \
    openmp openmpi python3-numpy

sudo -H -u archsci bash -c 'yaourt -S --needed --noconfirm sofa'

rm -rf /var/cache/pacman/pkg/* /home/archsci/temp
