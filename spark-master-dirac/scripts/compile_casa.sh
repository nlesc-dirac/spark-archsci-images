#!/bin/bash

# CASACORE
mkdir -p /opt/soft/casacore/data
cd /opt/soft/casacore/data
wget -c ftp://ftp.astron.nl/outgoing/Measures/WSRT_Measures.ztar
tar zxfv WSRT_Measures.ztar && rm -f WSRT_Measures.ztar


cd /tmp
rm -rf casacore_install
git clone --progress --verbose https://github.com/casacore/casacore.git casacore_install
cd casacore_install


mkdir build
cd build
cmake -DUSE_FFTW3=ON -DCMAKE_INSTALL_PREFIX=/opt/soft/casacore -DDATA_DIR=/opt/soft/casacore/data -DUSE_OPENMP=ON \
    -DUSE_HDF5=ON -DBUILD_PYTHON=ON -DUSE_THREADS=ON ..
make -j4
make install
