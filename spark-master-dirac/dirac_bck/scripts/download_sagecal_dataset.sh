#!/bin/bash

mkdir -p /opt/dirac/datasets && cd /opt/dirac/datasets
wget -c \
    https://raw.githubusercontent.com/nlesc-dirac/sagecal/master/test/sm.ms.tar && \
    tar -xf sm.ms.tar
