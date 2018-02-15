#!/bin/sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib:/opt/soft/casacore:/host/excon_src_v2/JAVA:/opt/hadoop-2.8.0/lib/native/

#export CLASSPATH=$(hadoop classpath --glob)
export CLASSPATH=$(hadoop classpath)
#echo $CLASSPATH

java -Djava.library.path=$(pwd) -jar Driver.jar -m /host/datasets/sm.ms


