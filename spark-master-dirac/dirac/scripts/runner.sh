#!/bin/sh

ids=$(hadoop fs -stat "%n" /opt/dirac/datasets/sm$1*)

for fname in $ids
do 
    echo $fname

    /opt/soft/spark/bin/spark-submit \
        --master spark://spark-master:6066 \
        --deploy-mode cluster  \
        --executor-cores 1 \
        --driver-class-path /opt/dirac/excon/JAVA \
        --driver-library-path /opt/dirac/excon/JAVA \
        --class Driver \
        /opt/dirac/excon/JAVA/Driver.jar -m  /opt/dirac/datasets/$fname

done






