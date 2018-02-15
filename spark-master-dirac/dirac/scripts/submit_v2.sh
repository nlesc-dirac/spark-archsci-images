#!/bin/sh


/spark/bin/spark-submit \
    --driver-class-path /opt/dirac/excon/JAVA \
    --class Driver \
    --master spark://spark-master:6066 \
    --deploy-mode cluster  \
    --executor-cores 1 \
    --driver-library-path /usr/lib \
    Driver.jar -m /host/datasets/sm.ms

#java -Djava.library.path=$(pwd) -cp . Driver -m /host/sm.ms

