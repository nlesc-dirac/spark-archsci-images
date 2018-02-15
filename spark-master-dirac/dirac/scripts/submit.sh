#!/bin/sh


/spark/bin/spark-submit --driver-class-path /host/excon_src/JAVA \
    --class Driver \
    --master spark://spark-master:6066 \
    --deploy-mode cluster  \
    --executor-cores 1 \
    --driver-library-path /usr/local/lib \
    Driver.jar
