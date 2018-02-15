#!/bin/bash

source /opt/soft/spark/conf/spark-env.sh

export SPARK_MASTER_HOST=sparkmaster
#export SPARK_MASTER_HOST=localhost
#export SPARK_MASTER_HOST=$(hostname)
#export SPARK_MASTER_HOST=$(ifconfig eth1 | grep "inet" | cut -d " " -f10)

touch $SPARK_MASTER_LOGFILE

$SPARK_HOME/sbin/start-master.sh \
    --host $SPARK_MASTER_HOST \
    --port $SPARK_MASTER_PORT \
    --webui-port $SPARK_MASTER_WEBUI_PORT \
    >> $SPARK_MASTER_LOGFILE &

sleep 5s

tail -f $SPARK_MASTER_LOGFILE
