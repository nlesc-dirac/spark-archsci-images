#!/usr/sbin/zsh

source /opt/soft/spark/conf/spark-env.sh

#export SPARK_WORKER_HOST=$(hostname)
export SPARK_WORKER_HOST=$(ifconfig eth1 | grep "inet" | cut -d " " -f10)

touch $SPARK_WORKER_LOGFILE

sleep 5s
$SPARK_HOME/sbin/start-slave.sh \
    --host $SPARK_WORKER_HOST \
    --webui-port $SPARK_WORKER_WEBUI_PORT \
    $SPARK_MASTER \
    >> $SPARK_WORKER_LOGFILE &

    #--host $SPARK_WORKER_HOST \

sleep 5s

tail -f $SPARK_WORKER_LOGFILE
