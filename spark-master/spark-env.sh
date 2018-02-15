#!/usr/bin/env bash

export JAVA_HOME=/usr/lib/jvm/default

export SPARK_DIST_CLASSPATH=$(/opt/soft/hadoop/bin/hadoop classpath)

export SPARK_MASTER_IP=localhost
export SPARK_LOCAL_IP=localhost

export SPARK_HOME=/opt/soft/spark
export SPARK_MASTER_PORT=7077
export SPARK_MASTER_WEBUI_PORT=8080
export SPARK_MASTER_LOGFILE=/opt/soft/spark/logs/spark_master.log

export SPARK_MASTER="spark://sparkmaster:7077"
export SPARK_WORKER_WEBUI_PORT=8081
export SPARK_WORKER_LOGFILE=/opt/soft/spark/logs/spark_worker.log

#export HADOOP_HOME=/opt/soft/hadoop
#export HADOOP_OPTS="-Djava.library.path=$HADOOP_HOME/lib"

export JAVA_HOME=/usr/lib/jvm/default
export SPARK_VERSION=2.2.1
export HADOOP_VERSION=2.7.5
export SPARK_HOME=/opt/soft/spark
export HADOOP_HOME=/opt/soft/hadoop
export HADOOP_NATIVE=/opt/soft/hadoop/lib/native/

export JAVA_LIBRARY_PATH=$JAVA_LIBRARY_PATH:$HADOOP_NATIVE
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HADOOP_NATIVE:/opt/soft/casacore/lib
export SPARK_YARN_USER_ENV="JAVA_LIBRARY_PATH=$JAVA_LIBRARY_PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
export HADOOP_COMMON_LIB_NATIVE_DIR=$HADOOP_HOME/lib/native
export HADOOP_OPTS="$HADOOP_OPTS -Djava.library.path=$HADOOP_HOME/lib/native"
