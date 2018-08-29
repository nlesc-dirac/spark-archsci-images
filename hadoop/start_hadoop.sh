#!/bin/bash

# source /etc/conf.d/hadoop
# . /etc/conf.d/hadoop


# echo "Start DFS"
# /opt/soft/hadoop/sbin/start-dfs.sh
#
# /opt/soft/hadoop/bin/hdfs dfs -mkdir /user
# /opt/soft/hadoop/bin/hdfs dfs -mkdir /user/dirac
#
# echo "Start YARN"
# /opt/soft/hadoop/sbin/start-yarn.sh



export PATH=$PATH:/usr/lib/jvm/java-8-openjdk/bin
export JAVA_HOME=/usr/lib/jvm/java-8-openjdk
export HADOOP_CLASSPATH=$(/opt/soft/hadoop/bin/hadoop classpath)
export PATH=$PATH:/opt/soft/hadoop/bin

echo
echo "Creating user directories"
hdfs dfs -mkdir /user
hdfs dfs -mkdir /user/dirac
hdfs dfs -put /opt/dirac/datasets /user/dirac/
hdfs dfs -put titanic.csv /user/dirac/
hdfs dfs -ls /user/dirac/

echo
echo "Show Java Processes"
jps

echo
echo "Get HDFS info"
hdfs dfsadmin -printTopology
hdfs dfsadmin -report

echo
echo "Active network services"
netstat -atup
