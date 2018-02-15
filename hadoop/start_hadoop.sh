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


echo
echo "Creating user directories"
/opt/soft/hadoop/bin/hdfs dfs -mkdir /user
/opt/soft/hadoop/bin/hdfs dfs -mkdir /user/dirac
/opt/soft/hadoop/bin/hdfs dfs -put /opt/data /user/dirac/
/opt/soft/hadoop/bin/hdfs dfs -put titanic.csv /user/dirac/

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
