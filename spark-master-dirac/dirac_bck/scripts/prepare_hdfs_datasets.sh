#!/bin/bash

#hdfs --help
#hdfs dfs --help

hdfs dfsadmin –report
hdfs dfs -mkdir /user
hdfs dfs –chmod –R 755 /user

hdfs dfs -mkdir /user/dirac
hdfs dfs –chown –R dirac:dirac

hdfs dfs -mkdir -p /user/dirac/datasets
hdfs dfs -put /opt/dirac/datasets /user/dirac/
hdfs dfs -ls /user/dirac/datasets

#hdfs dfs –rm /user/test/test.txt
