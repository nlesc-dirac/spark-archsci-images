hdfs dfs -mkdir -p /user/dirac
hdfs dfs -ls /user/dirac
hdfs dfs -put /opt/dirac/datasets/sm.ms_sf.tar /user/dirac
hdfs dfs -put /opt/dirac/datasets/sm.ms /user/dirac
hdfs dfs -ls /user/dirac
