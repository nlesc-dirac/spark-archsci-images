# About
Docker images for Apache Swarm cluster with HDFS. The images can be used to have a cluster for general use.
There are also [images](https://github.com/nlesc-dirac/spark-docker-swarm) for DIRAC project to run SageCal


# Images

## General 
| Name    | info   |
| ------- | --------- |
| spark-hadoop-base  | base image for Spark and Hadoop |
| spark-master   |  Spark master node image |
| spark-worker    | Spark worker (slave) node image |
| hadoop | Hadoop (HDFS)  image |


## Images for DIRAC project
| Name    | info   |
| ------- | --------- |
| spark-master-dirac | Spark master node with Casacore, BLAS, SageCal|
| spark-worker-dirac | Spark worker node with Casacore, BLAS, SageCal |

