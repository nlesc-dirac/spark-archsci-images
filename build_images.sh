#!/bin/sh

# https://docs.docker.com/engine/reference/commandline/build/#squash-an-images-layers-squash-experimental-only

#docker image prune -a

version=$(date +"%d%m%Y")

#export DOCKER_ID_USER=$(whoami)
#docker login

#for folder in $(ls -d */)
for folder in \
	spark-hadoop-base \
	hadoop \
	spark-master \
	spark-worker \
	spark-master-dirac \
	spark-worker-dirac
do
    name=${folder%%/}
    echo
    echo "==========================================="
    echo "Building: "$name
    echo "==========================================="


    cmd="docker build -t fdiblen/$name:$version -t fdiblen/$name:latest ./$name"
    echo
    echo "Running:"
    echo $cmd
    echo
    echo
    eval $cmd

    if [ $? -eq 0 ]
    then
      echo
      echo "Successfully built the $name image!"
    else
      echo
      echo "Could not built the $name image" >&2
      exit $?
    fi


    ### https://docs.docker.com/docker-cloud/builds/push-images/
    #docker push fdiblen/$name:$version
    #docker push fdiblen/$name:latest
done
