#!/bin/sh

java -Djava.library.path=$(pwd) -jar Driver.jar -m /opt/dirac/datasets/sm.ms

java -Djava.library.path=$(pwd) \
	-jar $(pwd)/Driver.jar \
	-m /opt/dirac/datasets/sm.ms

