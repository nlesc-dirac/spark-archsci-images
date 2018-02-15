FROM fdiblen/spark-master-dirac:latest
MAINTAINER fdiblen

USER root
ENV DISPLAY :0

RUN rm -f /etc/supervisor/conf.d/supervisord.conf
COPY conf/supervisord.conf  /etc/supervisor/conf.d/supervisord.conf

COPY start_spark_worker.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/*

#RUN mkdir -p /opt/soft/spark/logs

#RUN mv /opt/soft/spark/conf/spark-env.sh /opt/soft/spark/conf/spark-env.sh.or
#COPY spark-env.sh /opt/soft/spark/conf/

COPY spark_worker_service.sh /tmp/
RUN cat /tmp/spark_worker_service.sh >> /etc/supervisor/conf.d/supervisord.conf && \
    rm -f /tmp/spark_worker_service.sh

RUN sudo pacman --noconfirm -Scc && \
    sudo rm -rf /home/archsci/temp && \
    sudo rm -rf /home/archsci/.cache && \
    sudo rm -rf /var/cache/pacman/pkg/* \
    sudo rm -rf yaourt-tmp-*


ENV JAVA_HOME=/usr/lib/jvm/default
ENV SPARK_VERSION=2.2.1
ENV HADOOP_VERSION=2.7.5

ENV SPARK_HOME=/opt/soft/spark
ENV HADOOP_HOME=/opt/soft/hadoop

ENV HADOOP_NATIVE=/opt/soft/hadoop/lib/native/
ENV JAVA_LIBRARY_PATH=$JAVA_LIBRARY_PATH:$HADOOP_NATIVE
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HADOOP_NATIVE:/opt/soft/casacore/lib
ENV SPARK_YARN_USER_ENV="JAVA_LIBRARY_PATH=$JAVA_LIBRARY_PATH,LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

ENV HADOOP_COMMON_LIB_NATIVE_DIR=/opt/soft/hadoop/lib/native
ENV HADOOP_OPTS="$HADOOP_OPTS -Djava.library.path=$HADOOP_HOME/lib/native"


ENV SPARK_MASTER "spark://sparkmaster:7077"
ENV SPARK_WORKER_WEBUI_PORT 8081
ENV SPARK_WORKER_LOGFILE=/opt/soft/spark/logs/spark_worker.log


EXPOSE 22 8081

CMD ["sudo","/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]
