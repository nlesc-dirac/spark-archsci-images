FROM fdiblen/spark-master:latest
MAINTAINER fdiblen

USER root
ENV DISPLAY :0

#RUN rm -f /etc/supervisor/conf.d/supervisord.conf
#COPY conf/supervisord.conf  /etc/supervisor/conf.d/supervisord.conf

ADD scripts/install_dependencies.sh /tmp/
RUN chmod +x /tmp/install_dependencies.sh && \
    bash /tmp/install_dependencies.sh && \
    rm -f /tmp/install_dependencies.sh

ADD scripts/compile_casa.sh /tmp/
RUN chmod +x /tmp/compile_casa.sh && \
    bash /tmp/compile_casa.sh && \
    rm -f /tmp/compile_casa.sh

COPY dirac /opt/dirac

ADD scripts/compile_sagecal.sh /tmp/
RUN chmod +x /tmp/compile_sagecal.sh && \
    bash /tmp/compile_sagecal.sh && \
    rm -f /tmp/compile_sagecal.sh


#COPY start_spark_master.sh /usr/local/bin/
#RUN chmod +x /usr/local/bin/*

#RUN mkdir /opt/soft/spark/logs

#RUN mv /opt/soft/spark/conf/spark-env.sh /opt/soft/spark/conf/spark-env.sh.orig
#COPY spark-env.sh /opt/soft/spark/conf/

#COPY spark_master_service.sh /tmp/
#RUN cat /tmp/spark_master_service.sh >> /etc/supervisor/conf.d/supervisord.conf && \
#    rm -f /tmp/spark_master_service.sh

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

ENV SPARK_MASTER_PORT 7077
ENV SPARK_MASTER_WEBUI_PORT 8080
ENV SPARK_MASTER_LOGFILE=/opt/soft/spark/logs/spark_master.log


EXPOSE 22 8080 7077 6066

CMD ["sudo","/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]
